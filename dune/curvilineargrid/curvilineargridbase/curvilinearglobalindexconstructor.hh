#ifndef DUNE_CURVILINEARGLOBALINDEXCONSTRUCTOR_HH
#define DUNE_CURVILINEARGLOBALINDEXCONSTRUCTOR_HH



#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/loggingtimer.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>


namespace Dune
{


/** Generates Global Indices for Edges, Faces and Elements
*
* Algorithm:
* 1) Communicate neighbour ranks associated with each process boundary corner
* 2) Compute (provisional) neighbour ranks of PB edges and faces by intersection of ranks of associated PB corners
* 2.1) Sometimes, an entity does not exist on a neighbouring process, even though all associated PB corners are present
* This only happens if the (provisional) number of neighbors is larger than 1 (complicated PB entity),
*      because each PB entity must have at least 1 neighbour.
* 2.2) For each complicated PB entity (edge/face), communicate EdgeKeys and FaceKeys to all provisional neighbours
* 2.3) For each received key, reply to sender if such entity exists on this process or not
* 2.4) Remove neighbour ranks mapping to non-existing entities
*
* 3) Find ownership of each edge and face. A shared entity is owned by the process with lowest rank
* 4) Communicate number of edges and faces owned by each process to all
* 5) Locally enumerate all edges, faces and elements owned by this process. That is, to assign them a global index
* 5.1) Global index for edges starts at nVertexTotal+nEdgesOwnedBeforeMe.
* 5.2) Global index for faces starts at nVertexTotal+nEdgeTotal+nFacesOwnedBeforeMe.
* 5.3) Global index for elements starts at nVertexTotal+nEdgesTotal+nFacesTotal+nElementsOwnedBeforeMe. Note that each process owns all its elements since they are not shared.
* 6) Communicate missing edge and face globalIndices
* 6.1) By analysing entity neighbours, each process can compute how many how many global indices it needs to send and to receive to each other process
* 6.2) Each process sends to each neighbour the shared entity global indices enumerated by this process and receives those enumerated by the neighbour process
* 7) Fill in Global2Local maps. They are required for user functionality and for construction of GhostElements
*
* [TODO] Communication of corner neighbour ranks via allgather very inefficient. Try to find better algorithm
* [TODO] MinRank-Ownership paradigm non-uniform. If ever becomes bottleneck, replace by XORRank-Ownership
*
* */

template <class GridBase>
class CurvilinearGlobalIndexConstructor
{
    /* public types */
	typedef typename GridBase::ctype   ctype;
	static const int dimension = GridBase::dimension;
	typedef typename Dune::ReferenceElements<ctype, dimension>  ReferenceElements;

	typedef          GridBase                                   GridBaseType;
	typedef typename GridBase::GridStorageType                  GridStorageType;
    typedef typename GridBase::LoggingMessage                   LoggingMessage;
    typedef typename Dune::LoggingTimer<LoggingMessage>         LoggingTimer;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    typedef typename GridStorageType::Vertex                    Vertex;
    typedef typename GridStorageType::VertexStorage             VertexStorage;
    typedef typename GridStorageType::EdgeStorage               EdgeStorage;
    typedef typename GridStorageType::FaceStorage               FaceStorage;
    typedef typename GridStorageType::EntityStorage             EntityStorage;

    typedef typename GridStorageType::EdgeKey                   EdgeKey;
    typedef typename GridStorageType::FaceKey                   FaceKey;

    typedef std::map<EdgeKey, LocalIndexType>                   EdgeKey2EdgeIndexMap;
    typedef std::map<FaceKey, LocalIndexType>                   FaceKey2FaceIndexMap;
    typedef typename EdgeKey2EdgeIndexMap::iterator             EdgeMapIterator;
    typedef typename FaceKey2FaceIndexMap::iterator             FaceMapIterator;

    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;

    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;




    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

public:


    CurvilinearGlobalIndexConstructor(
    	GridStorageType & gridstorage,
    	GridBaseType & gridbase,
    	MPIHelper &mpihelper,

    	EdgeKey2EdgeIndexMap & edgeKey2LocalIndexMap,
    	FaceKey2FaceIndexMap & internalFaceKey2LocalIndexMap,
    	FaceKey2FaceIndexMap & domainBoundaryFaceKey2LocalIndexMap,
    	FaceKey2FaceIndexMap & processBoundaryFaceKey2LocalIndexMap) :
    	   	edgeKey2LocalIndexMap_(edgeKey2LocalIndexMap),
    	   	internalFaceKey2LocalIndexMap_(internalFaceKey2LocalIndexMap),
    	   	domainBoundaryFaceKey2LocalIndexMap_(domainBoundaryFaceKey2LocalIndexMap),
    	   	processBoundaryFaceKey2LocalIndexMap_(processBoundaryFaceKey2LocalIndexMap),

    	    gridstorage_(gridstorage),
    	    gridbase_(gridbase),
    	    mpihelper_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        // 1) Communicate process ranks associated with each process boundary corner
        globalCommunicateCornerNeighborRank();

    }


    void generateEdgeGlobalIndex()
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // 2) Compute neighbour ranks of PB edges and faces by intersection of ranks of associated PB corners
        //    Then eliminate non-existing entities generated this way
        // *************************************************************************
        globalComputeEdgeNeighborRank();


        // 3) Get edges and faces on this process that are not owned by this process
        // *************************************************************************
        LocalIndexSet edgeNonOwned;  // LocalIndex of edge not owned by this process

        for (Local2LocalIterator edgeIter  = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].begin();
        		                 edgeIter != gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end(); edgeIter++ )
        {
        	LocalIndexType thisEdgeLocalIndex = (*edgeIter).first;
        	LocalIndexType thisEdgePBIndex    = (*edgeIter).second;

            int edgeOwnerCandidateRank = gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisEdgePBIndex][0];
            if (edgeOwnerCandidateRank < rank_) { edgeNonOwned.insert(thisEdgeLocalIndex); }
        }

        int nEdgeOwned = gridstorage_.edge_.size() - edgeNonOwned.size();


        // 4) Communicate number of edges and faces owned by each process to all
        // *************************************************************************
        std::vector<int> edgesOnProcess(size_, 0);      // owned edges [rank]

        collective_comm.allgather (&nEdgeOwned, 1, reinterpret_cast<int*> (edgesOnProcess.data()));

        int edgesBeforeMe = 0;      // Sum(edgesOwned : rank < thisRank)

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	gridstorage_.nEntityTotal_[EDGE_CODIM] += edgesOnProcess[iProc];
            if (iProc < rank_)  { edgesBeforeMe += edgesOnProcess[iProc]; }
        }


        // 5) Enumerate all edges, faces and elements that you own
        // This process owns this edge if it is in the edgemap and if it is not in the non-owned edge map
        // *************************************************************************
        GlobalIndexType iEdgeGlobalId = edgesBeforeMe;
        for (LocalIndexType iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)
        {
            if (edgeNonOwned.find(iEdge) == edgeNonOwned.end())  { gridstorage_.edge_[iEdge].globalIndex = iEdgeGlobalId++; }
            else { LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: do not own edge localIndex=" + std::to_string(iEdge) + " of type=" + Dune::PartitionName(gridstorage_.edge_[iEdge].ptype)); }
        }


        // 6) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
        // *************************************************************************
        globalDistributeMissingEdgeGlobalIndex();


    }


    void generateFaceGlobalIndex()
    {
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();


        // 2) Compute neighbour ranks of PB edges and faces by intersection of ranks of associated PB corners
        //    Then eliminate non-existing entities generated this way
        // *************************************************************************
        globalComputeFaceNeighborRank();


        // 3) Get edges and faces on this process that are not owned by this process
        // *************************************************************************
        LocalIndexSet faceNonOwned;  // LocalIndex of face not owned by this process

        for (Local2LocalIterator faceIter = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].begin();
        		                 faceIter != gridstorage_.processBoundaryIndexMap_[FACE_CODIM].end(); faceIter++ )
        {
        	LocalIndexType thisFaceLocalIndex = (*faceIter).first;
        	LocalIndexType thisFacePBIndex = (*faceIter).second;
            int faceOwnerCandidateRank = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFacePBIndex][0];
            if (faceOwnerCandidateRank < rank_) { faceNonOwned.insert(thisFaceLocalIndex); }
        }

        int nFaceOwned = gridstorage_.face_.size() - faceNonOwned.size();


        // 4) Communicate number of edges and faces owned by each process to all
        // *************************************************************************
        std::vector<int> facesOnProcess(size_);      // owned faces [rank]

        collective_comm.allgather (&nFaceOwned, 1, reinterpret_cast<int*> (facesOnProcess.data()));

        int facesBeforeMe = 0;      // Sum(facesOwned : rank < thisRank)
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	gridstorage_.nEntityTotal_[FACE_CODIM] += facesOnProcess[iProc];
            if (iProc < rank_)  { facesBeforeMe += facesOnProcess[iProc]; }
        }


        // 5) Enumerate all edges, faces and elements that you own
        // *************************************************************************
        GlobalIndexType iFaceGlobalId = facesBeforeMe;

        // Faces that are not shared with other processes are automatically owned by this process
        for (FaceMapIterator faceIter = internalFaceKey2LocalIndexMap_.begin(); faceIter != internalFaceKey2LocalIndexMap_.end(); faceIter++)              { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.face_[localIndex].globalIndex = iFaceGlobalId++; }
        for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)  { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.face_[localIndex].globalIndex = iFaceGlobalId++; }

        // This process owns this face if it is in the processInternalMap and if it is not in the non-owned face map
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	LocalIndexType thisFaceLocalIndex = (*faceIter).second;
            if (faceNonOwned.find(thisFaceLocalIndex) == faceNonOwned.end())  {
            	gridstorage_.face_[thisFaceLocalIndex].globalIndex = iFaceGlobalId++;
            }
        }

        // 6) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
        // *************************************************************************
        globalDistributeMissingFaceGlobalIndex();
    }


    void generateElementGlobalIndex()
    {
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        int elementsOwned = gridstorage_.element_.size();

        // 4) Communicate number of edges and faces owned by each process to all
        // *************************************************************************
        std::vector<int> elementsOnProcess(size_);  // owned elements [rank]

        collective_comm.allgather (&elementsOwned, 1, reinterpret_cast<int*> (elementsOnProcess.data()));

        int elementsBeforeMe = 0;   // Sum(elementsOwned : rank < thisRank)

        for (int iProc = 0; iProc < size_; iProc++)
        {
            if (iProc < rank_)  { elementsBeforeMe += elementsOnProcess[iProc]; }
        }


        // 5) Enumerate all edges, faces and elements that you own
        // Enumerating elements is simply shifting the local index, since all elements on this process are owned by it
        // *************************************************************************
        for (int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = elementsBeforeMe + i; }


        // 6) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
        // *************************************************************************

        // There are no missing element global indices, because each process owns all of its interior elements
    }




protected:

    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGlobalIndices()
     * ***************************************************************************/

    /** Communicate the process ranks of neighbor processes for all process boundary vertices
     *
     * Algorithm:
     * 1) collective_comm.max() - find the maximal number of process boundary corners per process
     * 2) Loop over maximal number of process boundary corners per process
     * 2.1) collective_comm.allgather() - communicate a global index of your process boundary corner to all other processes
     * 2.2) If all process boundary corners of this process have already been communicated, create and communicate fake indices (negative)
     *      This is necessary to keep the protocol going until the last process communicates all its corners
     * 2.3) From received vertices, select ones that are on this process, and mark the sender as the neighbor
     *
     * [TODO] Algorithm possibly inefficient. No ideas how to improve at the moment
     * * Every process boundary corners is communicated to all processes, but can be used only by few
     * * All processes have to wait until the process with the largest number of corners finishes communicating, since they could receive sth from it
     *
     * [TODO] Currently uses negative LocalIndex values to pass fake vertex info. If we want a uint impl, need to have different def. for fake vertex
     *
     *
     * */

    void globalCommunicateCornerNeighborRank ()
    {
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started communicating corner process boundary neighbors");

        // 1) collective_comm.max() - find the maximal number of process boundary corners per process
        // ********************************************************

        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // Reserve memory for saving ranks associated to process boundary corners
        LocalIndexType thisProcessBoundarySize = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].size();

        LocalIndexType maxProcessBoundarySize = collective_comm.max(thisProcessBoundarySize);

        // 2) collective_comm.allgather() - communicate global index of your process boundary corner to all other processes
        // Repeat this process until every process has communicated all its corners
        // If you run out of corners, communicate fake corners
        // ********************************************************
        Local2LocalIterator procCornerIter = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].begin();

        for (LocalIndexType iCorner = 0; iCorner < maxProcessBoundarySize; iCorner++)
        {
        	Dune::LoggingMessage::writePatience("Determining boundary corner process neighbours...", iCorner, maxProcessBoundarySize);

            // If all process boundary corners have been sent start sending fake corners
        	GlobalIndexType thisCornerGlobalIndex = -1;
        	if (iCorner < thisProcessBoundarySize)
        	{
        		LocalIndexType thisCornerLocalIndex = (*(procCornerIter++)).first;
        		thisCornerGlobalIndex = gridstorage_.point_[thisCornerLocalIndex].globalIndex;
        	}

            // Communicate global indices
            std::vector<LocalIndexType> processBoundaryGlobalIndexSet (size_);
            collective_comm.allgather(&thisCornerGlobalIndex, 1, reinterpret_cast<LocalIndexType*> (processBoundaryGlobalIndexSet.data()) );

            // Loop over corners sent by other processes. If this corner present, note its sender rank
            for (int iProc = 0; iProc < size_; iProc++)
            {
                // Only consider non-fake corners sent by other processes
                if ((iProc != rank_) && (processBoundaryGlobalIndexSet[iProc] >= 0))
                {
                    // Attempt to find this corner global id among process boundary corners of this process
                	GlobalIndexType thisCornerGlobalIndex = processBoundaryGlobalIndexSet[iProc];
                	Global2LocalIterator tmpIter = gridstorage_.entityIndexMap_[VERTEX_CODIM].find(thisCornerGlobalIndex);

                    // If this corner is present, note its sender process
                    if (tmpIter != gridstorage_.entityIndexMap_[VERTEX_CODIM].end()) {
                    	LocalIndexType thisCornerLocalIndex = (*tmpIter).second;
                    	LocalIndexType thisCornerLocalPBIndex = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndex];
                    	gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisCornerLocalPBIndex].push_back(iProc);
                    }
                }
            }
        }

        // 3) Sort all neighbor rank sets, to accelerate set intersection algorithm in future
        // ********************************************************
        for (int i = 0; i < gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM].size(); i++)
        {
            std::sort(gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][i].begin(), gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][i].end());
        }


        // Testing output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: -- Process boundary corner";
        for (Local2LocalIterator cornerIter = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].begin(); cornerIter != gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].end(); cornerIter++)
        {
        	log_stream << " GlobalIndex=" << (*cornerIter).first;
        	log_stream << " has Neighbors=(" << Dune::VectorHelper::vector2string(gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][(*cornerIter).second]) << ")";
        }
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished corner process boundary neighbors");
    }


    /** Compute the process ranks of neighbor processes for all process boundary edges
     *
     * Algorithm:
     * 1) Loop over all edges in the edge map
     * 1.1) For each corner in the EdgeKey get associated neighbor ranks from provided vertex neighbor ranks
     * 1.2) Perform intersection on the two sets
     * 1.3) Following edge map write that intersection to the output array
     *
     * */
    void globalComputeEdgeNeighborRank()
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing edge process boundary neighbors");


        // For each process stores the set of edge indices local to this process, which are shared between more than two processes in total
    	typedef std::pair<LocalIndexType, EdgeKey> TmpEdgeData;

        std::vector<std::vector<TmpEdgeData > > neighborProcessComplicatedEdgePBLocalIndex(size_);

    	// 1) Compute neighbor ranks for each process boundary edge by intersecting neighbor ranks of its corners
    	// *************************************************************************************************************
        int edgeCount = 0;
        for (Local2LocalIterator edgeIter = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].begin();
        		                 edgeIter != gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end(); edgeIter++ )
        {
        	Dune::LoggingMessage::writePatience("Determining boundary edge process neighbours...", edgeCount++, gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size());

            LocalIndexType thisEdgeLocalIndex = (*edgeIter).first;
            LocalIndexType thisPBEdgeLocalIndex = (*edgeIter).second;

            // Get corners of the edge
            std::vector<LocalIndexType> thisCornerLocalIndices = gridbase_.entityCornerLocalIndex(EDGE_CODIM, thisEdgeLocalIndex);

            EdgeKey thisEdgeKey;
            thisEdgeKey.node0 = gridstorage_.point_[thisCornerLocalIndices[0]].globalIndex;
            thisEdgeKey.node1 = gridstorage_.point_[thisCornerLocalIndices[1]].globalIndex;
            thisEdgeKey.sort();

            // Get neighbor processes associated with each corner
            LocalIndexType thisVertexPBIndex0 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndices[0]];
            LocalIndexType thisVertexPBIndex1 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndices[1]];

            std::vector<int> corner0neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex0];
            std::vector<int> corner1neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex1];

            // Find neighbors common to both edge corners
            std::vector<int> edgeneighborset = Dune::VectorHelper::sortedSetIntersection(corner0neighborset, corner1neighborset);

            // Debug info
            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "localEdgeIndex=" << thisEdgeLocalIndex;
            log_stream << " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ")";
            //log_stream << "Neighbors[0]=(" << Dune::VectorHelper::vector2string(corner0neighborset) << ")";
            //log_stream << " Neighbors[1]=(" << Dune::VectorHelper::vector2string(corner1neighborset) << ")";
            log_stream << " Intersection=" << Dune::VectorHelper::vector2string(edgeneighborset);
            LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());


            int nEdgeNeighbor = edgeneighborset.size();
            if (nEdgeNeighbor < 1) {
            	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found no neighbor processes to an edge ");
            	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found no neighbor processes to an edge ");
            }
            else if (nEdgeNeighbor > 1)
            {
            	// Add a complicated edge for further verification
            	// Store only after verification
                TmpEdgeData thisPBEdgeData(thisPBEdgeLocalIndex, thisEdgeKey);
                for (int iEdge = 0; iEdge < edgeneighborset.size(); iEdge++)  {
                	neighborProcessComplicatedEdgePBLocalIndex[edgeneighborset[iEdge]].push_back(thisPBEdgeData);
                }
            } else
            {
                // Store the edge neighbor rank set
                edgeneighborset.swap(gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisPBEdgeLocalIndex]);
            }
        }


        // 2) Communicate to each process the number of complicated edges shared with it
        // *************************************************************************************************************
        std::vector<int> processNComplicatedEdgeRequested(size_);
        std::vector<int> processNComplicatedEdgeToSend(size_);
        for (int iProc = 0; iProc < size_; iProc++)  { processNComplicatedEdgeRequested[iProc] = neighborProcessComplicatedEdgePBLocalIndex[iProc].size(); }
        MPI_Alltoall(processNComplicatedEdgeRequested.data(), 1, MPI_INT, reinterpret_cast<int*>(processNComplicatedEdgeToSend.data()), 1, MPI_INT, comm);

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Total complicated edges per process =(" + Dune::VectorHelper::vector2string(processNComplicatedEdgeRequested) + ")");


        // 3) Communicate to each process the shared complicated edge EdgeKeys
        // *************************************************************************************************************
        int thisCommSize = 0;
        std::vector<int> processEdgeKeyRequested, sdispls;
        std::vector<int> processEdgeKeyToSend, rdispls;

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Sending complicated boundary edge keys to neighbour processes...", iProc, size_);

        	for (int iEdge = 0; iEdge < processNComplicatedEdgeRequested[iProc]; iEdge++)
        	{
        		EdgeKey thisEdgeKey = neighborProcessComplicatedEdgePBLocalIndex[iProc][iEdge].second;
        		processEdgeKeyRequested.push_back(thisEdgeKey.node0);
        		processEdgeKeyRequested.push_back(thisEdgeKey.node1);
        	}

        	processNComplicatedEdgeRequested[iProc] *= 2;
        	processNComplicatedEdgeToSend[iProc]    *= 2;
        	thisCommSize += processNComplicatedEdgeToSend[iProc];

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedEdgeRequested[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedEdgeToSend[iProc-1] );
        }

        processEdgeKeyToSend.resize(thisCommSize);
        MPI_Alltoallv (processEdgeKeyRequested.data(), processNComplicatedEdgeRequested.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processEdgeKeyToSend.data()), processNComplicatedEdgeToSend.data(), rdispls.data(), MPI_INT, comm );
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated complicated edge EdgeKeys");


        // 4) Communicate to each process whether requested edges exist on this process
        // *************************************************************************************************************

        // Note: now we communicate 1 int for every edge key requested, so send and recv switch places and are divided by 2

        thisCommSize = 0;
        std::vector<int> processEdgeExistToSend;      rdispls.clear();
        std::vector<int> processEdgeExistRequested;   sdispls.clear();

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Determining whether the received complicated edge candidates exist...", iProc, size_);

        	processNComplicatedEdgeRequested[iProc] /= 2;
        	processNComplicatedEdgeToSend[iProc] /= 2;
        	thisCommSize += processNComplicatedEdgeRequested[iProc];

        	for (int iEdge = 0; iEdge < processNComplicatedEdgeToSend[iProc]; iEdge++)
        	{
        		EdgeKey thisEdgeKey;
        		thisEdgeKey.node0 = processEdgeKeyToSend[iData++];
        		thisEdgeKey.node1 = processEdgeKeyToSend[iData++];

        		bool isReal = (edgeKey2LocalIndexMap_.find(thisEdgeKey) != edgeKey2LocalIndexMap_.end());
        		processEdgeExistToSend.push_back( isReal ? 1 : 0 );
        	}

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedEdgeToSend[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedEdgeRequested[iProc-1] );
        }

        processEdgeExistRequested.resize(thisCommSize, 0);
        MPI_Alltoallv (processEdgeExistToSend.data(), processNComplicatedEdgeToSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processEdgeExistRequested.data()), processNComplicatedEdgeRequested.data(), rdispls.data(), MPI_INT, comm );
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated if requested EdgeKeys correspond to real edges");


        // 5) Fill in correct neighbors for complicated edges
        // *************************************************************************************************************
        iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Storing correct complicated edge neighbours...", iProc, size_);

        	for (int iEdge = 0; iEdge < processNComplicatedEdgeRequested[iProc]; iEdge++)
        	{
        		bool isReal = (processEdgeExistRequested[iData++] == 1);
        		LocalIndexType thisEdgePBLocalIndex = neighborProcessComplicatedEdgePBLocalIndex[iProc][iEdge].first;
        		if (isReal)  { gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisEdgePBLocalIndex].push_back(iProc); }

        		std::stringstream log_stream;
        		log_stream << " complicated edge PBIndex=" << thisEdgePBLocalIndex << " marked as real=" << isReal << " by process " << iProc;
        		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
        	}
        }
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing edge process boundary neighbors");


        // 6) Sort all edge neighbor rank sets
        // *************************************************************************************************************
        for (int iEdge = 0; iEdge < gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].size(); iEdge++)
        {
        	Dune::LoggingMessage::writePatience("Sorting edge neighbour rank sets...", iEdge, gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].size());
        	std::sort(gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][iEdge].begin(), gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][iEdge].end());
        }


    }


    /** Compute the process ranks of neighbor processes for all process boundary faces
     *
     *  Algorithm:
     *  1) Loop over all process boundary faces in the face map
     *  1.1) For each face corner in the FaceKey get associated neighbor ranks from provided vertex neighbor ranks
     *  1.2) Perform intersection on the three sets
     *  1.3) Ideally the intersection should result in one single rank, which is this face's neighbor. Otherwise throw error
     *
     * */
    void globalComputeFaceNeighborRank()
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing face process boundary neighbors");


        // For each process stores the set of face indices local to this process, which are shared between more than two processes in total
    	typedef std::pair<LocalIndexType, FaceKey> TmpFaceData;
        std::vector<std::vector<TmpFaceData > > neighborProcessComplicatedFaceLocalIndex(size_);

        int faceCount = 0;
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++ )
        {
        	Dune::LoggingMessage::writePatience("Determining set of complicated faces...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

            // Get corners of the face
            FaceKey thisFaceKey = (*faceIter).first;
            LocalIndexType thisFaceLocalIndex = (*faceIter).second;

            // Get neighbor processes associated with each corner
            LocalIndexType thisVertexLocalIndex0 = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisFaceKey.node0];
            LocalIndexType thisVertexLocalIndex1 = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisFaceKey.node1];
            LocalIndexType thisVertexLocalIndex2 = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisFaceKey.node2];

            LocalIndexType thisVertexPBIndex0 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisVertexLocalIndex0];
            LocalIndexType thisVertexPBIndex1 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisVertexLocalIndex1];
            LocalIndexType thisVertexPBIndex2 = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisVertexLocalIndex2];

            std::vector<int> corner0neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex0];
            std::vector<int> corner1neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex1];
            std::vector<int> corner2neighborset = gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][thisVertexPBIndex2];

            // Find neighbors common to all 3 face corners. Need to intersect sets twice
            std::vector<int> faceneighborset;
            faceneighborset = Dune::VectorHelper::sortedSetIntersection(corner0neighborset, corner1neighborset);
            faceneighborset = Dune::VectorHelper::sortedSetIntersection(faceneighborset,    corner2neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "localFaceIndex=" << thisFaceLocalIndex;
            log_stream << " FaceKey=(" << thisFaceKey.node0 << "," << thisFaceKey.node1 << "," << thisFaceKey.node2 << ")";
            //log_stream << "Neighbors[0]=(" << Dune::VectorHelper::vector2string(corner0neighborset) << ")";
            //log_stream << " Neighbors[1]=(" << Dune::VectorHelper::vector2string(corner1neighborset) << ")";
            //log_stream << " Neighbors[2]=(" << Dune::VectorHelper::vector2string(corner2neighborset) << ")";
            log_stream << " Intersection=(" << Dune::VectorHelper::vector2string(faceneighborset) << ")";
            LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

            int nFaceNeighbor = faceneighborset.size();

            if (nFaceNeighbor < 1) {
            	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found face with neighbor nProcess=" + std::to_string(nFaceNeighbor));
            	DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected number of neighbor processes to a face");
            }
            else if (nFaceNeighbor > 1)
            {
              	// Add a complicated face for further verification
              	// Store only after verification
            	TmpFaceData thisPBFaceData(thisFaceLocalIndex, thisFaceKey);
                for (int iFace = 0; iFace < nFaceNeighbor; iFace++)
                {
                  	neighborProcessComplicatedFaceLocalIndex[faceneighborset[iFace]].push_back(thisPBFaceData);
                }
            } else
            {
                // Store the face neighbor rank. Face is only allowed to have exactly one neighbor
            	LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][thisFaceLocalIndex];
            	gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].push_back(faceneighborset[0]);
            }
        }


        // 2) Communicate to each process the number of complicated faces shared with it
        // *************************************************************************************************************
        std::vector<int> processNComplicatedFaceRequested(size_);
        std::vector<int> processNComplicatedFaceToSend(size_);
        for (int iProc = 0; iProc < size_; iProc++)  { processNComplicatedFaceRequested[iProc] = neighborProcessComplicatedFaceLocalIndex[iProc].size(); }
        MPI_Alltoall(processNComplicatedFaceRequested.data(), 1, MPI_INT, reinterpret_cast<int*>(processNComplicatedFaceToSend.data()), 1, MPI_INT, comm);
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Complicated faces per process sent=( " + Dune::VectorHelper::vector2string(processNComplicatedFaceRequested) + ") received =(" + Dune::VectorHelper::vector2string(processNComplicatedFaceToSend) + ")");


        // 3) Communicate to each process the shared complicated face FaceKeys
        // *************************************************************************************************************
        int thisCommSize = 0;
        std::vector<int> processFaceKeyRequested, sdispls;
        std::vector<int> processFaceKeyToSend, rdispls;

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Communicating complicated face keys to neighbour processors...", iProc, size_);

        	for (int iFace = 0; iFace < processNComplicatedFaceRequested[iProc]; iFace++)
        	{
        		FaceKey thisFaceKey = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].second;
        		processFaceKeyRequested.push_back(thisFaceKey.node0);
        		processFaceKeyRequested.push_back(thisFaceKey.node1);
        		processFaceKeyRequested.push_back(thisFaceKey.node2);
        	}

        	processNComplicatedFaceRequested[iProc] *= 3;
        	processNComplicatedFaceToSend[iProc]    *= 3;
        	thisCommSize += processNComplicatedFaceToSend[iProc];

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedFaceRequested[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedFaceToSend[iProc-1] );
        }

        processFaceKeyToSend.resize(thisCommSize);
        MPI_Alltoallv (processFaceKeyRequested.data(), processNComplicatedFaceRequested.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processFaceKeyToSend.data()), processNComplicatedFaceToSend.data(), rdispls.data(), MPI_INT, comm );
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated complicated face FaceKeys");

        //std::cout << "process_" << rank_ << "stage 3) sendcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceRequested) << " recvcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceToSend) <<" send=" << Dune::VectorHelper::vector2string(processFaceKeyRequested) << " recv=" << Dune::VectorHelper::vector2string(processFaceKeyToSend) << std::endl;


        // 4) Communicate to each process whether requested faces exist on this process
        // *************************************************************************************************************

        // Note: now we communicate 1 int for every face key requested, so send and recv switch places and are divided by 3

        std::vector<int> processFaceExistToSend;                                          rdispls.clear();
        std::vector<int> processFaceExistRequested(processFaceKeyRequested.size() / 3);   sdispls.clear();

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Determining whether the received complicated face candidates exist...", iProc, size_);

        	processNComplicatedFaceRequested[iProc] /= 3;
        	processNComplicatedFaceToSend[iProc] /= 3;

        	for (int iFace = 0; iFace < processNComplicatedFaceToSend[iProc]; iFace++)
        	{
        		FaceKey thisFaceKey;
        		thisFaceKey.node0 = processFaceKeyToSend[iData++];
        		thisFaceKey.node1 = processFaceKeyToSend[iData++];
        		thisFaceKey.node2 = processFaceKeyToSend[iData++];

        		bool isReal = (processBoundaryFaceKey2LocalIndexMap_.find(thisFaceKey) != processBoundaryFaceKey2LocalIndexMap_.end());
        		processFaceExistToSend.push_back( isReal ? 1 : 0 );
        	}

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + processNComplicatedFaceToSend[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + processNComplicatedFaceRequested[iProc-1] );
        }
        MPI_Alltoallv (processFaceExistToSend.data(), processNComplicatedFaceToSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(processFaceExistRequested.data()), processNComplicatedFaceRequested.data(), rdispls.data(), MPI_INT, comm );
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated correspondence of requested FaceKeys correspond to real faces");

        //std::cout << "process_" << rank_ << "stage 4) sendcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceToSend) << " recvcounts=" << Dune::VectorHelper::vector2string(processNComplicatedFaceRequested) <<" send=" << Dune::VectorHelper::vector2string(processFaceExistToSend) << " recv=" << Dune::VectorHelper::vector2string(processFaceExistRequested) << std::endl;


        // 5) Fill in correct neighbors for complicated faces
        // *************************************************************************************************************
        iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Storing correct complicated face neighbours...", iProc, size_);

        	for (int iFace = 0; iFace < processNComplicatedFaceRequested[iProc]; iFace++)
        	{
        		bool isReal = (processFaceExistRequested[iData++] == 1);
        		LocalIndexType thisFaceLocalIndex = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].first;
        		FaceKey thisFaceKey = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].second;

        		std::stringstream log_stream;
        		log_stream << " complicated face LocalIndex=" << thisFaceLocalIndex;
        		log_stream << " FaceKey=(" << thisFaceKey.node0 << "," << thisFaceKey.node1 << "," << thisFaceKey.node2 << ")";
        		log_stream << " marked as real=" << isReal << " by process " << iProc;
        		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

        		if (isReal)
        		{
        			LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][thisFaceLocalIndex];
        			int nNeighborAlready = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].size();

        			// If the face neighbor has already been assigned, this face has more than 1 real neighbor process, which is impossible
        			if (nNeighborAlready != 0)
        			{
                    	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found face with neighbor more than two even after cross-check");
                    	DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected number of neighbor processes to a face");
        			}

        			gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].push_back(iProc);
        		}
        	}
        }


        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing face process boundary neighbors");
    }


    /** Communicates all process boundary face global Id's to the neighbors if owned
     *
     * Algorithm:
     *
     * 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
     * 1.1) Assemble a global index array to send to each process
     * 1.2) Note how many faces will be received from each process
     * 2) Assemble one big send array from small arrays (FaceKey + globalIndex)
     * 3) MPI_Alltoallv - communicate this array
     * 4) Save global indices for non-owned faces. Find the exact face by using the communicated FaceKey
     *
     * Optimization Proposal:
     * In principle communication of the FaceKey is not necessary. Instead, the natural FaceKey "<" operator
     * can be used to sort all communicated faces, thus allowing the receiving process to "figure out" what are
     * the faces sent to it by sorting its own faces.
     * 1) Sort all process boundary faces wrt FaceKey
     * 2) Fill arrays to send according to this sorted order
     * 3) Make map for each process from rank & received face to local face index accoridng to the sorted FaceKey order
     * 3) Communicate only globalIndices
     * 4) Use constructed map to fill in received global Indices
     *
     * Will decrease the global communication at the expense of increasing local computation time
     *
     * */
    void globalDistributeMissingFaceGlobalIndex()
    {
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing face GlobalIndices");

        typedef std::pair<FaceKey, GlobalIndexType>  FaceInfo;
        std::vector< std::vector< FaceInfo > > facesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_FACEINFO = 4;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which faces are missing");


        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        int faceCount = 0;
        for (FaceMapIterator iter = processBoundaryFaceKey2LocalIndexMap_.begin(); iter != processBoundaryFaceKey2LocalIndexMap_.end(); iter++)
        {
        	Dune::LoggingMessage::writePatience("Determining faces with missing global indices...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

            LocalIndexType localFaceIndex = (*iter).second;
            LocalIndexType localFacePBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][localFaceIndex];
            int neighborRank = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][localFacePBIndex][0];

            // If the neighbor of this face has lower rank, then add it to recv, else to send
            if (neighborRank < rank_)  { recvcounts[neighborRank] += N_INTEGER_FACEINFO; }
            else
            {
            	GlobalIndexType thisGlobalIndex = gridstorage_.face_[localFaceIndex].globalIndex;
                facesToSend[neighborRank].push_back(FaceInfo((*iter).first, thisGlobalIndex ));
            }
        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Sending missing face global indices...", iProc, size_);

            sendcounts[iProc] = facesToSend[iProc].size() * N_INTEGER_FACEINFO;
            totalRecvSize += recvcounts[iProc];
            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + sendcounts[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + recvcounts[iProc-1] );

            for (int j = 0; j < facesToSend[iProc].size(); j++)
            {
                sendbuf.push_back(facesToSend[iProc][j].second);
                sendbuf.push_back(facesToSend[iProc][j].first.node0);
                sendbuf.push_back(facesToSend[iProc][j].first.node1);
                sendbuf.push_back(facesToSend[iProc][j].first.node2);
            }

        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + Dune::VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + Dune::VectorHelper::vector2string(recvcounts) + ")");

        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Marking received face global indices...", iProc, size_);

            int nThisFaceInfo = recvcounts[iProc] / N_INTEGER_FACEINFO;

            for (int iFace = 0; iFace < nThisFaceInfo; iFace++)
            {
                FaceKey thisKey;
                GlobalIndexType thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];
                thisKey.node2 = recvbuf[iData++];

                FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.find(thisKey);

                if (faceIter == processBoundaryFaceKey2LocalIndexMap_.end()) {
                	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated FaceKey does not correspond to any face on this process");
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated FaceKey does not correspond to any face on this process ");
                }
                else
                {
                    LocalIndexType localFaceIndex = (*faceIter).second;
                    gridstorage_.face_[localFaceIndex].globalIndex = thisGlobalId;
                }
            }
        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing face GlobalIndices");
    }


    /** Communicates all process boundary face global indices to the neighbors if owned
     *
     * Algorithm:
     *
     * 1) Loop over all process boundary edges, split edges into ones to be sent and to be received
     * 1.1) If this edge rank lower than all other neighbor ranks, note to send it to all neighbors,
     *      Otherwise note which neighbor to receive it from
     * 1.1) Assemble a global index array to send to each process
     * 1.2) Note how many edges will be received from each process
     * 2) Assemble one big send array from small arrays (EdgeKey + globalIndex)
     * 3) MPI_Alltoallv - communicate this array
     * 4) Save global indices for non-owned edges. Find the exact edge by using the communicated EdgeKey
     *
     * Optimization Proposal:
     * Same as for FaceGlobalIndex algorithm
     *
     * */
    void globalDistributeMissingEdgeGlobalIndex()
    {
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing edge GlobalIndices");

        typedef std::pair<EdgeKey, GlobalIndexType>  EdgeInfo;
        std::vector< std::vector< EdgeInfo > > edgesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_EDGEINFO = 3;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which edges are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        int edgeCount = 0;
        for (Local2LocalIterator iter  = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].begin();
        		                 iter != gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end(); iter++)
        {
        	Dune::LoggingMessage::writePatience("Determining edges with missing global indices...", edgeCount++, gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size());

            LocalIndexType localEdgeIndex = (*iter).first;
            LocalIndexType localEdgePBIndex = (*iter).second;

            // Construct EdgeKey
            std::vector<LocalIndexType> thisCornerLocalIndices =  gridbase_.entityCornerLocalIndex(EDGE_CODIM, localEdgeIndex);
            EdgeKey thisEdgeKey;
            thisEdgeKey.node0 = gridstorage_.point_[thisCornerLocalIndices[0]].globalIndex;
            thisEdgeKey.node1 = gridstorage_.point_[thisCornerLocalIndices[1]].globalIndex;
            thisEdgeKey.sort();


            int candidateOwnerRank = gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex][0];

            //std::cout << "process_" << rank_ <<  " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ") localIndex=" << localEdgeIndex <<  std::endl;

            // If the one of the neighbors of this edge has lower rank, then note one more received edge from that process
            // else note to send it to all other neighbors
            if (candidateOwnerRank < rank_)  { recvcounts[candidateOwnerRank] += N_INTEGER_EDGEINFO; }
            else
            {
            	GlobalIndexType thisGlobalIndex = gridstorage_.edge_[localEdgeIndex].globalIndex;

                EdgeInfo thisEdgeInfo(thisEdgeKey, thisGlobalIndex);

                for (int iNeighbor = 0; iNeighbor < gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex].size(); iNeighbor++)
                {
                    int thisNeighborRank = gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex][iNeighbor];
                    edgesToSend[thisNeighborRank].push_back(thisEdgeInfo);
                };
            }
        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int i = 0; i < size_; i++)
        {
        	Dune::LoggingMessage::writePatience("Sending missing global edge indices to neighbours...", i, size_);

            sendcounts[i] = edgesToSend[i].size() * N_INTEGER_EDGEINFO;
            totalRecvSize += recvcounts[i];
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );

            for (int j = 0; j < edgesToSend[i].size(); j++)
            {
                sendbuf.push_back(edgesToSend[i][j].second);
                sendbuf.push_back(edgesToSend[i][j].first.node0);
                sendbuf.push_back(edgesToSend[i][j].first.node1);
            }

        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + Dune::VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + Dune::VectorHelper::vector2string(recvcounts) + ")");



        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);   // There are 3 integers per sent FaceInfo
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	Dune::LoggingMessage::writePatience("Marking received edge global indices...", iProc, size_);

            int nThisEdgeInfo = recvcounts[iProc] / N_INTEGER_EDGEINFO;

            for (int iEdge = 0; iEdge < nThisEdgeInfo; iEdge++)
            {
                EdgeKey thisKey;
                GlobalIndexType thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];

                EdgeMapIterator edgeIter = edgeKey2LocalIndexMap_.find(thisKey);

                if (edgeIter == edgeKey2LocalIndexMap_.end()) {
                	std::stringstream log_str;
                	log_str << "CurvilinearGridConstructor: Communicated EdgeKey (" << thisKey.node0 << ", " << thisKey.node1 << ") does not correspond to any edge on this process";
                	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, log_str.str());
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated EdgeKey does not correspond to any edge on this process "); }
                else
                {
                    LocalIndexType localEdgeIndex = (*edgeIter).second;
                    gridstorage_.edge_[localEdgeIndex].globalIndex = thisGlobalId;
                }
            }
        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing edge GlobalIndices");
    }


private:

    EdgeKey2EdgeIndexMap & edgeKey2LocalIndexMap_;
    FaceKey2FaceIndexMap & internalFaceKey2LocalIndexMap_;
    FaceKey2FaceIndexMap & domainBoundaryFaceKey2LocalIndexMap_;
    FaceKey2FaceIndexMap & processBoundaryFaceKey2LocalIndexMap_;

    // Curvilinear Grid Storage Class
    GridStorageType & gridstorage_;

    // Reference to Curvilinear Grid Base - necessary for OCTree construction
    GridBaseType & gridbase_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;


};

} // Namespace Dune


#endif // DUNE_CURVILINEARGLOBALINDEXCONSTRUCTOR_HH
