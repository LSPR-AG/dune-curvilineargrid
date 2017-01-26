#ifndef DUNE_CURVILINEARGLOBALINDEXCONSTRUCTOR_HH
#define DUNE_CURVILINEARGLOBALINDEXCONSTRUCTOR_HH



#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/loggingtimer.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/utility/allcommunication.hh>


namespace Dune
{

namespace CurvGrid {

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
    typedef typename GridBase::LoggingTimer                     LoggingTimer;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    typedef typename GridStorageType::GlobalCoordinate                    GlobalCoordinate;
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
    	    mpihelper_(mpihelper),
			allcomm_(mpihelper)
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
            else { LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: do not own edge localIndex=" + std::to_string(iEdge) + " of type=" + Dune::PartitionName(gridstorage_.edge_[iEdge].ptype)); }
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

        for (const auto & faceIndexPair : gridstorage_.processBoundaryIndexMap_[FACE_CODIM])
        {
        	LocalIndexType thisFaceLocalIndex = faceIndexPair.first;
        	LocalIndexType thisFacePBIndex = faceIndexPair.second;
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
        for (const auto & faceIter : internalFaceKey2LocalIndexMap_)					{ LocalIndexType localIndex = faceIter.second;  gridstorage_.face_[localIndex].globalIndex = iFaceGlobalId++; }
        for (const auto & faceIter : domainBoundaryFaceKey2LocalIndexMap_)	{ LocalIndexType localIndex = faceIter.second;  gridstorage_.face_[localIndex].globalIndex = iFaceGlobalId++; }

        // This process owns this face if it is in the processInternalMap and if it is not in the non-owned face map
        for (const auto & faceIter : processBoundaryFaceKey2LocalIndexMap_)
        {
        	LocalIndexType thisFaceLocalIndex = faceIter.second;
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
        for (unsigned int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = elementsBeforeMe + i; }


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
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started communicating corner process boundary neighbors");

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
        	LoggingMessage::writePatience("Determining boundary corner process neighbours...", iCorner, maxProcessBoundarySize);

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
        for (auto && neighborRankVec : gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM])
        {
            std::sort(neighborRankVec.begin(), neighborRankVec.end());
        }


        // Testing output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: -- Process boundary corner";
        for (Local2LocalIterator cornerIter = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].begin(); cornerIter != gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].end(); cornerIter++)
        {
        	log_stream << " GlobalIndex=" << (*cornerIter).first;
        	log_stream << " has Neighbors=(" << VectorHelper::vector2string(gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM][(*cornerIter).second]) << ")";
        }
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished corner process boundary neighbors");
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
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing edge process boundary neighbors");


        // For each process stores the set of edge indices local to this process, which are shared between more than two processes in total
    	typedef std::pair<LocalIndexType, EdgeKey> TmpEdgeData;

        std::vector<std::vector<TmpEdgeData > > neighborProcessComplicatedEdgePBLocalIndex(size_);

    	// 1) Compute neighbor ranks for each process boundary edge by intersecting neighbor ranks of its corners
    	// *************************************************************************************************************
        int edgeCount = 0;
        for (const auto & edgeIndexPair : gridstorage_.processBoundaryIndexMap_[EDGE_CODIM])
        {
        	LoggingMessage::writePatience("Determining boundary edge process neighbours...", edgeCount++, gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size());

            LocalIndexType thisEdgeLocalIndex = edgeIndexPair.first;
            LocalIndexType thisPBEdgeLocalIndex = edgeIndexPair.second;

            // Get corners of the edge
            std::vector<LocalIndexType> thisCornerLocalIndices = gridbase_.entity().cornerLocalIndex(EDGE_CODIM, thisEdgeLocalIndex);

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
            std::vector<int> edgeneighborset = VectorHelper::sortedSetIntersection(corner0neighborset, corner1neighborset);

            // Debug info
            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "localEdgeIndex=" << thisEdgeLocalIndex;
            log_stream << " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ")";
            //log_stream << "Neighbors[0]=(" << VectorHelper::vector2string(corner0neighborset) << ")";
            //log_stream << " Neighbors[1]=(" << VectorHelper::vector2string(corner1neighborset) << ")";
            log_stream << " Intersection=" << VectorHelper::vector2string(edgeneighborset);
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());


            int nEdgeNeighbor = edgeneighborset.size();
            if (nEdgeNeighbor < 1) {
            	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found no neighbor processes to an edge ");
            	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found no neighbor processes to an edge ");
            }
            else if (nEdgeNeighbor > 1)
            {
            	// Add a complicated edge for further verification
            	// Store only after verification
                TmpEdgeData thisPBEdgeData(thisPBEdgeLocalIndex, thisEdgeKey);
                for (unsigned int iEdge = 0; iEdge < edgeneighborset.size(); iEdge++)  {
                	neighborProcessComplicatedEdgePBLocalIndex[edgeneighborset[iEdge]].push_back(thisPBEdgeData);
                }
            } else
            {
                // Store the edge neighbor rank set
                edgeneighborset.swap(gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisPBEdgeLocalIndex]);
            }
        }

        // 3) Communicate to each process the shared complicated edge EdgeKeys
        // *************************************************************************************************************
        int thisCommSize = 0;
        std::vector<int> processEdgeKeyRequested;
        std::vector<int> processNComplicatedEdgeVerticesRequested(size_);

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Sending complicated boundary edge keys to neighbour processes...", iProc, size_);
        	processNComplicatedEdgeVerticesRequested[iProc] = neighborProcessComplicatedEdgePBLocalIndex[iProc].size();

        	for (int iEdge = 0; iEdge < processNComplicatedEdgeVerticesRequested[iProc]; iEdge++)
        	{
        		EdgeKey thisEdgeKey = neighborProcessComplicatedEdgePBLocalIndex[iProc][iEdge].second;
        		processEdgeKeyRequested.push_back(thisEdgeKey.node0);
        		processEdgeKeyRequested.push_back(thisEdgeKey.node1);
        	}
        	processNComplicatedEdgeVerticesRequested[iProc] *= 2;	// Two vertices per edge
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Total complicated edges per process =(" + VectorHelper::vector2string(processNComplicatedEdgeVerticesRequested) + ")");

        std::vector<int> processEdgeKeyToSend;
        std::vector<int> processNComplicatedEdgeVerticesToSend;
        allcomm_.all2all(processEdgeKeyRequested, processNComplicatedEdgeVerticesRequested, processEdgeKeyToSend, processNComplicatedEdgeVerticesToSend);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated complicated edge EdgeKeys");


        // 4) Communicate to each process whether requested edges exist on this process
        // *************************************************************************************************************

        // Note: now we communicate 1 int for every edge key requested, so send and recv switch places and are divided by 2


        std::vector<int> processNComplicatedEdgesToSend(size_);
        std::vector<int> processEdgeExistToSend;

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Determining whether the received complicated edge candidates exist...", iProc, size_);

        	processNComplicatedEdgesToSend[iProc] = processNComplicatedEdgeVerticesToSend[iProc] / 2;

        	for (int iEdge = 0; iEdge < processNComplicatedEdgesToSend[iProc]; iEdge++)
        	{
        		EdgeKey thisEdgeKey;
        		thisEdgeKey.node0 = processEdgeKeyToSend[iData++];
        		thisEdgeKey.node1 = processEdgeKeyToSend[iData++];

        		bool isReal = (edgeKey2LocalIndexMap_.find(thisEdgeKey) != edgeKey2LocalIndexMap_.end());
        		processEdgeExistToSend.push_back( isReal ? 1 : 0 );
        	}
        }

        std::vector<int> processEdgeExistRequested;
        std::vector<int> processNComplicatedEdgesRequested;
        allcomm_.all2all(processEdgeExistToSend, processNComplicatedEdgesToSend, processEdgeExistRequested, processNComplicatedEdgesRequested);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated if requested EdgeKeys correspond to real edges");

        // Check consistency of send-receive
        assert(processNComplicatedEdgeVerticesRequested.size() == processNComplicatedEdgesRequested.size());
        for (unsigned int i = 0; i < processNComplicatedEdgeVerticesRequested.size(); i++) {
        	if (processNComplicatedEdgeVerticesRequested[i] != 2 * processNComplicatedEdgesRequested[i]) {
        		std::cout << rank_ << " bug: " << VectorHelper::vector2string(processNComplicatedEdgeVerticesRequested)
        			<< " ::: " << VectorHelper::vector2string(processNComplicatedEdgesRequested) << std::endl;
        	}

        	assert(processNComplicatedEdgeVerticesRequested[i] == 2 * processNComplicatedEdgesRequested[i]);
        }


        // 5) Fill in correct neighbors for complicated edges
        // *************************************************************************************************************
        iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Storing correct complicated edge neighbours...", iProc, size_);

        	for (int iEdge = 0; iEdge < processNComplicatedEdgesRequested[iProc]; iEdge++)
        	{
        		assert(iData < processEdgeExistRequested.size());
        		assert(iEdge < neighborProcessComplicatedEdgePBLocalIndex[iProc].size());


        		bool isReal = (processEdgeExistRequested[iData++] == 1);
        		LocalIndexType thisEdgePBLocalIndex = neighborProcessComplicatedEdgePBLocalIndex[iProc][iEdge].first;

        		assert(thisEdgePBLocalIndex < gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].size());

        		if (isReal)  { gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][thisEdgePBLocalIndex].push_back(iProc); }

        		std::stringstream log_stream;
        		log_stream << " complicated edge PBIndex=" << thisEdgePBLocalIndex << " marked as real=" << isReal << " by process " << iProc;
        		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
        	}
        }
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing edge process boundary neighbors");


        // 6) Sort all edge neighbor rank sets
        // *************************************************************************************************************
        for (unsigned int iEdge = 0; iEdge < gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].size(); iEdge++)
        {
        	LoggingMessage::writePatience("Sorting edge neighbour rank sets...", iEdge, gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].size());
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
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started computing face process boundary neighbors");


        // For each process stores the set of face indices local to this process, which are shared between more than two processes in total
    	typedef std::pair<LocalIndexType, FaceKey> TmpFaceData;
        std::vector<std::vector<TmpFaceData > > neighborProcessComplicatedFaceLocalIndex(size_);

        int faceCount = 0;
        for (const auto & faceKeyIndexPair : processBoundaryFaceKey2LocalIndexMap_)
        {
        	LoggingMessage::writePatience("Determining set of complicated faces...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

            // Get corners of the face
            FaceKey thisFaceKey = faceKeyIndexPair.first;
            LocalIndexType thisFaceLocalIndex = faceKeyIndexPair.second;

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
            faceneighborset = VectorHelper::sortedSetIntersection(corner0neighborset, corner1neighborset);
            faceneighborset = VectorHelper::sortedSetIntersection(faceneighborset,    corner2neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: --";
            log_stream << "localFaceIndex=" << thisFaceLocalIndex;
            log_stream << " FaceKey=(" << thisFaceKey.node0 << "," << thisFaceKey.node1 << "," << thisFaceKey.node2 << ")";
            //log_stream << "Neighbors[0]=(" << VectorHelper::vector2string(corner0neighborset) << ")";
            //log_stream << " Neighbors[1]=(" << VectorHelper::vector2string(corner1neighborset) << ")";
            //log_stream << " Neighbors[2]=(" << VectorHelper::vector2string(corner2neighborset) << ")";
            log_stream << " Intersection=(" << VectorHelper::vector2string(faceneighborset) << ")";
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

            int nFaceNeighbor = faceneighborset.size();

            if (nFaceNeighbor < 1) {
            	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found face with neighbor nProcess=" + std::to_string(nFaceNeighbor));
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


        // 3) Communicate to each process the shared complicated face FaceKeys
        // *************************************************************************************************************
        std::vector<int> processFaceKeyRequested;
        std::vector<int> processNComplicatedFaceVertexRequested(size_);

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Communicating complicated face keys to neighbour processors...", iProc, size_);

        	processNComplicatedFaceVertexRequested[iProc] = neighborProcessComplicatedFaceLocalIndex[iProc].size();

        	for (int iFace = 0; iFace < processNComplicatedFaceVertexRequested[iProc]; iFace++)
        	{
        		FaceKey thisFaceKey = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].second;
        		processFaceKeyRequested.push_back(thisFaceKey.node0);
        		processFaceKeyRequested.push_back(thisFaceKey.node1);
        		processFaceKeyRequested.push_back(thisFaceKey.node2);
        	}
        	processNComplicatedFaceVertexRequested[iProc] *= 3;
        }

        std::vector<int> processFaceKeyToSend;
        std::vector<int> processNComplicatedFaceToSend;
        allcomm_.all2all(processFaceKeyRequested, processNComplicatedFaceVertexRequested, processFaceKeyToSend, processNComplicatedFaceToSend);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated complicated face FaceKeys");

        //std::cout << "process_" << rank_ << "stage 3) sendcounts=" << VectorHelper::vector2string(processNComplicatedFaceVertexRequested) << " recvcounts=" << VectorHelper::vector2string(processNComplicatedFaceToSend) <<" send=" << VectorHelper::vector2string(processFaceKeyRequested) << " recv=" << VectorHelper::vector2string(processFaceKeyToSend) << std::endl;


        // 4) Communicate to each process whether requested faces exist on this process
        // *************************************************************************************************************

        // Note: now we communicate 1 int for every face key requested, so send and recv switch places and are divided by 3

        std::vector<int> processFaceExistToSend;

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Determining whether the received complicated face candidates exist...", iProc, size_);

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
        }

        std::vector<int> processFaceExistRequested;
        std::vector<int> processNComplicatedFaceRequested;
        allcomm_.all2all(processFaceExistToSend, processNComplicatedFaceToSend, processFaceExistRequested, processNComplicatedFaceRequested);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated correspondence of requested FaceKeys correspond to real faces");

        // Self-test
        for (unsigned int iProc = 0; iProc < size_; iProc++) { assert(processNComplicatedFaceVertexRequested[iProc] == 3 * processNComplicatedFaceRequested[iProc]); }
        assert(processFaceExistRequested.size() == processFaceKeyRequested.size() / 3);

        //std::cout << "process_" << rank_ << "stage 4) sendcounts=" << VectorHelper::vector2string(processNComplicatedFaceToSend) << " recvcounts=" << VectorHelper::vector2string(processNComplicatedFaceVertexRequested) <<" send=" << VectorHelper::vector2string(processFaceExistToSend) << " recv=" << VectorHelper::vector2string(processFaceExistRequested) << std::endl;


        // 5) Fill in correct neighbors for complicated faces
        // *************************************************************************************************************
        iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Storing correct complicated face neighbours...", iProc, size_);

        	for (int iFace = 0; iFace < processNComplicatedFaceRequested[iProc]; iFace++)
        	{
        		bool isReal = (processFaceExistRequested[iData++] == 1);
        		LocalIndexType thisFaceLocalIndex = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].first;
        		FaceKey thisFaceKey = neighborProcessComplicatedFaceLocalIndex[iProc][iFace].second;

        		std::stringstream log_stream;
        		log_stream << " complicated face LocalIndex=" << thisFaceLocalIndex;
        		log_stream << " FaceKey=(" << thisFaceKey.node0 << "," << thisFaceKey.node1 << "," << thisFaceKey.node2 << ")";
        		log_stream << " marked as real=" << isReal << " by process " << iProc;
        		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

        		if (isReal)
        		{
        			LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][thisFaceLocalIndex];
        			int nNeighborAlready = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].size();

        			// If the face neighbor has already been assigned, this face has more than 1 real neighbor process, which is impossible
        			if (nNeighborAlready != 0)
        			{
                    	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Found face with neighbor more than two even after cross-check");
                    	DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected number of neighbor processes to a face");
        			}

        			gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex].push_back(iProc);
        		}
        	}
        }


        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished computing face process boundary neighbors");
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
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing face GlobalIndices");

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();

        typedef std::pair<FaceKey, GlobalIndexType>  FaceInfo;
        std::vector< std::vector< FaceInfo > > facesToSend (size_);

        int N_INTEGER_FACEINFO = 4;
        std::vector<int> sendbuf, sendcounts(size_);
        std::vector<int> recvcountsExpected(size_);

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which faces are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        int faceCount = 0;
        for (const auto & faceIndexPair : processBoundaryFaceKey2LocalIndexMap_)
        {
        	LoggingMessage::writePatience("Determining faces with missing global indices...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

            LocalIndexType localFaceIndex = faceIndexPair.second;
            LocalIndexType localFacePBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM][localFaceIndex];

            assert(gridstorage_.PB2PBNeighborRank_[FACE_CODIM][localFacePBIndex].size() == 1); // Exactly 1 neighbor process for each PB face
            int neighborRank = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][localFacePBIndex][0];

            // If the neighbor of this face has lower rank, then add it to recv, else to send
            if (neighborRank < rank_)  { recvcountsExpected[neighborRank] += N_INTEGER_FACEINFO; }
            else
            {
            	GlobalIndexType thisGlobalIndex = gridstorage_.face_[localFaceIndex].globalIndex;
                facesToSend[neighborRank].push_back(FaceInfo(faceIndexPair.first, thisGlobalIndex ));
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Sending missing face global indices...", iProc, size_);

            sendcounts[iProc] = facesToSend[iProc].size() * N_INTEGER_FACEINFO;

            for (unsigned int j = 0; j < facesToSend[iProc].size(); j++)
            {
                sendbuf.push_back(facesToSend[iProc][j].second);
                sendbuf.push_back(facesToSend[iProc][j].first.node0);
                sendbuf.push_back(facesToSend[iProc][j].first.node1);
                sendbuf.push_back(facesToSend[iProc][j].first.node2);
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + VectorHelper::vector2string(recvcountsExpected) + ")");

        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        std::vector<int> recvbuf;
        std::vector<int> recvcounts;
        allcomm_.all2all(sendbuf, sendcounts, recvbuf, recvcounts);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

        // Self-check
        assert(recvcountsExpected.size() == recvcounts.size());
        for (unsigned int i = 0; i < recvcountsExpected.size(); i++) { assert(recvcountsExpected[i] == recvcounts[i]); }



        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Marking received face global indices...", iProc, size_);

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
                	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated FaceKey does not correspond to any face on this process");
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated FaceKey does not correspond to any face on this process ");
                }
                else
                {
                    LocalIndexType localFaceIndex = (*faceIter).second;
                    gridstorage_.face_[localFaceIndex].globalIndex = thisGlobalId;
                }
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing face GlobalIndices");
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
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started distributing missing edge GlobalIndices");

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();

        typedef std::pair<EdgeKey, GlobalIndexType>  EdgeInfo;
        std::vector< std::vector< EdgeInfo > > edgesToSend (size_);

        int N_INTEGER_EDGEINFO = 3;
        std::vector<int> sendbuf, sendcounts(size_);
        std::vector<int> recvcounts(size_);

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Checking which edges are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        int edgeCount = 0;
        for (const auto & edgeIndexPair : gridstorage_.processBoundaryIndexMap_[EDGE_CODIM])
        {
        	LoggingMessage::writePatience("Determining edges with missing global indices...", edgeCount++, gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size());

            LocalIndexType localEdgeIndex = edgeIndexPair.first;
            LocalIndexType localEdgePBIndex = edgeIndexPair.second;

            // Construct EdgeKey
            std::vector<LocalIndexType> thisCornerLocalIndices =  gridbase_.entity().cornerLocalIndex(EDGE_CODIM, localEdgeIndex);
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

                for (const auto & thisNeighborRank : gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][localEdgePBIndex]) {
                    edgesToSend[thisNeighborRank].push_back(thisEdgeInfo);
                };
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int i = 0; i < size_; i++)
        {
        	LoggingMessage::writePatience("Sending missing global edge indices to neighbours...", i, size_);

            sendcounts[i] = edgesToSend[i].size() * N_INTEGER_EDGEINFO;

            for (unsigned int j = 0; j < edgesToSend[i].size(); j++)
            {
                sendbuf.push_back(edgesToSend[i][j].second);
                sendbuf.push_back(edgesToSend[i][j].first.node0);
                sendbuf.push_back(edgesToSend[i][j].first.node1);
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Sending  sendcounts=(" + VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + VectorHelper::vector2string(recvcounts) + ")");



        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        std::vector<int> recvbuf;   // There are 3 integers per sent FaceInfo
        std::vector<int> recvcounts2;
        allcomm_.all2all(sendbuf, sendcounts, recvbuf, recvcounts2);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Extracting missing indices");

        // Self-check
        assert(recvcounts.size() == recvcounts2.size());
        for (unsigned int i = 0; i < recvcounts.size(); i++) { assert(recvcounts[i] == recvcounts2[i]); }


        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	LoggingMessage::writePatience("Marking received edge global indices...", iProc, size_);

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
                	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, log_str.str());
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated EdgeKey does not correspond to any edge on this process "); }
                else
                {
                    LocalIndexType localEdgeIndex = (*edgeIter).second;
                    gridstorage_.edge_[localEdgeIndex].globalIndex = thisGlobalId;
                }
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished distributing missing edge GlobalIndices");
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

    // MPI Communication wrapper
    AllCommunication allcomm_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;


};

} // Namespace CurvGrid

} // Namespace Dune


#endif // DUNE_CURVILINEARGLOBALINDEXCONSTRUCTOR_HH
