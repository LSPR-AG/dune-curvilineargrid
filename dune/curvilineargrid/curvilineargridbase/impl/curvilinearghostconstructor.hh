/***************************************************************************
                          tetmesh.h  -  description
                             -------------------
    begin                : Mon Dec 15 2003
    copyright            : (C) 2003 by Roman Geus
    email                : roman.geus@psi.ch
    edited by            : Hua Guo, Sep 2010
***************************************************************************/

/***************************************************************************
                          curvilineargridbase.hh
                             -------------------
    begin                : Tue Nov 25 2014
    copyright            : (C) 2014 by Aleksejs Fomins, LSPR AG
    description          : Upgraded the mesh to curvilinear grid
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEARGHOSTCONSTRUCTOR_HH
#define DUNE_CURVILINEARGHOSTCONSTRUCTOR_HH

#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/impl/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearperiodicconstructor.hh>

#include <dune/curvilineargrid/utility/allcommunication.hh>


namespace Dune {

namespace CurvGrid {

// Forward declaration
//template<class ct, int cdim, bool isCached>
//class CurvilinearGridStorage;


template <class GridBase>
class CurvilinearGhostConstructor {
public:

    /* public types */
	static const int dimension = GridBase::dimension;

    typedef typename GridBase::GridStorageType                  GridStorageType;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
    typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    // Containers
    typedef typename GridStorageType::GlobalCoordinate                    GlobalCoordinate;
    typedef typename GridStorageType::VertexStorage             VertexStorage;
    typedef typename GridStorageType::EdgeStorage               EdgeStorage;
    typedef typename GridStorageType::FaceStorage               FaceStorage;
    typedef typename GridStorageType::EntityStorage             EntityStorage;

    // Local and global maps and their iterators
    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    // Periodic constructions
    typedef CurvilinearPeriodicConstructor<GridBase>  PeriodicConstr;
    typedef typename PeriodicConstr::GlobalIndexPeriodicScatterDataMap  GlobalIndexPeriodicScatterDataMap;

    // MPI Wrapper




public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
	CurvilinearGhostConstructor(
    		GridStorageType & gridstorage,
    		MPIHelper &mpihelper,
			const PeriodicConstr * periodicConstr = nullptr
	) :
        gridstorage_(gridstorage),
        mpihelper_(mpihelper),
		allcomm_(mpihelper),
		periodicConstr_(periodicConstr),
		nPeriodicGhostNeighbor_(0),
		nProcessGhost_(0)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "Initialized CurvilinearGhostConstructor");
    }


    /** Communicates the Ghost Elements
     *
     * Prerequisites:
     * * Requires all processBoundaries to know neighboring process rank
     * * Requires existence of face global index
     *
     * Aspects:
     * * Ghost elements can have different interpolation order
     *
     * Algorithm:
     * 1) Compute number of process boundaries shared with each process, as well as and the exact interpolation orders of Ghost Elements
     * 1.1) Note that a ghost element can have more than 1 process boundary face associated to it, so a set of faces needs to be stored
     * 1.2) Note that for this reason it is not possible to know in advance how many ghosts will be received from a neighbor,
     *      so it has to be communicated
     * 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the Ghost Elements it is going to receive
     * 3) MPI_alltoallv - Package and send element globalIndex + elementPhysicalTag + all interpolatory Vertex global indices
     * 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
     * 5) MPI_alltoall - Communicates to each process the number of missing vertices out of the ones it had communicated.
     * 5.1) Then communicate the globalIndices's of all missing vertices
     * 5.2) Then communicate the vertex coordinates corresponding to received global indices
     * 6) Distrubute vertex coordinates and add received coordinates to the mesh
     *
     *
     * TODO: Adapt for the case with no ghost elements whatsoever - then the 2nd neighbor of face should be another periodic face
     * TODO: Only supports tetrahedral ghost elements at the moment
     *
     * */
    void generate()
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Generating Ghost Elements");

        thisProcessNeighborGhostLocalIndex_.resize(size_, std::vector<int>() );
        thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_.resize(size_);
        neighborProcessGhostInterpOrder_.resize(size_, std::vector<int>() );
        neighborProcessNAssociatedFace_.resize(size_, std::vector<int>() );


        // 1) Compute number of process boundaries shared with each process, as well as and the exact interpolation orders of Ghost Elements
        // *************************************************************************************
        ghostComputeLocalNeighborGhostElements();


        // 2) MPI_alltoallv - communicate to each process the number of Ghost elements it is going to receive
        // Then for each element communicate the associated ghost element interpolation order and associated number of process boundary faces
        // *************************************************************************************
        ghostDistributePreliminaries();


        // 3) MPI_alltoallv - Package element globalIndex + elementPhysicalTag + all associated face global indices + all interpVertex globalIds
        // *************************************************************************************
        std::vector<int> packageGhostElementData;
        std::vector<int> packageGhostElementCounts;
        ghostDistributeGhostElements(packageGhostElementData, packageGhostElementCounts);

        thisProcessNeighborGhostLocalIndex_.clear();
        thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_.clear();


        // 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
        // *************************************************************************************
        std::vector<std::vector<int> > missingVertices (size_);
        ghostInsertGhostElements(packageGhostElementData, packageGhostElementCounts, missingVertices);

        neighborProcessGhostInterpOrder_.clear();
        neighborProcessNAssociatedFace_.clear();
        packageGhostElementData.clear();


        // 5) Communicates to each process the number of missing vertices out of the ones it had communicated
        // Then communicate the globalId's of all missing vertices
        // *************************************************************************************
        std::vector<int> packageMissingVertexGlobalIndices;
        std::vector<int> verticesRequestedByThis;
        std::vector<int> verticesToSendByThis;

        ghostCommunicateMissingVertexGlobalIndices(missingVertices, packageMissingVertexGlobalIndices, verticesRequestedByThis, verticesToSendByThis);


        // 6) Distrubute vertex coordinates and add received coordinates to the mesh
        // Note that verticesRequestedByThis and verticesToSendByThis arrays are in the reversed order: this time we send as much as we have received before, and receive as much as we have sent before
        // *************************************************************************************
        ghostCommunicateMissingVertexCoordinates (missingVertices, packageMissingVertexGlobalIndices, verticesToSendByThis, verticesRequestedByThis);

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished Generating Ghost Elements");


        // ConsistencyDump - for debug
		// std::stringstream log_str;
		// log_str << "Consistency Dump nVertex=" << gridstorage_.point_.size() << " with mapsize=" << gridstorage_.entityIndexMap_[VERTEX_CODIM].size() << std::endl;
		// for (int i = 0; i < gridstorage_.point_.size(); i++)  { log_str << "  -- localIndex=" << i << " globalIndex=" << gridstorage_.point_[i].globalIndex << " coord=(" << gridstorage_.point_[i].coord << ")" << std::endl;  }
		//LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_str.str());
    }



protected:


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGhostElements()
     * ***************************************************************************/


    /** For each neighbor process, compute the set of internal elements that will be sent to it as ghost elements
     *
     * Algorithm:  (No communication needed)
     *
     * 1) Loop over all Periodic faces, if they exist
     * 1.1) Verify the face structural type was originally a domain boundary, and set it to periodic boundary
     * 1.2) Check if periodic neighbor is local. If yes, it has already been treated in the periodic constructor and can safely be ignored
     * 1.3) Otherwise, add this face as subentity for interior element to be communicated as ghost
     * 2) Loop over all PB faces
     * 2.1) Verify the face structural type is periodic boundary
     * 2.2) Add this face as subentity for interior element to be communicated as ghost
     * 3) Add provisional ghost elements for communication
     *
     *
     * */
    void ghostComputeLocalNeighborGhostElements()
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Computing to-be ghost elements locally");

        typedef typename std::map<LocalIndexType, std::vector<InternalIndexType> > TmpMap;
        typedef typename std::map<LocalIndexType, std::vector<InternalIndexType> >::iterator TmpMapIterator;

        std::vector<std::map<LocalIndexType, LocalIndexType> > boundaryInternalLocalIndex2TmpIndex(size_);
        std::vector< std::vector<InternalIndexType> > boundaryInternalAssocSubentityIndex;


        // Determine boundary internal elements associated with each periodic boundary, if there are at all periodic boundaries
        if (gridstorage_.periodicCuboidDimensions_.size() > 0) {
        	assert(periodicConstr_ != nullptr);  // One must provide the periodic map if periodicity is used

        	//const GlobalIndexPeriodicScatterDataMap & periodicMap = periodicConstr_->map();

			for (auto const & faceIter: periodicConstr_->map())
			{
				// Get global and local face index
				GlobalIndexType thisFaceGlobalIndex = faceIter.first;
				auto thisFaceIndexIter = gridstorage_.entityIndexMap_[FACE_CODIM].find(thisFaceGlobalIndex);
				assert(thisFaceIndexIter != gridstorage_.entityIndexMap_[FACE_CODIM].end());
				LocalIndexType thisFaceLocalIndex = thisFaceIndexIter->second;

				// Get neighbor face index and rank
				int thisNeighborFaceRank = faceIter.second.ownerRank_;
				GlobalIndexType thisNeighborFaceGInd = faceIter.second.gind_;

				// Prepare the face for global communication of non-local periodic neighbors
				auto thisNeighborFaceIndexIter = gridstorage_.entityIndexMap_[FACE_CODIM].find(thisNeighborFaceGInd);
				if (thisNeighborFaceIndexIter == gridstorage_.entityIndexMap_[FACE_CODIM].end()) {
					LocalIndexType  thisGhostLocalIndex = gridstorage_.face_[thisFaceLocalIndex].element1Index;
					InternalIndexType thisFaceSubentityIndex = gridstorage_.face_[thisFaceLocalIndex].element1SubentityIndex;

		            assert(thisGhostLocalIndex >= 0);
		            assert(thisFaceSubentityIndex >= 0);

					std::vector<InternalIndexType> assocFaceSubentityIndex;

					auto tmpIter = boundaryInternalLocalIndex2TmpIndex[thisNeighborFaceRank].find(thisGhostLocalIndex);
					if (tmpIter != boundaryInternalLocalIndex2TmpIndex[thisNeighborFaceRank].end()) {
						boundaryInternalAssocSubentityIndex[tmpIter->second].push_back(thisFaceSubentityIndex);
					} else {
						nPeriodicGhostNeighbor_++;
						boundaryInternalAssocSubentityIndex.push_back(std::vector<InternalIndexType>(1, thisFaceSubentityIndex));
						boundaryInternalLocalIndex2TmpIndex[thisNeighborFaceRank][thisGhostLocalIndex] = boundaryInternalAssocSubentityIndex.size() - 1;
					}
				}
			}
        }


        // Determine boundary internal elements associated with each periodic boundary
        for (auto const & faceIter : gridstorage_.processBoundaryIndexMap_[FACE_CODIM])
        {
            LocalIndexType thisFaceLocalIndex = faceIter.first;
            LocalIndexType thisFaceLocalPBIndex = faceIter.second;
            int thisNeighborFaceRank = gridstorage_.PB2PBNeighborRank_[FACE_CODIM][thisFaceLocalPBIndex][0];
            GlobalIndexType thisFaceGlobalIndex = gridstorage_.face_[thisFaceLocalIndex].globalIndex;
            LocalIndexType  thisGhostLocalIndex = gridstorage_.face_[thisFaceLocalIndex].element1Index;
            InternalIndexType thisFaceSubentityIndex = gridstorage_.face_[thisFaceLocalIndex].element1SubentityIndex;

            assert(thisGhostLocalIndex >= 0);
            assert(thisFaceSubentityIndex >= 0);

            // Check that the face is marked as interior
            assert(gridstorage_.face_[thisFaceLocalIndex].ptype == Dune::PartitionType::BorderEntity);

            std::vector<InternalIndexType> assocFaceSubentityIndex;

            auto tmpIter = boundaryInternalLocalIndex2TmpIndex[thisNeighborFaceRank].find(thisGhostLocalIndex);
            if (tmpIter != boundaryInternalLocalIndex2TmpIndex[thisNeighborFaceRank].end()) {
            	boundaryInternalAssocSubentityIndex[tmpIter->second].push_back(thisFaceSubentityIndex);
            } else {
            	nProcessGhost_++;

            	boundaryInternalAssocSubentityIndex.push_back(std::vector<InternalIndexType>(1, thisFaceSubentityIndex));
            	boundaryInternalLocalIndex2TmpIndex[thisNeighborFaceRank][thisGhostLocalIndex] = boundaryInternalAssocSubentityIndex.size() - 1;
            }
        }

        std::cout << "On rank " << rank_
				<< ", periodic ghosts " << nPeriodicGhostNeighbor_
				<< ", process ghosts " << nProcessGhost_
				<< std::endl;


        // Associate the prospective Ghosts with respective neighbor processes
        for (int iProc = 0; iProc < size_; iProc++) {
        	for (auto && gelemIter : boundaryInternalLocalIndex2TmpIndex[iProc]) {
        		assert(gelemIter.second < boundaryInternalAssocSubentityIndex.size());

        		thisProcessNeighborGhostLocalIndex_[iProc].push_back(gelemIter.first);
        		thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_[iProc].push_back(boundaryInternalAssocSubentityIndex[gelemIter.second]);
        	}
        }
    }


    /** Communicate to each process a list of interpolation orders of the ghost elements it is going to receive
     *
     * Algorithm:
     *
     * 1) Communicate number of ghosts to send to each other process
     * 2) For each ghost element communicate its interpolation order and the number of associated faces
     *
     *
     * Optimization Proposal:
     * Same as for FaceGlobalIndex algorithm
     *
     * */
    void ghostDistributePreliminaries()
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicate ghost element interpolatory orders and neighbor face count");


    	// 2) For each ghost element communicate its interpolation order and the number of associated faces
    	// *****************************************************************
        std::vector<int> sendbuf, sendcounts;
        std::vector<int> recvbuf, recvcounts;


        for (int iProc = 0; iProc < size_; iProc++)
        {
        	for (int iGhost = 0; iGhost < thisProcessNeighborGhostLocalIndex_[iProc].size(); iGhost++)
        	{
            	LocalIndexType thisElemLocalIndex = thisProcessNeighborGhostLocalIndex_[iProc][iGhost];
            	int thisElemPBNeighbors = thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_[iProc][iGhost].size();
            	sendbuf.push_back(gridstorage_.element_[thisElemLocalIndex].interpOrder);
            	sendbuf.push_back(thisElemPBNeighbors);
        	}
            sendcounts.push_back(2 * thisProcessNeighborGhostLocalIndex_[iProc].size());
        }

        // Perform MPI_Alltoallv
        // [FIXME] Use neighbor communication, it is already implemented in allcomm_
        allcomm_.all2all(sendbuf, sendcounts, recvbuf, recvcounts);

    	// 3) Parse the received data.
    	// *****************************************************************
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++) {
        	assert(recvcounts[iProc] % 2 == 0);  // 2 objects sent for each face
            for (int iElem = 0; iElem < recvcounts[iProc] / 2; iElem++) {
            	assert(iData+1 < recvbuf.size());
            	neighborProcessGhostInterpOrder_[iProc].push_back(recvbuf[iData++]);
            	neighborProcessNAssociatedFace_[iProc].push_back(recvbuf[iData++]);
            }
        }
    }


    /** Communicate to all ghost element information except of the explicit vertex coordinates
     * MPI_alltoallv - Package
     *      - element globalIndex
     *      - elementPhysicalTag
     *      - all associated PBFace subentity indices
     *      - all interpVertex globalIndices
     *      - all edge globalIndices
     *      - all face globalIndices
     *
     * Optimization Proposal:
     * - It is possible not to send the associated face global index, and instead figure out the order of elements
     * according to the order of associated FaceKeys. The gain however is very small
     * - It is possible not to send all edges and vertices, but only those that are not present on the element,
     * based on subentity index of the face. For high orders communication gain is significant. However, complicated
     * pre and post processing required to figure out the exact missing entities and then reassemble all together.
     * Further complication because there may be more than 1 face associated to an element.
     * Further complication because for periodic faces need to send everything
     *
     * TODO: If planning to use with non-tetrahedral meshes, need to pass element type as well
     *
     * */
    void ghostDistributeGhostElements(
    		std::vector<int> & recvPackageGhostElementData,
			std::vector<int> & recvPackageGhostElementCounts)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicate ghost element data");

    	std::vector<int> sendPackageGhostElementData;
    	std::vector<int> sendcounts(size_, 0);

    	// [TODO] Swap explicit integers with corresponding geometryType functions
        int nSubentityEdge = 6;
        int nSubentityFace = 4;

        for (int iProc = 0; iProc < size_; iProc++)
        {
            // Assemble the array to send
            // *****************************************************************
            for (unsigned int iElem = 0; iElem < thisProcessNeighborGhostLocalIndex_[iProc].size(); iElem++)
            {
                LocalIndexType ghostElementLocalIndex = thisProcessNeighborGhostLocalIndex_[iProc][iElem];

                // Compute SendPackage size
                //*****************************************
                int nThisNeighborPBFace = thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_[iProc][iElem].size();
                int thisDofNum = gridstorage_.element_[ghostElementLocalIndex].vertexIndexSet.size();
                sendcounts[iProc] += 2 + nThisNeighborPBFace + nSubentityFace + nSubentityEdge + thisDofNum;

                // Package element data
                //*****************************************
                sendPackageGhostElementData.push_back(gridstorage_.element_[ghostElementLocalIndex].globalIndex);
                sendPackageGhostElementData.push_back(gridstorage_.element_[ghostElementLocalIndex].physicalTag);

                // Package subentity PBFace Internal Indices
                //*****************************************
                for (int iFace = 0; iFace < nThisNeighborPBFace; iFace++)
                {
                	sendPackageGhostElementData.push_back(thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_[iProc][iElem][iFace]);
                }

                // Package subentity Face Global Indices
                //*****************************************
                if (gridstorage_.elementSubentityCodim1_[ghostElementLocalIndex].size() != nSubentityFace)  {
                	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: elementSubentityCodim1_ - an element has unexpected number of subentity faces");
                	DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: elementSubentityCodim1_ - an element has unexpected number of subentity faces");
                }
                for (int iFace = 0; iFace < nSubentityFace; iFace++)
                {
                	LocalIndexType localFaceIndex = gridstorage_.elementSubentityCodim1_[ghostElementLocalIndex][iFace];
                	//LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " /////////////// element=" + std::to_string(ghostElementLocalIndex) + " sending face=" + std::to_string(localFaceIndex));
                	sendPackageGhostElementData.push_back(gridstorage_.face_[localFaceIndex].globalIndex);
                }

                // Package subentity Edge Global Indices
                //*****************************************
                if (gridstorage_.elementSubentityCodim2_[ghostElementLocalIndex].size() != nSubentityEdge)  {
                	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: elementSubentityCodim2_ - an element has unexpected number of subentity edges");
                	DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: elementSubentityCodim1_ - an element has unexpected number of subentity edges");
                }
                for (int iEdge = 0; iEdge < nSubentityEdge; iEdge++)
                {
                	LocalIndexType  localEdgeIndex = gridstorage_.elementSubentityCodim2_[ghostElementLocalIndex][iEdge];
                	GlobalIndexType globalEdgeIndex = gridstorage_.edge_[localEdgeIndex].globalIndex;
                	//LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " /////////////// element=" + std::to_string(ghostElementLocalIndex) + " sending edge localIndex=" + std::to_string(localEdgeIndex) + " globalIndex=" + std::to_string(globalEdgeIndex));

                	sendPackageGhostElementData.push_back(globalEdgeIndex);
                }

                // Package subentity vertex data
                //*****************************************
                for (int iDof = 0; iDof < thisDofNum; iDof++) {
                	LocalIndexType localVertexIndex = gridstorage_.element_[ghostElementLocalIndex].vertexIndexSet[iDof];
                	sendPackageGhostElementData.push_back(gridstorage_.point_[localVertexIndex].globalIndex);
                }
            }
        }

        //std::cout << "CommRequest: send="<< VectorHelper::vector2string(sendPackageGhostElementData) << " recvsize=" << totalRecvSize << std::endl;

        allcomm_.all2all(sendPackageGhostElementData, sendcounts, recvPackageGhostElementData, recvPackageGhostElementCounts);

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicated ghost elements sendcounts=(" + VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + VectorHelper::vector2string(recvPackageGhostElementCounts) + ")" );
    }


    /** Add received elements to the mesh. For each vertex global index, find if coordinate is already present on this process
     *  If not, mark this vertex as a missing vertex for further communication.
     *
     *  1) Add all data on this ghost element to ghost element array
     *  2) Map global ghost element index to local ghost element array index
     *  3) Add local ghost element index as a neighbor to corresponding process boundary face
     *  3.1) Check if the global index of the face is found on this process -> processBoundary, or periodicMap -> periodicBoundary
     *
     *
     * */
    void ghostInsertGhostElements (
            const std::vector<int> & packageGhostElementData,
			const std::vector<int> & packageGhostElementCounts,
            std::vector<std::vector<GlobalIndexType> > & missingVertices)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Inserting communicated ghost elements");

    	// [TODO] Replace dodgy periodicity check with a member function of gridbase, e.g. if (gridbase_.havePeriodic()) {...}
    	bool havePeriodic = gridstorage_.periodicCuboidDimensions_.size() > 0;

        // [TODO] Replace by ReferenceElement routines
        int nSubentityEdge = 6;
        int nSubentityFace = 4;

        Dune::GeometryType meshGeometryType=Dune::GeometryTypes::simplex(dimension);

        int countMarkedPeriodicFaces = 0;
        int countMarkedAllBoundaryFaces = 0;

        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
        	std::set<GlobalIndexType> missingVerticesFromThisProcess;

            for (unsigned int iGhost = 0; iGhost < neighborProcessGhostInterpOrder_[iProc].size(); iGhost++)
            {
            	// Unpack element data
            	// ***********************************************************************
                EntityStorage ghostElement;
                ghostElement.geometryType = meshGeometryType;
                ghostElement.globalIndex = packageGhostElementData[iData++];
                ghostElement.ptype       = Dune::PartitionType::GhostEntity;
                ghostElement.interpOrder = neighborProcessGhostInterpOrder_[iProc][iGhost];
                ghostElement.physicalTag = packageGhostElementData[iData++];

                std::vector<InternalIndexType> ghostCommunicatingFaceSubentityIndexSet;
                for (int iFace = 0; iFace < neighborProcessNAssociatedFace_[iProc][iGhost]; iFace++)
                {
                	int faceElementInternalIndex = packageGhostElementData[iData++];
                	assert((faceElementInternalIndex>= 0)&&(faceElementInternalIndex < nSubentityFace));
                	ghostCommunicatingFaceSubentityIndexSet.push_back(faceElementInternalIndex);
                }


                // Create the ghost element and insert it into global map
                // ***********************************************************************
                LocalIndexType ghostElementLocalIndex = gridstorage_.element_.size();
                gridstorage_.entityIndexMap_[ELEMENT_CODIM][ghostElement.globalIndex] = ghostElementLocalIndex;


                // Unpack face data. Create ghost faces, note them as the ghost element subentities
                // ***********************************************************************
                gridstorage_.elementSubentityCodim1_.push_back(std::vector<LocalIndexType>());
                for (int iFace = 0; iFace < nSubentityFace; iFace++)
                {
                	GlobalIndexType ghostFaceGlobalIndex = packageGhostElementData[iData++];
                	LocalIndexType thisFaceLocalIndex;
                	Global2LocalIterator ghostFaceIter = gridstorage_.entityIndexMap_[FACE_CODIM].find(ghostFaceGlobalIndex);

                	// If the face already exists, then it is one of the process boundaries, otherwise it is a new ghost face
                	if (ghostFaceIter == gridstorage_.entityIndexMap_[FACE_CODIM].end())
                	{
                		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Face globalIndex=" + std::to_string(ghostFaceGlobalIndex) + " is new");

                    	FaceStorage ghostFace;
                    	ghostFace.geometryType=Dune::GeometryTypes::triangle;
                    	ghostFace.globalIndex              = ghostFaceGlobalIndex;
                    	ghostFace.ptype                    = Dune::PartitionType::GhostEntity;
                    	ghostFace.boundaryType             = GridStorageType::FaceBoundaryType::None;

                    	ghostFace.element1Index            = ghostElementLocalIndex;
                    	ghostFace.element2Index            = -1;      // Neighbor unknown, except for shared face, which is assigned later
                    	ghostFace.element1SubentityIndex   = iFace;  // Faces should be communicated in the correct subentity order
                    	ghostFace.element2SubentityIndex   = -1;  // Neighbor unknown, except for shared face, which is assigned later
                    	ghostFace.physicalTag              = -1;      // By convention, out of all ghost entities, only element tags are communicated. Update when necessary

                    	thisFaceLocalIndex = gridstorage_.face_.size();
                    	gridstorage_.face_.push_back(ghostFace);
                    	gridstorage_.entityIndexMap_[FACE_CODIM][ghostFaceGlobalIndex] = thisFaceLocalIndex;
                	} else
                	{
                		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Face globalIndex=" + std::to_string(ghostFaceGlobalIndex) + " already exists");
                		thisFaceLocalIndex = ghostFaceIter->second;
                	}

                	// Note this face as the ghost element subentity
                	gridstorage_.elementSubentityCodim1_[ghostElementLocalIndex].push_back(thisFaceLocalIndex);
                }


                // Unpack edge data. Create ghost edges, note them as the ghost element subentities
                // ***********************************************************************
                gridstorage_.elementSubentityCodim2_.push_back(std::vector<LocalIndexType>());
                for (int iEdge = 0; iEdge < nSubentityEdge; iEdge++)
                {
                	GlobalIndexType thisEdgeGlobalIndex = packageGhostElementData[iData++];
                	LocalIndexType thisEdgeLocalIndex;
                	Global2LocalIterator thisEdgeIter = gridstorage_.entityIndexMap_[EDGE_CODIM].find(thisEdgeGlobalIndex);

                	// If the edge already exists, then it is one of the process boundaries, otherwise it is a new ghost edge
                	if ( thisEdgeIter == gridstorage_.entityIndexMap_[EDGE_CODIM].end() )
                	{
                		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Edge globalIndex=" + std::to_string(thisEdgeGlobalIndex) + " is new");

                    	EdgeStorage thisEdge;
                		thisEdge.globalIndex      = thisEdgeGlobalIndex;
                		thisEdge.ptype            = Dune::PartitionType::GhostEntity;
                		thisEdge.elementIndex     = ghostElementLocalIndex;
                		thisEdge.subentityIndex   = iEdge;  // Edges should be communicated in the correct subentity order

                    	thisEdgeLocalIndex = gridstorage_.edge_.size();
                    	gridstorage_.edge_.push_back(thisEdge);
                    	gridstorage_.entityIndexMap_[EDGE_CODIM][thisEdgeGlobalIndex] = thisEdgeLocalIndex;
                	} else
                	{
                		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Edge globalIndex=" + std::to_string(thisEdgeGlobalIndex) + " already exists");
                		thisEdgeLocalIndex= (*thisEdgeIter).second;
                	}

                	// Note this face as the ghost element subentity
                	gridstorage_.elementSubentityCodim2_[ghostElementLocalIndex].push_back(thisEdgeLocalIndex);
                }


                // Read DoF of this element
                // ***********************************************************************
                int thisElementDof = CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, ghostElement.interpOrder);
                for (int iDof = 0; iDof < thisElementDof; iDof++)
                {
                	GlobalIndexType thisVertexGlobalIndex = packageGhostElementData[iData++];
                    Global2LocalIterator vertexIter = gridstorage_.entityIndexMap_[VERTEX_CODIM].find(thisVertexGlobalIndex);

                    // If this vertex already exists, just reuse it. Otherwise, create new vertex, and later request its coordinate
                    if (vertexIter != gridstorage_.entityIndexMap_[VERTEX_CODIM].end()) {
                    	LocalIndexType thisVertexLocalIndex = (*vertexIter).second;
                        ghostElement.vertexIndexSet.push_back(thisVertexLocalIndex);

                        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Vertex already on this process GlobalIndex=" + std::to_string(thisVertexGlobalIndex));
                    }
                    else
                    {
                    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Ghost Vertex missing on this process GlobalIndex=" + std::to_string(thisVertexGlobalIndex));

                        // Create a new vertex with local index pointing to the end of current vertex array
                        LocalIndexType localVertexIndex = gridstorage_.point_.size();
                    	ghostElement.vertexIndexSet.push_back(localVertexIndex);

                        // Insert the fake vertex into the mesh
                        insertFakeVertex(thisVertexGlobalIndex, Dune::PartitionType::GhostEntity);

                        // Note that this vertex needs communicating
                        missingVerticesFromThisProcess.insert(thisVertexGlobalIndex);
                    }
                }

                // Add the element to the array
                gridstorage_.element_.push_back(ghostElement);

                // Loop over shared faces (process boundaries and periodic)
                // Associate all relevant boundary faces with this ghost elements
                // ***********************************************************************
                for (const auto & ghostFaceSubIndex : ghostCommunicatingFaceSubentityIndexSet)
                {
                	LocalIndexType    ghostFaceLocalIndex = gridstorage_.elementSubentityCodim1_[ghostElementLocalIndex][ghostFaceSubIndex];
                	GlobalIndexType   ghostFaceGlobalIndexInNeighbor = gridstorage_.face_[ghostFaceLocalIndex].globalIndex;

                	// For process boundaries, the global index of shared face is the same as seen from both neighbors
                	GlobalIndexType sharedFaceGlobalIndexThisProcess = ghostFaceGlobalIndexInNeighbor;

                	// Check if face is periodic, in this case the global index of the shared face on this process is different than that on the sender
                	bool isFacePeriodic = false;
                	if (havePeriodic) {
                    	auto facePeriodicIter = periodicConstr_->mapInverse().find(ghostFaceGlobalIndexInNeighbor);
                    	isFacePeriodic =  (facePeriodicIter != periodicConstr_->mapInverse().end());

                    	countMarkedAllBoundaryFaces++;

                    	if (isFacePeriodic) {
                    		countMarkedPeriodicFaces++;

                    		sharedFaceGlobalIndexThisProcess = facePeriodicIter->second.gind_;

                    		// For communicated periodic faces the opposite face should not be local to this process
                    		// Note that periodic faces local to this process are processed separately in the first routine of this class
                    		// NOTE: Nope, it already is local because we have just inserted it above when creating the ghost face
                    		//auto faceG2LiterNeighbor = gridstorage_.entityIndexMap_[FACE_CODIM].find(ghostFaceGlobalIndexInNeighbor);
                    		//assert(faceG2LiterNeighbor == gridstorage_.entityIndexMap_[FACE_CODIM].end());
                    	}
                	}

                	// Find local index of shared face as seen from inside
                	auto faceG2LiterThis = gridstorage_.entityIndexMap_[FACE_CODIM].find(sharedFaceGlobalIndexThisProcess);
                	assert(faceG2LiterThis != gridstorage_.entityIndexMap_[FACE_CODIM].end());  // The face local to this process should have a local index
                	GlobalIndexType localFaceLocalIndex = faceG2LiterThis->second;

                	if (isFacePeriodic) {
                    	GlobalIndexType localElementLocalIndex = gridstorage_.face_[localFaceLocalIndex].element1Index;
                    	InternalIndexType localFaceSubIndex = gridstorage_.face_[localFaceLocalIndex].element1SubentityIndex;

                    	// Associate the ghost element face with the local element (Needed for DataHandle)
                    	// IMPORTANT NOTE: Only a periodic ghost face can have a defined neighbor
                    	//   Regular ghosts do not have neighbors, and the periodic boundary is already shared by interior and ghost, so it contains all the information the ghost needs in this case
                    	gridstorage_.face_[ghostFaceLocalIndex].element2Index = localElementLocalIndex;
                    	gridstorage_.face_[ghostFaceLocalIndex].element2SubentityIndex = localFaceSubIndex;  // Note the subentity index in outside
                	} else {
                		// If this is a process boundary face, then the local index of the shared face should be the same as seen from interior and ghost elements
                    	if (ghostFaceLocalIndex != localFaceLocalIndex) {
    						LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Error: Received Ghost process boundary face not found among faces of this process");
    						DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: Received Ghost process boundary face not found among faces of this process");
                        }
                	}

                	// Associate the local face with the ghost element
                	gridstorage_.face_[localFaceLocalIndex].element2Index = ghostElementLocalIndex;
                	gridstorage_.face_[localFaceLocalIndex].element2SubentityIndex = ghostFaceSubIndex;  // Note the subentity index in outside
                }
            }

            // Rewrite missing vertices from a set to an array, such that the iteration order of elements is always the same
            typedef typename std::set<GlobalIndexType>::iterator  TmpSetIterator;
            for (TmpSetIterator tmpIter = missingVerticesFromThisProcess.begin(); tmpIter != missingVerticesFromThisProcess.end(); tmpIter++)
            {
            	missingVertices[iProc].push_back(*tmpIter);
            }
        }

        std::cout << "On rank " << rank_ << " periodic ghosts communicated " << countMarkedPeriodicFaces
        		<< " out of total " << countMarkedAllBoundaryFaces << std::endl;
    }


    /** Communicate to each process the number of missing vertices out of the ones it had provided with Ghost Elements
     *  Then communicate the globalIndices of all missing vertices
     *
     * */
    void ghostCommunicateMissingVertexGlobalIndices(
    		const std::vector<std::vector<GlobalIndexType> > & missingVertices,
            std::vector<int> & packageMissingVertexGlobalIndices,
            std::vector<int> & verticesRequestedByThis,
            std::vector<int> & verticesToSendByThis)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicating ghost element missing vertex indices");

        // 4.2) MPI_alltoallv - tell each process the list of global Indices of coordinates you want from it
        // Cleanup
        std::vector<int> sendcounts(size_, 0);
        std::vector<int> recvcounts;
        std::vector<int> missingVertexGlobalIndexRequested;

        for (int iProc = 0; iProc < size_; iProc++)
        {
            sendcounts[iProc] = missingVertices[iProc].size();

            for (int iVert = 0; iVert < sendcounts[iProc]; iVert++)  {
            	missingVertexGlobalIndexRequested.push_back(missingVertices[iProc][iVert]);
            }
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Missing vertex numbers= " + VectorHelper::vector2string(missingVertexGlobalIndexRequested) + " sendcounts=(" + VectorHelper::vector2string(sendcounts) + ") recvcounts=(" + VectorHelper::vector2string(recvcounts) + ")" );

        allcomm_.all2all(missingVertexGlobalIndexRequested, sendcounts, packageMissingVertexGlobalIndices, recvcounts);

        // We will require the information about requested and sent vertices when we communicate the coordinates
        verticesRequestedByThis.swap(sendcounts);
        verticesToSendByThis.swap(recvcounts);
    }


    // Distrubute vertex coordinates and add received coordinates to the mesh
    void ghostCommunicateMissingVertexCoordinates (
            const std::vector<std::vector<GlobalIndexType> > & missingVertices,
			const std::vector<int> & packageMissingVertexGlobalIndices,
			const std::vector<int> & verticesToSendByThis,
			const std::vector<int> & verticesToReceiveByThis  // [TODO] PERIODIC: Assert this vec is equivalent to the send vector in the following iteration
    )
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Communicating ghost element missing vertex coordinates");

        // 4.3) MPI_alltoallv - package list of globalId+coordinate for each process and send it
        std::vector<int> sendcounts(size_);
        std::vector<int> recvcounts(size_);
        std::vector<double> recvbuf, sendbuf;

        int iData = 0;

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: verticesToSend=(" + VectorHelper::vector2string(verticesToSendByThis) + ") vertices to receive=(" + VectorHelper::vector2string(verticesToReceiveByThis) + ")" + " missingGlobalIndicesPackage=(" + VectorHelper::vector2string(packageMissingVertexGlobalIndices) + ")" );


        for (int iProc = 0; iProc < size_; iProc++)
        {
            // Go through all vertices requested from this process. Package coordinates
            for (int j = 0; j < verticesToSendByThis[iProc]; j++)
            {
            	GlobalIndexType thisVertexGlobalIndex = packageMissingVertexGlobalIndices[iData++];
            	LocalIndexType thisVertexLocalIndex = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisVertexGlobalIndex];

                GlobalCoordinate p = gridstorage_.point_[thisVertexLocalIndex].coord;
                for (int iDim = 0; iDim < dimension; iDim++)  { sendbuf.push_back(p[iDim]); }

                //std::cout << "process_" << rank_ << " sending to process " << iProc << " a requested vertex " << thisVertexGlobalIndex << " with coord " << p << std::endl;
            }

            // We communicate (coord = 3 doubles) for each sent/received vertex
            // We now receive the amount we sent before, and send the amount we received before
            sendcounts[iProc] = 3 * verticesToSendByThis[iProc];
        }

        allcomm_.all2all(sendbuf, sendcounts, recvbuf, recvcounts);

        //MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        //MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE, reinterpret_cast<double*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_DOUBLE, comm );

        // Assign coordinates to all missing vertices
        iData = 0;

        for (int iProc = 0; iProc < size_; iProc++)
        {
        	for (unsigned int iVert = 0; iVert < missingVertices[iProc].size(); iVert++)
            {
                GlobalCoordinate thisCoord;
                thisCoord[0] = recvbuf[iData++];
                thisCoord[1] = recvbuf[iData++];
                thisCoord[2] = recvbuf[iData++];

                GlobalIndexType thisVertexGlobalIndex = missingVertices[iProc][iVert];
                LocalIndexType thisVertexLocalIndex = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisVertexGlobalIndex];
                gridstorage_.point_[thisVertexLocalIndex].coord = thisCoord;

                //std::cout << "process_" << rank_ << " receiving requested vertex " << thisVertexGlobalIndex << " with coord " << thisCoord << std::endl;
            }
        }
    }


    void insertFakeVertex(GlobalIndexType globalIndex, Dune::PartitionType ptype = Dune::PartitionType::InteriorEntity)
    {
        VertexStorage point;
        point.globalIndex = globalIndex;
        point.ptype = ptype;

        gridstorage_.entityIndexMap_[VERTEX_CODIM][globalIndex] = gridstorage_.point_.size();
        gridstorage_.point_.push_back(point);
    }


private: // Private members

    // For each other process mark set of local indices of elements of this process which will become ghost elements
    std::vector< std::vector<LocalIndexType> > thisProcessNeighborGhostLocalIndex_;
    // For each other process mark set of associated boundary faces (process and periodic) of elements of this process which will become ghost elements
    std::vector< std::vector<std::vector<LocalIndexType> > > thisProcessNeighborGhostBoundaryFaceSubentityIndexSet_;
    // For each other process stores the set of interpolation orders of Ghost Elements that process wishes to communicate to this process
    std::vector< std::vector<InterpolatoryOrderType> > neighborProcessGhostInterpOrder_;
    // For each other process stores the set of numbers of PBFaces associated with each ghost element it is planning to send to this process
    std::vector< std::vector<int> > neighborProcessNAssociatedFace_;

    // MPI Communication wrapper
    AllCommunication allcomm_;


    // Map received from periodic boundary constructor
    const PeriodicConstr * periodicConstr_;

    int nPeriodicGhostNeighbor_;  // Count the number of periodic ghosts that are not located on this process
    int nProcessGhost_;  // Count the number of periodic ghosts that are not located on this process

    // Curvilinear Grid Storage Class
    GridStorageType & gridstorage_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARGHOSTCONSTRUCTOR_HH
