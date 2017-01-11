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

#ifndef DUNE_CURVILINEARPOSTCONSTRUCTOR_HH
#define DUNE_CURVILINEARPOSTCONSTRUCTOR_HH

#include <limits>
#include <map>
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

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>

#include <dune/curvilineargrid/utility/allcommunication.hh>


namespace Dune {

namespace CurvGrid {


// Forward declaration
//template <class ct, int cdim, bool isCached>
//class CurvilinearGridStorage;

//template <class ct, int cdim, bool isCached>
//class CurvilinearGridBase;



template <class GridBase>
class CurvilinearPostConstructor {
public:

    /* public types */
	typedef typename GridBase::ctype   ctype;
	static const int dimension = GridBase::dimension;

	typedef GridBase                            GridBaseType;
	typedef typename GridBase::GridStorageType  GridStorageType;

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

    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;

    typedef typename GridStorageType::EntityNeighborRankVector  EntityNeighborRankVector;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    // Face Structural Type
    static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;
    static const unsigned int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;
    static const unsigned int INTERIOR_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::InteriorBoundary;
    static const unsigned int PERIODIC_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::PeriodicBoundary;

public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearPostConstructor(
    		GridStorageType & gridstorage,
    		GridBaseType & gridbase,
    		MPIHelper &mpihelper) :
        gridstorage_(gridstorage),
        gridbase_(gridbase),
        mpihelper_(mpihelper),
		allcomm_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "Initialized CurvilinearPostConstructor");
    }


    // Generate
    // Loop over all corners of the mesh (including ghost) and add them to the interator set
    // Do this by iterating over all corners of all elements
    void generateCornerIndex()
    {
    	// Loop over all elements
    	for (unsigned int iElem = 0; iElem < gridstorage_.element_.size(); iElem++)
    	{
    		LoggingMessage::writePatience("Generating unique corner indices...", iElem, gridstorage_.element_.size());

        	// Get LocalCornerIndices
    		std::vector<LocalIndexType> thisCornerIndex = gridbase_.entityCornerLocalIndex(ELEMENT_CODIM, iElem);

        	// Loop over LocalCornerIndices
    		for (unsigned int iCorner = 0; iCorner < thisCornerIndex.size(); iCorner++)
    		{
    			// If this corner has not been added yet, add it
    			LocalIndexType thisIndex = thisCornerIndex[iCorner];
    			if (gridstorage_.cornerIndexMap_.find(thisIndex) == gridstorage_.cornerIndexMap_.end())
    			{
    				int thisCornerUniqueIndex = gridstorage_.cornerIndexMap_.size();
    				gridstorage_.cornerIndexMap_[thisIndex] = thisCornerUniqueIndex;
    				gridstorage_.cornerIndexMapRev_[thisCornerUniqueIndex] = thisIndex;
    			}
    		}
    	}
    }




    // Create iterator lists for quick iteration over mesh entities of a given type
    // NOTE: Must only insert corners, not other interpolatory vertices
    void generateIteratorSets()
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: Started generating iterator lists");

    	// Initialize all index sets
    	// ***************************************************************
    	gridstorage_.faceDomainBoundaryIndexSet_ = LocalIndexSet();
    	for (unsigned int iCodim = 0; iCodim <= dimension; iCodim++)
    	{
            gridstorage_.entityAllIndexSet_[iCodim] = LocalIndexSet();
            gridstorage_.entityInternalIndexSet_[iCodim] = LocalIndexSet();
            gridstorage_.entityProcessBoundaryIndexSet_[iCodim] = LocalIndexSet();
            gridstorage_.entityGhostIndexSet_[iCodim] = LocalIndexSet();

        	gridstorage_.entityDuneInteriorIndexSet_[iCodim] = LocalIndexSet();
        	gridstorage_.entityDuneInteriorBorderIndexSet_[iCodim] = LocalIndexSet();

    	}

    	// Generate iterator lists
    	// ***************************************************************
    	int cornerCount = 0;
        for (const auto & cornerIndexPair : gridstorage_.cornerIndexMap_) {
        	LoggingMessage::writePatience("Generating corner iterator list...", cornerCount++, gridstorage_.cornerIndexMap_.size());
        	fillPartitionIterator(VERTEX_CODIM, cornerIndexPair.first, gridstorage_.point_[cornerIndexPair.first].ptype);
        }

        for (LocalIndexType iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)         {
        	LoggingMessage::writePatience("Generating edge iterator list...", iEdge, gridstorage_.edge_.size());
        	fillPartitionIterator(EDGE_CODIM,    iEdge, gridstorage_.edge_[iEdge].ptype);
        }

        for (LocalIndexType iFace = 0; iFace < gridstorage_.face_.size(); iFace++)         {
        	LoggingMessage::writePatience("Generating face iterator list...", iFace, gridstorage_.face_.size());
        	fillPartitionIterator(FACE_CODIM,    iFace, gridstorage_.face_[iFace].ptype, gridstorage_.face_[iFace].boundaryType );
        }

        for (LocalIndexType iElem = 0; iElem < gridstorage_.element_.size(); iElem++)      {
        	LoggingMessage::writePatience("Generating element iterator list...", iElem, gridstorage_.element_.size());
        	fillPartitionIterator(ELEMENT_CODIM, iElem, gridstorage_.element_[iElem].ptype);
        }

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: Finished generating iterator lists");

    }


    /** \brief To minimise communication time, it is advantageous for ghost entities and
     *  internal entities next to boundary to have their own LocalIndex
     *
     *  Algorithm:
     *  1) Iterate over all PB faces
     *  1.1) Iterate over all subentities of internal neighbour (including itself),
     *       add them to map unless they already exist
     *  1.2) For each subentity, note neighbour rank of this PB Face
     *  1.3) Repeat the same for Ghost Neighbors
     *
     *  1) Find all Internal and Ghost entities associated with this PB Face
     *  2) Add entities to the map unless they are already added.
     *  3) Mark PB face neighbor rank on all subset entities
     *
     *  */
    void generateCommunicationMaps()
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: Started generating communication maps");

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// Resize all neighbor rank arrays, except the boundary internals, since we do not know their lengths
		// Boundary Internals will be resized as its number is calculated
		for (int iCodim = 0; iCodim <= dimension; iCodim++)
		{
			int nPB = (iCodim == 0) ? 0 : gridbase_.nEntity(iCodim, Dune::PartitionType::BorderEntity);
			int nG  =                     gridbase_.nEntity(iCodim, Dune::PartitionType::GhostEntity);

			// PB2PB already exists, so need not resize it
			gridstorage_.PB2GNeighborRank_[iCodim].resize(nPB);
			gridstorage_.G2BIPBNeighborRank_[iCodim].resize(nG);
			gridstorage_.G2GNeighborRank_[iCodim].resize(nG);
		}


		// Local part of the communication map construction
		//   - Construct communication maps I -> G and G->I
		//   - Contribute locally known part to PB->G
		computeIG(Dune::PartitionType::BorderEntity, NO_BOUNDARY_TYPE);
		computeIG(Dune::PartitionType::InteriorEntity, PERIODIC_BOUNDARY_TYPE);

		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: Finished generating communication maps");


	    /** \brief After this function, all entities that can be communicated over should be
	     * associated an array of ranks of all other processes over which this entity is shared
	     *
	     *
	     *  After the generateCommunicationMaps(), the neighbour ranks that could be computed locally
	     *  have already been added to the corresponding maps. That is, the protocols
	     *    * PB->PB
	     *    * I -> G
	     *    * G -> I
	     *  now have sufficient information, and the protocol
	     *    * PB->G
	     *  has candidates on all processes, but is incomplete
	     *
	     *  So this function must enable protocols
	     *    * PB->G, G->PB, G->G
	     *
	     *
	     *  Algorithm:
	     *  1) Iterate over all PB entities
	     *  1.1) divide provisional PB->G set by PB->PB set to see which provisional PB->G are new
	     *  1.2) Mark number of real PB-G candidates for each process
	     *  2) For all PB entities having non-zero PB->G, communicate G to all neighbouring PB
	     *
	     *  3) For all PB append received G by using union on them - This completes PB->G (hopefully)
	     *  4) For all PB entities having non-zero PB->G, communicate self to all G of (PB->G)
	     *  5) For all G append received PB by using union on them - This completes G->PB (hopefully)
	     *
	     *  6) For all PB entities having non-zero PB->G, communicate to all G all remaining G
	     *  6.1) Optimization - do this only if you are lowest rank among all PB-neighbors
	     *  6.2) Further optimization - do this only if you are modulus-rank among all PB-neighbors
	     *  7) For all G append received G by using union on them - This completes G->G (hopefully)
	     *
	     *
	     *  */
    	for (int iCodim = 0; iCodim <= dimension; iCodim++)
    	{
    		LoggingMessage::writePatience("Communicating neighbour ranks for different communication protocols...", iCodim, dimension + 1);

    		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: Started communicating entity ranks for codim=" + std::to_string(iCodim));


    		//1) Compute true PB-G candidates, communicate their number to PB neighbor entities
    		//2) Communicate PB-G candidates and fill them on the receiving end
    		// ************************************************************************************
    	    communicatePBG(iCodim, comm);


    	    //4) For all PB entities having non-zero PB->G, communicate self to all G of (PB->G)
    	    //5) For all G append received PB by using union on them - This completes G->PB (hopefully)
    		// ************************************************************************************
    	    communicateGPB(iCodim, comm);

    	    //6) For all PB entities having non-zero PB->G, communicate to all G all remaining G
    	    //7) For all G append received G by using union on them - This completes G->G (hopefully)
    	    // ************************************************************************************
    	    communicateGG(iCodim, comm);

    	    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: Finished communicating entity ranks for codim=" + std::to_string(iCodim));
    	}
    }


protected:





    // ************************************************************************************
	// Auxiliary methods for Iterator list generation
	// ************************************************************************************


    // Fill all LocalIndexSets necessary to iterate over the grid
    // Note: Must not contain all interpolatory vertices, only corners
    void fillPartitionIterator (int codim, LocalIndexType localIndex, PartitionType pitype, StructuralType bordertype = NO_BOUNDARY_TYPE)
    {
    	// Check if the entity is of a valid type
    	gridbase_.assertValidCodimStructuralType(codim, pitype);

    	// All-set includes entities of all structural types
        gridstorage_.entityAllIndexSet_[codim].insert(localIndex);

    	switch (pitype)
    	{
    	case Dune::PartitionType::InteriorEntity  :
    	{
        	if (bordertype == DOMAIN_BOUNDARY_TYPE) {
        		assert(codim == FACE_CODIM);  // Only faces can be boundarySegments
        		gridstorage_.faceDomainBoundaryIndexSet_.insert(localIndex);
        	} else if (bordertype == PERIODIC_BOUNDARY_TYPE) {
        		assert(codim == FACE_CODIM);  // Only faces can be boundarySegments
        		gridstorage_.facePeriodicBoundaryIndexSet_.insert(localIndex);
        		gridstorage_.faceDomainBoundaryIndexSet_.insert(localIndex);  // Periodic boundaries are also domain boundaries
        	} else if (bordertype == INTERIOR_BOUNDARY_TYPE) {
        		assert(codim == FACE_CODIM);  // Only faces can be boundarySegments
        		gridstorage_.faceInteriorBoundaryIndexSet_.insert(localIndex);
        	} else {
        		// Consider a face internal only if it is not a boundary segment
        		// Consider all other codim interior entities internal
        		gridstorage_.entityInternalIndexSet_[codim].insert(localIndex);
        	}

    		gridstorage_.entityDuneInteriorIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorBorderIndexSet_[codim].insert(localIndex);
    	} break;
    	case Dune::PartitionType::BorderEntity    :
    	{
    		gridstorage_.entityProcessBoundaryIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorBorderIndexSet_[codim].insert(localIndex);
    	} break;
    	case Dune::PartitionType::GhostEntity     :
    	{
    		gridstorage_.entityGhostIndexSet_[codim].insert(localIndex);
    	} break;

    	}
    }


    // Extracts the locally available data to construct the protocols I -> G, G ->I and, partially, PB->G
    // This is possible, because the PB<->PB is generated during global index construction, and
    //    I -> G and G->I are effectively equivalent to PB<->PB, given that ghost elements exist
    void computeIG(Dune::PartitionType ptype, StructuralType boundarytype) {

    	// [FIXME] GET DATA FROM REF SUBENTITY SIZE
    	// [FIXME] GET DATA INDIVIDUALLY FOR EACH ELEMENT
		int nEdgeTriangle = 3;
		int nVertexTriangle = 3;
		int nTriangleTet = 4;
		int nEdgeTet = 6;
		int nVertexTet = 4;

		bool isPeriodic = boundarytype == PERIODIC_BOUNDARY_TYPE;

		// [FIXME] Iterate also over periodic
		int faceCount = 0;
		const LocalIndexSet & faceIndexSet = gridbase_.entityIndexSetSelect(FACE_CODIM , ptype, boundarytype);
		for (const auto & insideFaceLocalIndex : faceIndexSet)
		{
			LoggingMessage::writePatience("Generating communication maps...", faceCount++, gridbase_.nEntity(FACE_CODIM, Dune::PartitionType::BorderEntity));

			// Verify that the outer neighbor of the process boundary exists
			if (!gridbase_.checkFaceOuterNeighbor(insideFaceLocalIndex)) {
				std::stringstream logstr;
				logstr << "Supposed process boundary " << insideFaceLocalIndex << " does not have initialised outer neighbor" << std::endl;
				DUNE_THROW(Dune::IOError, logstr.str());
			}

			LocalIndexType insideElementLocalIndex = gridbase_.faceNeighbor(insideFaceLocalIndex, 0);
			LocalIndexType outsideElementLocalIndex = gridbase_.faceNeighbor(insideFaceLocalIndex, 1);

			Dune::PartitionType insideElementPartitionType = gridbase_.entityPartitionType(ELEMENT_CODIM, insideElementLocalIndex);
			Dune::PartitionType outsideElementPartitionType = gridbase_.entityPartitionType(ELEMENT_CODIM, outsideElementLocalIndex);

			// Avoid communication between local ghost pairs. It is completely unnecessary - any would be communicated from entity to itself on the same process
			bool isInteriorGhost = isPeriodic && (outsideElementPartitionType == Dune::PartitionType::InteriorEntity);
			if (!isInteriorGhost) {

				// Pair local face and element indices, as seen from both neighbors of the face
				// For PB it is the same face,
				// In periodic case these are different faces
				typedef std::pair<LocalIndexType, LocalIndexType> FaceElemLocalIndexPair;
				std::vector<FaceElemLocalIndexPair> neighborLocalIndex {
					FaceElemLocalIndexPair(insideElementLocalIndex, insideFaceLocalIndex),
					FaceElemLocalIndexPair(outsideElementLocalIndex, insideFaceLocalIndex)
				};


				if ((insideElementPartitionType != Dune::PartitionType::InteriorEntity) ||
						(outsideElementPartitionType != Dune::PartitionType::GhostEntity))
				{
					std::cerr << rank_ << " I->G communication constructor for a communication set with ptype="
							<< Dune::PartitionName(ptype) << ", boundarytype="
							<< boundarytype
							<< " found unexpected ptype pair=("
							<< Dune::PartitionName(insideElementPartitionType) << ", "
							<< Dune::PartitionName(outsideElementPartitionType) << ") " << std::endl;
					DUNE_THROW(Dune::IOError, "Unexpected partition type pair");
				}


				if (isPeriodic) {
					// NOTE: The inside face stores subindex info. Outside (ghost) face does not store anything
					InternalIndexType faceSubIndexInOuter = gridbase_.faceSubIndexInNeighbor(insideFaceLocalIndex, 1);
					neighborLocalIndex[1].second = gridbase_.subentityLocalIndex(outsideElementLocalIndex, ELEMENT_CODIM, FACE_CODIM, faceSubIndexInOuter);
				}

				// Find neighbor process rank
				int thisFaceNeighborRank = gridbase_.commEntityNeighborRankSet(FACE_CODIM, insideFaceLocalIndex, ptype, ptype, boundarytype)[0];


				// 1) Find all Internal and Ghost entities associated with this boundary face
				// **********************************************************************

				for (const auto & neighborLocalIndexPair : neighborLocalIndex) {

					LocalIndexType thisElementLocalIndex = neighborLocalIndexPair.first;
					LocalIndexType thisFaceLocalIndex = neighborLocalIndexPair.second;

					Dune::PartitionType thisElementPType = gridbase_.entityPartitionType(ELEMENT_CODIM, thisElementLocalIndex);

					// Find local indices of all subentities of the PB face, as seen from this neighbor element
					std::vector<LocalIndexType> thisFaceSubentityLocalIndex[4];
					thisFaceSubentityLocalIndex[FACE_CODIM].push_back(thisFaceLocalIndex);
					for (int i = 0; i < nEdgeTriangle; i++)    { thisFaceSubentityLocalIndex[EDGE_CODIM].push_back  (gridbase_.subentityLocalIndex(thisFaceLocalIndex, FACE_CODIM, EDGE_CODIM, i)); }
					for (int i = 0; i < nVertexTriangle; i++)  { thisFaceSubentityLocalIndex[VERTEX_CODIM].push_back(gridbase_.subentityLocalIndex(thisFaceLocalIndex, FACE_CODIM, VERTEX_CODIM, i)); }

					// Find local indices of all subentities of this element
					std::vector<LocalIndexType> thisElementSubentityLocalIndex[4];
					thisElementSubentityLocalIndex[ELEMENT_CODIM].push_back(thisElementLocalIndex);
					for (int i = 0; i < nTriangleTet; i++)  { thisElementSubentityLocalIndex[FACE_CODIM].push_back  (gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, FACE_CODIM, i)); }
					for (int i = 0; i < nEdgeTet; i++)      { thisElementSubentityLocalIndex[EDGE_CODIM].push_back  (gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, EDGE_CODIM, i)); }
					for (int i = 0; i < nVertexTet; i++)    { thisElementSubentityLocalIndex[VERTEX_CODIM].push_back(gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, VERTEX_CODIM, i)); }

					for (int iCodim = 0; iCodim <= dimension; iCodim++)
					{
						// Sort subentity index vectors
						std::sort(thisFaceSubentityLocalIndex[iCodim].begin(), thisFaceSubentityLocalIndex[iCodim].end());
						std::sort(thisElementSubentityLocalIndex[iCodim].begin(), thisElementSubentityLocalIndex[iCodim].end());

						// Subtract face subentity set from element subentity set, such that only internal and ghost subentities are left
						thisElementSubentityLocalIndex[iCodim] = VectorHelper::sortedSetComplement(thisElementSubentityLocalIndex[iCodim], thisFaceSubentityLocalIndex[iCodim]);

						// 2) Add entities to the map unless they are already added.
						// **********************************************************************

						Local2LocalMap & thisLocalCommMap = gridbase_.selectCommMap(iCodim, thisElementPType, NO_BOUNDARY_TYPE);


						for (unsigned int iEntity = 0; iEntity < thisElementSubentityLocalIndex[iCodim].size(); iEntity++ )
						{
							LocalIndexType thisEntityLocalIndex = thisElementSubentityLocalIndex[iCodim][iEntity];
							Dune::PartitionType thisEntityType = gridbase_.entityPartitionType(iCodim, thisEntityLocalIndex);
							StructuralType thisBoundaryType = iCodim == FACE_CODIM ?
									gridbase_.faceBoundaryType(thisEntityLocalIndex) :
									NO_BOUNDARY_TYPE;

							std::stringstream log_str;
							log_str << "CurvilinearPostConstructor: ---- Iterating over subentity ";
							log_str << " codim=" << iCodim;
							log_str << " subentityNo=" << iEntity;
							log_str << " localindex =" << thisEntityLocalIndex;
							log_str << " ptype =" << Dune::PartitionName(thisEntityType);
							log_str << " boundaryType =" << thisBoundaryType;
							LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_str.str());


							// If this is a PB Entity, that happens to not be on the face, it is possible
							// that it is on a face of a different process, then it is a PB->G link
							// Otherwise, this entity is either internal or ghost, and it contributes to one
							// of the new maps
							// [TODO] NOTE: Periodic <-> Ghost protocol not implemented at the moment
							if (thisEntityType == Dune::PartitionType::BorderEntity)
							{
								std::vector<int> & thisEntityPBNeighbors = gridbase_.commEntityNeighborRankSet(iCodim, thisEntityLocalIndex, Dune::PartitionType::BorderEntity, Dune::PartitionType::BorderEntity, NO_BOUNDARY_TYPE);
								bool isNewRank = !VectorHelper::isInside(thisEntityPBNeighbors, thisFaceNeighborRank);
								//std::cout << "testing isInside vector=(" << VectorHelper::vector2string(thisEntityPBNeighbors) << ") elem=" << thisFaceNeighborRank << " isNew=" << isNewRank << std::endl;

								// If this rank is not already in PB-PB of this entity, then it must be in PB-G
								if (isNewRank)
								{
									LocalIndexType thisEntityPBIndex = gridstorage_.processBoundaryIndexMap_[iCodim][thisEntityLocalIndex];
									gridstorage_.PB2GNeighborRank_[iCodim][thisEntityPBIndex].push_back(thisFaceNeighborRank);
								}
							} else {
								Local2LocalIterator thisIter = thisLocalCommMap.find(thisEntityLocalIndex);
								LocalIndexType thisEntitySubsetIndex;

								if (thisIter == thisLocalCommMap.end())
								{
									thisEntitySubsetIndex = thisLocalCommMap.size();
									thisLocalCommMap[thisEntityLocalIndex] = thisEntitySubsetIndex;

									// Enlarge the neighbour rank vector such that it has enough entries for all entities of subset
									if (thisElementPType == Dune::PartitionType::InteriorEntity)  {
										gridstorage_.BI2GNeighborRank_[iCodim].resize(thisEntitySubsetIndex + 1);
									}
								} else {
									thisEntitySubsetIndex = (*thisIter).second;
								}

								// 3) Mark PB face neighbor rank on all subset entities
								// Note: No communication of interior to interior
								// **********************************************************************
								if (thisElementPType == Dune::PartitionType::InteriorEntity)  {
									gridstorage_.BI2GNeighborRank_[iCodim][thisEntitySubsetIndex].push_back(thisFaceNeighborRank);
								} else {
									gridstorage_.G2BIPBNeighborRank_[iCodim][thisEntitySubsetIndex].push_back(thisFaceNeighborRank);
								}
							}
						}
					}
				}
			}
		}


		// Compactify I -> G, as this array will not be accessed at a later stage
		// Compactification means deleting repeating elements, because each neighbor rank set should only contain each of its neighbour ranks once
		// The reason for them being added more than one time, is because a boundary entity can be assigned a neighbour rank by more than one PB face
		// [TODO] Not sure if this is actually necessary at all. Are there repeating entities possible?
		for (int iCodim = 0; iCodim <= dimension; iCodim++) {
			for (auto && entityNbRanks : gridstorage_.BI2GNeighborRank_[iCodim])  { VectorHelper::compactify(entityNbRanks); }
			for (auto && entityNbRanks : gridstorage_.G2BIPBNeighborRank_[iCodim])  { VectorHelper::compactify(entityNbRanks); }
		}


    }



    // ************************************************************************************
	// Auxiliary methods for communicating communication entity neighbor ranks
	// ************************************************************************************

    /* \brief Construction of Process Boundary -> Ghost communication protocol
     *
     * [TODO]  Implement for periodic, whenever becomes necessary
     * [TODO] Possibly set complement unnecessary, since already only adding neighbor if it is not a neighbor
     *
     */
    void communicatePBG(int codim, MPI_Comm comm)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: -- Started ProcessBoundary-Ghost communication construction");

		//1) Iterate over all PB entities, find entities with non-zero PB-G candidates
		//   Mark these entities to send to all their PB neighbors
		// ************************************************************************************

    	// Entity number and candidate number per process
    	std::vector<int> nPBGEntitySend(size_);
    	std::vector<int> nPBGRankPerProcessSend(size_);

		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			// PB2G may contain repeating entities, and is not sorted, need to compactify
			VectorHelper::compactify(gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex]);

			//1.1) divide provisional PB->G set by PB->PB set to see which provisional PB->G are new
			gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex] = VectorHelper::sortedSetComplement(
				gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex],
				gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex]
			);

			// 1.2) Mark number of real PB-G candidates for each process
			int PBPBSize = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int realPBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();

			if (realPBGSize > 0)
			{
    			for (int i = 0; i < PBPBSize; i++)
    			{
    				int thisPBNeighborRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][i];
    				nPBGEntitySend[thisPBNeighborRank]++;
    				nPBGRankPerProcessSend[thisPBNeighborRank] += realPBGSize;
    			}
			}
		}


		//2) For each PB-G entity, communicate its global index
		// as well as communicate how many candidates it is going to send
		// ************************************************************************************

		// Construct displacements
		std::vector<int> displSendTmp2(size_);
		int sendSize2 = 0;

		for (int i = 0; i < size_; i++)
		{
			displSendTmp2[i] = (i == 0) ? 0 : displSendTmp2[i-1] + nPBGEntitySend[i-1];
			sendSize2 += nPBGEntitySend[i];
		}


		// Fill candidate number send arrays
    	std::vector<int> nPBGRankPerEntitySend(sendSize2);
    	std::vector<int> PBGEntityGlobalIndexSend(sendSize2);

		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalIndex = pbIndexPair.first;
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			// 1.2) Mark number of real PB-G candidates for each process
			int PBPBSize = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int realPBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();

			if (realPBGSize > 0)
			{
    			for (int i = 0; i < PBPBSize; i++)
    			{
    				int thisPBNeighborRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][i];
    				int thisTmpIndex = displSendTmp2[thisPBNeighborRank]++;

    				if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, PBGEntityGlobalIndexSend[thisTmpIndex]))
    				{
    					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found by its local index");
    				}
    				nPBGRankPerEntitySend[thisTmpIndex] = realPBGSize;
    			}
			}
		}


		// Communicate PBG entity global indices
//		MPI_Alltoallv (
//			                       PBGEntityGlobalIndexSend.data(),  nPBGEntitySend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(PBGEntityGlobalIndexRecv.data()), nPBGEntityRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);

    	std::vector<int> nPBGEntityRecv;
    	std::vector<int> PBGEntityGlobalIndexRecv;

		allcomm_.all2all(PBGEntityGlobalIndexSend, nPBGEntitySend, PBGEntityGlobalIndexRecv, nPBGEntityRecv);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated global indices");


		// Communicate candidate rank numbers
//		MPI_Alltoallv (
//				                   nPBGRankPerEntitySend.data(),  nPBGEntitySend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(nPBGRankPerEntityRecv.data()), nPBGEntityRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);


		std::vector<int> nPBGEntityRecv2;
		std::vector<int> nPBGRankPerEntityRecv;

		allcomm_.all2all(nPBGRankPerEntitySend, nPBGEntitySend, nPBGRankPerEntityRecv, nPBGEntityRecv2);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of ranks per entity");



		//3) Communicate PB-G candidate ranks
		// ************************************************************************************

		std::vector<int> displSendTmp3(size_);
		int sendSize3 = 0;

		for (int i = 0; i < size_; i++)
		{
			displSendTmp3[i] = (i == 0) ? 0 : displSendTmp3[i-1] + nPBGRankPerProcessSend[i-1];
			sendSize3 += nPBGRankPerProcessSend[i];
		}

		// Fill in communication arrays
		std::vector<int> neighborPBGRankSetSend(sendSize3);
		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			// Send all provisional PB-G ranks to all PB neighbors of this entity
			int PBPBSize = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int realPBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();


			for (int i = 0; i < realPBGSize; i++)
			{
				int thisPBGRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];

				for (int j = 0; j < PBPBSize; j++)
				{
    				int thisPBPBRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][j];
    				int tmpIndex = displSendTmp3[thisPBPBRank]++;

    				neighborPBGRankSetSend[tmpIndex] = thisPBGRank;
    			}
			}
		}

//		MPI_Alltoallv (
//			                       neighborPBGRankSetSend.data(),  nPBGRankPerProcessSend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(neighborPBGRankSetRecv.data()), nPBGRankPerProcessRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);

		std::vector<int> nPBGRankPerProcessRecv;
		std::vector<int> neighborPBGRankSetRecv;
		allcomm_.all2all(neighborPBGRankSetSend, nPBGRankPerProcessSend, neighborPBGRankSetRecv, nPBGRankPerProcessRecv);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated ranks");


		//4) Fill in
		// ************************************************************************************
		int iRankData = 0;
		int iEntityData = 0;

		for (int i = 0; i < size_; i++)
		{
			for (int j = 0; j < nPBGEntityRecv[i]; j++)
			{
				int nRankPerEntity = nPBGRankPerEntityRecv[iEntityData];
				GlobalIndexType thisEntityGlobalIndex = PBGEntityGlobalIndexRecv[iEntityData];
				iEntityData++;

				// Get local index corresponding to the communicated global index, check that it exists
				LocalIndexType thisEntityLocalIndex;
				if (!gridbase_.findEntityLocalIndex(codim, thisEntityGlobalIndex, thisEntityLocalIndex)) {
					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found corresponding to communicated global index");
				}

				// Check if the structural type of the received entity is PB
				StructuralType thisEntityType = gridbase_.entityPartitionType(codim, thisEntityLocalIndex);
				assert(thisEntityType == Dune::PartitionType::BorderEntity);

				LocalIndexType thisEntityLocalPBIndex = gridstorage_.processBoundaryIndexMap_[codim][thisEntityLocalIndex];
				for (int k = 0; k < nRankPerEntity; k++)
				{
					int thisNeighborRank = neighborPBGRankSetRecv[iRankData++];

					if (abs(thisNeighborRank) >= size_)
					{
						LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Error: Unexpected received rank=" + std::to_string(thisNeighborRank));
						assert(abs(thisNeighborRank) < size_);
					}

					gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].push_back(thisNeighborRank);

				}

			}
		}
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Filled in received data");


		//5) Compactify PB-G arrays (sort and eliminate repeating)
		// ************************************************************************************

		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			// Compactify the neighbor ranks
			VectorHelper::compactify(gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex]);

			// divide PB->G set by PB->PB set
			// This is just to make sure that no PB-PB links were picked up in the process
			gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex] = VectorHelper::sortedSetComplement(
				gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex],
				gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex]
			);
		}

		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: -- Finished ProcessBoundary-Ghost communication construction");
    }


    /** \brief For each PB entity that has non-zero PB-G, communicate its globalIndex and rank of self
     *
     * Algorithm:
     * 1) Communicate number of G-PB entities to send to each process
     * 2) Communicate global indices for each G-PB entity
     * 3) On receiving end, mark sender's rank on all received G-PB
     *
     * [TODO]  Implement for periodic, whenever becomes necessary
     *
     * */
    void communicateGPB(int codim, MPI_Comm comm)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: -- Started Ghost-ProcessBoundary communication construction");

		// 1) Communicate number of G-PB entities to send to each process
		// ********************************************************************
    	std::vector<int> nGPBEntitySend(size_);


		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			for (int i = 0; i < PBGSize; i++)
			{
				int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
				nGPBEntitySend[thisGNeighborRank]++;
			}
		}
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of Ghost candidates per process");


		// 2) Communicate global indices for each G-PB entity
		// ********************************************************************

		// Create displacement arrays
		std::vector<int> displSendTmp2(size_);
		int sendSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSendTmp2[i] = (i == 0) ? 0 : displSendTmp2[i-1] + nGPBEntitySend[i-1];
			sendSize += nGPBEntitySend[i];
		}

		// fill global index array
    	std::vector<int> GPBEntityGlobalIndexSend(sendSize);

		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalIndex = pbIndexPair.first;
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			for (int i = 0; i < PBGSize; i++)
			{
				int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
				int thisTmpIndex = displSendTmp2[thisGNeighborRank]++;

				// Check if element with this local index exists at all, otherwise bug in the map
				if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, GPBEntityGlobalIndexSend[thisTmpIndex])) {
					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found by its local index");
				}
			}
		}

		// Communicate
//		MPI_Alltoallv (
//				                   GPBEntityGlobalIndexSend.data(),  nGPBEntitySend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(GPBEntityGlobalIndexRecv.data()), nGPBEntityRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);

		std::vector<int> nGPBEntityRecv;
		std::vector<int> GPBEntityGlobalIndexRecv;

		allcomm_.all2all(GPBEntityGlobalIndexSend, nGPBEntitySend, GPBEntityGlobalIndexRecv, nGPBEntityRecv);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Sent own rank to all ghosts");


		// 3) On receiving end, mark sender's rank on all received G-PB
		// ********************************************************************
		int iData = 0;

		for (int i = 0; i < size_; i++)
		{
			for (int j = 0; j < nGPBEntityRecv[i]; j++)
			{
				GlobalIndexType thisEntityGlobalIndex = GPBEntityGlobalIndexRecv[iData++];

				// Get local index corresponding to the communicated global index, check that it exists
				LocalIndexType thisEntityLocalIndex;
				if (!gridbase_.findEntityLocalIndex(codim, thisEntityGlobalIndex, thisEntityLocalIndex)) {
					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found corresponding to communicated global index");
				}

				// Check if the structural type of the received entity is Ghost
				Dune::PartitionType thisEntityType = gridbase_.entityPartitionType(codim, thisEntityLocalIndex);
				if (thisEntityType != Dune::PartitionType::GhostEntity)
				{
					std::stringstream logstr;
					logstr << "CurvilinearPostConstructor: --   codim=" + codim;
					logstr << " entity globalIndex=" + thisEntityGlobalIndex;
					logstr << " expected type=" + Dune::PartitionName(Dune::PartitionType::GhostEntity);
					logstr << " received=" + Dune::PartitionName(thisEntityType);

					LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, logstr.str());
					assert(thisEntityType == Dune::PartitionType::GhostEntity);
				}

				LocalIndexType thisEntityLocalGhostIndex = gridstorage_.ghostIndexMap_[codim][thisEntityLocalIndex];
				gridstorage_.G2BIPBNeighborRank_[codim][thisEntityLocalGhostIndex].push_back(i);
			}
		}
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   filled received data");


		// 4) Compactify G-PB arrays (sort and eliminate repeating)
		// ********************************************************************
		for (const auto & ghostIndexPair : gridstorage_.ghostIndexMap_[codim])
		{
			LocalIndexType thisEntityGhostLocalIndex = ghostIndexPair.second;

			//! Compactify the neighbor ranks
			//! \note no need to set-divide here, since only BI and PB were communicated
			VectorHelper::compactify(gridstorage_.G2BIPBNeighborRank_[codim][thisEntityGhostLocalIndex]);
		}

		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: -- Finished Ghost-ProcessBoundary communication construction");
    }


    /** \brief For each PB entity that has non-zero PB-G, and is lowest rank among all its neighbor PB,
     *   communicate to each its neighbour G all other G on the list
     *
     * [TODO]  Implement for periodic, whenever becomes necessary
     * [TODO]  Seems like this procedure could be compactified in 1 loop since there is no feedback
     *
     * */
    void communicateGG(int codim, MPI_Comm comm)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: -- Started Ghost-Ghost communication construction");

    	// 1) Communicate number of Ghost entities this process will communicate to each
    	// It will communicate to each ghost the rest of the ghosts neighboring this entity
    	// It will only communicate if this process owns this entity
    	// ********************************************************************************

    	std::vector<int> nGGEntitySend(size_);
    	std::vector<int> nGGRanksPerProcessSend(size_);

    	// Loop over all PB
		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			//! Assume that all ghost neighbors of this PB are on the other processes
			int nGhostNeighbors = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int candidateRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][0];

			assert(candidateRank != rank_);  // This map should not point to self

			bool hasComm = true;
			hasComm &= (nGhostNeighbors > 1);    // Only communicate to Ghost if at least 2 ghosts share this
			hasComm &= (candidateRank > rank_);  // Only communicate if you own this entity

			if (hasComm)
			{
				int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
				for (int i = 0; i < PBGSize; i++)
				{
					int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
					nGGEntitySend[thisGNeighborRank]++;
					nGGRanksPerProcessSend[thisGNeighborRank] += PBGSize - 1;
				}
			}
		}
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of ghost entities");


    	// 2) Communicate global indices of elements, as well as
		// number of ranks to communicate for each element
    	// ********************************************************************************

		// Construct displacements
		std::vector<int> displSendTmp2(size_);
		int sendSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSendTmp2[i] = (i == 0) ? 0 : displSendTmp2[i-1] + nGGEntitySend[i-1];
			sendSize += nGGEntitySend[i];
		}


		// Fill candidate number send arrays
    	std::vector<int> nGGRankPerEntitySend(sendSize);

    	std::vector<int> GGEntityGlobalIndexSend(sendSize);
		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalIndex = pbIndexPair.first;
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			int nGhostNeighbors = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int candidateRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][0];

			bool hasComm = true;
			hasComm &= (nGhostNeighbors > 1);    // Only communicate to Ghost if at least 2 ghosts share this
			hasComm &= (candidateRank > rank_);  // Only communicate if you own this entity

			if (hasComm)
			{
				int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
    			for (int i = 0; i < PBGSize; i++)
    			{
    				int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
    				int thisTmpIndex = displSendTmp2[thisGNeighborRank]++;

    				// Check if element with this local index exists at all, otherwise bug in the map
    				if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, GGEntityGlobalIndexSend[thisTmpIndex])) {
    					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found by its local index");
    				}

    				// To each ghost neighbour we are going to communicate all other ghost neighbour ranks
    				nGGRankPerEntitySend[thisTmpIndex] = nGhostNeighbors - 1;
    			}
			}
		}


		// Communicate PBG entity global indices
//		MPI_Alltoallv (
//				                   GGEntityGlobalIndexSend.data(),  nGGEntitySend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(GGEntityGlobalIndexRecv.data()), nGGEntityRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);

    	std::vector<int> nGGEntityRecv1;
    	std::vector<int> GGEntityGlobalIndexRecv;

		allcomm_.all2all(GGEntityGlobalIndexSend, nGGEntitySend, GGEntityGlobalIndexRecv, nGGEntityRecv1);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated global indices");


		// Communicate candidate rank numbers
//		MPI_Alltoallv (
//				                   nGGRankPerEntitySend.data(),  nGGEntitySend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(nGGRankPerEntityRecv.data()), nGGEntityRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);

		std::vector<int> nGGEntityRecv2;
		std::vector<int> nGGRankPerEntityRecv;

		allcomm_.all2all(nGGRankPerEntitySend, nGGEntitySend, nGGRankPerEntityRecv, nGGEntityRecv2);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated data sizes");



		//3) Communicate G-G candidate ranks
		// ************************************************************************************

		int sendSize3 = 0;
		std::vector<int> displSendTmp3(size_);

		for (int i = 0; i < size_; i++)
		{
			displSendTmp3[i] = (i == 0) ? 0 : displSendTmp3[i-1] + nGGRanksPerProcessSend[i-1];
			sendSize3 += nGGRanksPerProcessSend[i];
		}


		// Fill in communication arrays
		std::vector<int> neighborGGRankSetSend(sendSize3);
		for (const auto & pbIndexPair : gridstorage_.processBoundaryIndexMap_[codim])
		{
			LocalIndexType thisEntityLocalPBIndex = pbIndexPair.second;

			int nGhostNeighbors = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int candidateRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][0];

			bool hasComm = true;
			hasComm &= (nGhostNeighbors > 1);    // Only communicate to Ghost if at least 2 ghosts share this
			hasComm &= (candidateRank > rank_);  // Only communicate if you own this entity


			if (hasComm)
			{
				// Send to all Ghost neighbors of this PB the ranks of all Ghosts except itself
				int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
				for (int i = 0; i < PBGSize; i++)
				{
					int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
					for (int j = 0; j < PBGSize; j++)
					{
						int thisGNeighborRankSend = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][j];
						if (thisGNeighborRank != thisGNeighborRankSend)
						{
							neighborGGRankSetSend[displSendTmp3[thisGNeighborRank]++] = thisGNeighborRankSend;
						}
					}
				}
			}
		}


//		MPI_Alltoallv (
//				                   neighborGGRankSetSend.data(),  nGGRanksPerProcessSend.data(), displSend.data(), MPI_INT,
//			reinterpret_cast<int*>(neighborGGRankSetRecv.data()), nGGRanksPerProcessRecv.data(), displRecv.data(), MPI_INT,
//			comm
//		);

		std::vector<int> nGGRanksPerProcessRecv;
		std::vector<int> neighborGGRankSetRecv;

		allcomm_.all2all(neighborGGRankSetSend, nGGRanksPerProcessSend, neighborGGRankSetRecv, nGGRanksPerProcessRecv);
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated ranks");


		//4) Fill in
		// ************************************************************************************
		int iRankData = 0;
		int iEntityData = 0;

		for (int i = 0; i < size_; i++)
		{
			for (int j = 0; j < nGGEntityRecv1[i]; j++)
			{
				int nRankPerEntity = nGGRankPerEntityRecv[iEntityData];
				GlobalIndexType thisEntityGlobalIndex = GGEntityGlobalIndexRecv[iEntityData];
				iEntityData++;

				// Get local index corresponding to the communicated global index, check that it exists
				LocalIndexType thisEntityLocalIndex;
				if (!gridbase_.findEntityLocalIndex(codim, thisEntityGlobalIndex, thisEntityLocalIndex)) {
					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found corresponding to communicated global index");
				}

				// Check if the structural type of the received entity is Ghost
				StructuralType thisEntityType = gridbase_.entityPartitionType(codim, thisEntityLocalIndex);
				assert(thisEntityType == Dune::PartitionType::GhostEntity);

				LocalIndexType thisEntityGhostLocalIndex = gridstorage_.ghostIndexMap_[codim][thisEntityLocalIndex];

				for (int k = 0; k < nRankPerEntity; k++)
				{
					int thisNeighborRank = neighborGGRankSetRecv[iRankData++];
					assert((thisNeighborRank >= 0) && (thisNeighborRank < size_) );
					gridstorage_.G2GNeighborRank_[codim][thisEntityGhostLocalIndex].push_back(thisNeighborRank);
				}
			}
		}
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Filled in the received data");


		//5) Compactify G-G arrays (sort and eliminate repeating)
		// ************************************************************************************
		for (const auto & ghostIndexPair : gridstorage_.ghostIndexMap_[codim])
		{
			LocalIndexType thisEntityGhostLocalIndex = ghostIndexPair.second;

			// Compactify the neighbor ranks
			VectorHelper::compactify(gridstorage_.G2GNeighborRank_[codim][thisEntityGhostLocalIndex]);

			//! No need to perform division, as it is assumed that only ghost entities were communicated,
			//! if all the previous steps were done correctly
		}


		// For debugging purposes
		std::vector<int> nNeighborG2G;
		for (int iGhost = 0; iGhost < gridstorage_.G2GNeighborRank_[codim].size(); iGhost++)  { nNeighborG2G.push_back(gridstorage_.G2GNeighborRank_[codim][iGhost].size()); }
		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: --   Number of neighbors for ghost entities=" + VectorHelper::vector2string(nNeighborG2G));

		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearPostConstructor: -- Finished Ghost-Ghost communication construction");
    }





private: // Private members

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

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARPOSTCONSTRUCTOR_HH
