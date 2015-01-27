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

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>


namespace Dune {


template <class ct, int cdim>
class CurvilinearPostConstructor {
public:

    /* public types */
    typedef Dune::CurvilinearGridStorage<ct, cdim>        GridStorageType;
    typedef Dune::CurvilinearGridBase<ct, cdim>           GridBaseType;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
    typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    // Containers
    typedef typename GridStorageType::Vertex                    Vertex;
    typedef typename GridStorageType::VertexStorage             VertexStorage;
    typedef typename GridStorageType::EdgeStorage               EdgeStorage;
    typedef typename GridStorageType::FaceStorage               FaceStorage;
    typedef typename GridStorageType::EntityStorage             EntityStorage;

    // Local and global maps and their iterators
    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;

    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;

    typedef typename GridStorageType::EntityNeighborRankVector  EntityNeighborRankVector;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    // Face Structural Type
    static const unsigned int DomainBoundaryType   = GridStorageType::PartitionType::DomainBoundary;
    static const unsigned int ProcessBoundaryType  = GridStorageType::PartitionType::ProcessBoundary;
    static const unsigned int InternalType         = GridStorageType::PartitionType::Internal;
    static const unsigned int GhostType            = GridStorageType::PartitionType::Ghost;

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV       = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG  = Dune::LoggingMessage::Category::DEBUG;
    static const unsigned int LOG_CATEGORY_ERROR  = Dune::LoggingMessage::Category::ERROR;

public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearPostConstructor(
    		bool verbose,
    		bool processVerbose,
    		GridStorageType & gridstorage,
    		GridBaseType & gridbase,
    		MPIHelper &mpihelper ) :
        verbose_(verbose),
        processVerbose_(processVerbose),
        gridstorage_(gridstorage),
        gridbase_(gridbase),
        mpihelper_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "Initialized CurvilinearPostConstructor");
    }


    // Create iterator lists for quick iteration over mesh entities of a given type
    // NOTE: Must only insert corners, not other interpolatory vertices
    void generateIteratorSets()
    {
        fillPartitionIteratorCorner();

        for (LocalIndexType iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)     { fillPartitionIterator(EDGE_CODIM, iEdge, gridstorage_.edge_[iEdge].structuralType); }

        for (LocalIndexType iFace = 0; iFace < gridstorage_.face_.size(); iFace++)     { fillPartitionIterator(FACE_CODIM, iFace, gridstorage_.face_[iFace].structuralType); }

        for (LocalIndexType iElem = 0; iElem < gridstorage_.element_.size(); iElem++)  { fillPartitionIterator(ELEMENT_CODIM, iElem, gridstorage_.element_[iElem].structuralType); }

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
		IndexSetIterator iterB = gridbase_.entityIndexBegin(FACE_CODIM , ProcessBoundaryType);
		IndexSetIterator iterE = gridbase_.entityIndexEnd(FACE_CODIM , ProcessBoundaryType);

		int nEdgeTriangle = 3;
		int nVertexTriangle = 3;
		int nTriangleTet = 4;
		int nEdgeTet = 6;
		int nVertexTet = 4;

		StructuralType tmpTypes = { InternalType, GhostType };

		// Resize the PB->G comm map, since we already know its size, it is the number of PB entities
		for (int iCodim = 1; iCodim <= cdim; iCodim++)
		{
			gridstorage_.PB2GNeighborRank_[iCodim].resize(gridstorage_.PB2PBNeighborRank_[iCodim].size());
		}




		// Iterate over all PB faces
		for (IndexSetIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisFaceLocalIndex = *iter;
			int thisFaceNeighborRank = gridbase_.processBoundaryNeighborRankSet(FACE_CODIM, thisFaceLocalIndex)[0];

			// 1) Find all Internal and Ghost entities associated with this PB Face
			// **********************************************************************

			// Fill in the sets associated with this face
			std::vector<LocalIndexType> faceSubentityIndex[4];
			faceSubentityIndex[1].push_back(thisFaceLocalIndex);
			for (int i = 0; i < nEdgeTriangle; i++)    { faceSubentityIndex[2].push_back(gridbase_.subentityLocalIndex(thisFaceLocalIndex, FACE_CODIM, EDGE_CODIM, i)); }
			for (int i = 0; i < nVertexTriangle; i++)  { faceSubentityIndex[3].push_back(gridbase_.subentityLocalIndex(thisFaceLocalIndex, FACE_CODIM, VERTEX_CODIM, i)); }

			// Fill in the sets associated with internal element
			std::vector<LocalIndexType> faceNeighborSubentityIndex[2][4];
			for (int iTmp = 0; iTmp <= 1; iTmp++)
			{
				LocalIndexType thisElementLocalIndex = gridbase_.faceNeighbor(thisFaceLocalIndex, iTmp);

				faceNeighborSubentityIndex[iTmp][0].push_back(thisElementLocalIndex);
				for (int i = 0; i < nTriangleTet; i++)  { faceNeighborSubentityIndex[iTmp][1].push_back(gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, FACE_CODIM, i)); }
				for (int i = 0; i < nEdgeTet; i++)      { faceNeighborSubentityIndex[iTmp][2].push_back(gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, EDGE_CODIM, i)); }
				for (int i = 0; i < nVertexTet; i++)    { faceNeighborSubentityIndex[iTmp][3].push_back(gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, VERTEX_CODIM, i)); }


				for (int iCodim = 1; iCodim <= cdim; iCodim++)
				{
					// Subtract face subentity set from element subentity set, such that only internal and ghost subentities are left
					faceNeighborSubentityIndex[iTmp][iCodim] = Dune::VectorHelper::sortedSetComplement(faceNeighborSubentityIndex[iTmp][iCodim], faceSubentityIndex[iCodim]);

					// 2) Add entities to the map unless they are already added.
					// **********************************************************************
					Local2LocalMap & thisLocalMap = selectCommMap(iCodim, tmpTypes[iTmp]);

					for (int iEntity = 0; iEntity < faceNeighborSubentityIndex[iTmp][iCodim].size(); iEntity++ )
					{
						LocalIndexType thisEntityLocalIndex = faceNeighborSubentityIndex[iTmp][iCodim][iEntity];
						StructuralType thisEntityType = gridbase_.entityStructuralType(iCodim, thisEntityLocalIndex);


						// If this is a PB Entity, that happens to not be on the face, it is possible
						// that it is on a face of a different process, then it is a PB->G link
						// Otherwise, this entity is either internal or ghost, and it contributest to one
						// of the new maps
						if (thisEntityType == ProcessBoundaryType)
						{
							LocalIndexType thisEntityPBIndex = gridstorage_.processBoundaryIndexMap_[iCodim][thisEntityLocalIndex];
							gridstorage_.PB2GNeighborRank_[iCodim][thisEntityPBIndex].push_back(thisFaceNeighborRank);
						}
						else
						{
							Local2LocalIterator thisIter = thisLocalMap.find(thisEntityLocalIndex);
							LocalIndexType thisEntitySubsetIndex;

							if (thisIter == thisLocalMap.end())
							{
								thisEntitySubsetIndex = thisLocalMap.size();
								thisLocalMap[thisEntityLocalIndex] = thisEntitySubsetIndex;

								// Enlarge the neighbour rank vector such that it has enough entries for all entities of subset
								if (iTmp == 0)  { gridstorage_.BI2GNeighborRank_[iCodim].push_back(std::vector<int>()); }
								else            { gridstorage_.G2BIPBNeighborRank_[iCodim].push_back(std::vector<int>()); }

							} else {
								thisEntitySubsetIndex = (*thisIter).second;
							}

							// 3) Mark PB face neighbor rank on all subset entities
							// **********************************************************************
							if (iTmp == 0)  { gridstorage_.BI2GNeighborRank_[iCodim][thisEntitySubsetIndex].push_back(thisFaceNeighborRank); }
							else            { gridstorage_.G2BIPBNeighborRank_[iCodim][thisEntitySubsetIndex].push_back(thisFaceNeighborRank); }
						}
					}
				}
			}
		}
    }


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
    void communicateCommunicationEntityNeighborRanks()
    {

    }




protected:

    Local2LocalMap & selectCommMap(int codim, StructuralType structtype)
    {
    	switch (structtype)
    	{
    	case InternalType          :  return gridstorage_.boundaryInternalEntityIndexMap_[codim];  break;
    	case ProcessBoundaryType   :  return gridstorage_.processBoundaryIndexMap_[codim];         break;
    	case GhostType             :  return gridstorage_.ghostIndexMap_[codim];                   break;

    	}
    }


    EntityNeighborRankVector & selectCommRankVector(int codim, StructuralType structSend, StructuralType structRecv)
    {
    	// Can only communicate over these 3 PartitionTypes
    	assert((structRecv == InternalType)||(structRecv == ProcessBoundaryType)||(structRecv == GhostType));

    	switch (structSend)
    	{
    	case InternalType          :   // Internal -> Ghost protocol
    	{
    		assert(structRecv == GhostType);
    		return gridstorage_.BI2GNeighborRank_[codim];
    	} break;
    	case ProcessBoundaryType   :   // PB -> PB and PB -> Ghost protocols
    	{
    		assert((structRecv == ProcessBoundaryType)||(structRecv == GhostType));
    		if (structRecv == ProcessBoundaryType)  { return gridstorage_.PB2PBNeighborRank_[codim]; }
    		if (structRecv == GhostType)            { return gridstorage_.PB2GNeighborRank_[codim]; }
    	} break;
    	case GhostType             :   // Ghost -> (Internal & PB) and Ghost -> Ghost protocols
    	{
    		assert((structRecv == InternalType)||(structRecv == ProcessBoundaryType)||(structRecv == GhostType));
    		if (structRecv == InternalType)         { return gridstorage_.G2BIPBNeighborRank_[codim]; }
    		if (structRecv == ProcessBoundaryType)  { return gridstorage_.G2BIPBNeighborRank_[codim]; }
    		if (structRecv == GhostType)            { return gridstorage_.G2GNeighborRank_[codim]; }
    	} break;

    	}
    }


    // Loop over all corners of the mesh (including ghost) and add them to the interator set
    // Do this by iterating over all corners of all elements
    void fillPartitionIteratorCorner()
    {
    	// Loop over all elements
    	for (LocalIndexType iElem = 0; iElem < gridstorage_.element_.size(); iElem++)
    	{
        	// Get LocalCornerIndices
    		std::vector<LocalIndexType> thisCornerIndex =  gridbase_.entityCornerLocalIndex(VERTEX_CODIM, iElem);

        	// Loop over LocalCornerIndices
    		for (int iCorner = 0; iCorner < thisCornerIndex.size(); iCorner++)
    		{
    			// If this corner has not been added yet, add it
    			LocalIndexType thisIndex = thisCornerIndex[iCorner];
    			if (gridstorage_.entityAllIndexSet_[VERTEX_CODIM].find(thisIndex) == gridstorage_.entityAllIndexSet_[VERTEX_CODIM].end())
    			{
    				fillPartitionIterator (VERTEX_CODIM, thisIndex, gridstorage_.point_[thisIndex].structuralType);
    			}
    		}
    	}
    }


    // Fill all LocalIndexSets necessary to iterate over the grid
    // Note: Must not contain all interpolatory vertices, only corners
    void fillPartitionIterator (int codim, LocalIndexType localIndex, StructuralType structtype)
    {
    	// Check if the entity is of a valid type
    	gridbase_.assertValidCodimStructuralType(codim, structtype);

    	// All-set includes entities of all structural types
        gridstorage_.entityAllIndexSet_[codim].insert(localIndex);

    	switch (structtype)
    	{
    	case InternalType          :
    	{
    		gridstorage_.entityInternalIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorBorderIndexSet_[codim].insert(localIndex);
    	} break;
    	case ProcessBoundaryType   :
    	{
    		gridstorage_.entityProcessBoundaryIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorBorderIndexSet_[codim].insert(localIndex);
    	} break;
    	case DomainBoundaryType    :
    	{
    		gridstorage_.entityDomainBoundaryIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorIndexSet_[codim].insert(localIndex);
    		gridstorage_.entityDuneInteriorBorderIndexSet_[codim].insert(localIndex);
    		break;
    	}
    	case GhostType             :
    	{
    		gridstorage_.entityGhostIndexSet_[codim].insert(localIndex);
    	} break;

    	}
    }




private: // Private members

    bool verbose_;
    bool processVerbose_;

    // Curvilinear Grid Storage Class
    GridStorageType & gridstorage_;

    // Reference to Curvilinear Grid Base - necessary for OCTree construction
    GridBaseType & gridbase_;


    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace Dune

#endif  // DUNE_CURVILINEARPOSTCONSTRUCTOR_HH
