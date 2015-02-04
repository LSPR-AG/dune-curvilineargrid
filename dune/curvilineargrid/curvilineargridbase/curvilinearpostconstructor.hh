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
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Started generating iterator lists");

        fillPartitionIteratorCorner();

        for (LocalIndexType iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)     { fillPartitionIterator(EDGE_CODIM, iEdge, gridstorage_.edge_[iEdge].structuralType); }

        for (LocalIndexType iFace = 0; iFace < gridstorage_.face_.size(); iFace++)     { fillPartitionIterator(FACE_CODIM, iFace, gridstorage_.face_[iFace].structuralType); }

        for (LocalIndexType iElem = 0; iElem < gridstorage_.element_.size(); iElem++)  { fillPartitionIterator(ELEMENT_CODIM, iElem, gridstorage_.element_[iElem].structuralType); }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Finished generating iterator lists");

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
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Started generating communication maps");

		IndexSetIterator iterB = gridbase_.entityIndexBegin(FACE_CODIM , ProcessBoundaryType);
		IndexSetIterator iterE = gridbase_.entityIndexEnd(FACE_CODIM , ProcessBoundaryType);

		int nEdgeTriangle = 3;
		int nVertexTriangle = 3;
		int nTriangleTet = 4;
		int nEdgeTet = 6;
		int nVertexTet = 4;

		StructuralType tmpTypes[2] = { InternalType, GhostType };

		// Resize all neighbor rank arrays, except the boundary internals, since we do not know their lengths
		// Boundary Internals will be resized as its number is calculated
		for (int iCodim = 0; iCodim <= cdim; iCodim++)
		{
			int nPB = (iCodim == 0) ? 0 : gridbase_.nEntity(iCodim, ProcessBoundaryType);
			int nG  = gridbase_.nEntity(iCodim, GhostType);

			// PB2PB already exists, so need not resize it
			gridstorage_.PB2GNeighborRank_[iCodim].resize(nPB);
			gridstorage_.G2BIPBNeighborRank_[iCodim].resize(nG);
			gridstorage_.G2GNeighborRank_[iCodim].resize(nG);
		}




		// Iterate over all PB faces
		for (IndexSetIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisFaceLocalIndex = *iter;
			int thisFaceNeighborRank = gridbase_.commEntityNeighborRankSet(FACE_CODIM, thisFaceLocalIndex, ProcessBoundaryType, ProcessBoundaryType)[0];

			// 1) Find all Internal and Ghost entities associated with this PB Face
			// **********************************************************************

			// Fill in the sets associated with this face
			std::vector<LocalIndexType> faceSubentityIndex[4];
			faceSubentityIndex[1].push_back(thisFaceLocalIndex);
			for (int i = 0; i < nEdgeTriangle; i++)    { faceSubentityIndex[EDGE_CODIM].push_back  (gridbase_.subentityLocalIndex(thisFaceLocalIndex, FACE_CODIM, EDGE_CODIM, i)); }
			for (int i = 0; i < nVertexTriangle; i++)  { faceSubentityIndex[VERTEX_CODIM].push_back(gridbase_.subentityLocalIndex(thisFaceLocalIndex, FACE_CODIM, VERTEX_CODIM, i)); }

			// Fill in the sets associated with internal element
			std::vector<LocalIndexType> faceNeighborSubentityIndex[2][4];
			for (int iTmp = 0; iTmp <= 1; iTmp++)
			{
				// Gets either internal or ghost element depending on iTmp
				LocalIndexType thisElementLocalIndex = gridbase_.faceNeighbor(thisFaceLocalIndex, iTmp);

				faceNeighborSubentityIndex[iTmp][ELEMENT_CODIM].push_back(thisElementLocalIndex);
				for (int i = 0; i < nTriangleTet; i++)  { faceNeighborSubentityIndex[iTmp][FACE_CODIM].push_back  (gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, FACE_CODIM, i)); }
				for (int i = 0; i < nEdgeTet; i++)      { faceNeighborSubentityIndex[iTmp][EDGE_CODIM].push_back  (gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, EDGE_CODIM, i)); }
				for (int i = 0; i < nVertexTet; i++)    { faceNeighborSubentityIndex[iTmp][VERTEX_CODIM].push_back(gridbase_.subentityLocalIndex(thisElementLocalIndex, ELEMENT_CODIM, VERTEX_CODIM, i)); }


				for (int iCodim = 0; iCodim <= cdim; iCodim++)
				{
					// Sort subentity vectors
					std::sort(faceSubentityIndex[iCodim].begin(), faceSubentityIndex[iCodim].end());
					std::sort(faceNeighborSubentityIndex[iTmp][iCodim].begin(), faceNeighborSubentityIndex[iTmp][iCodim].end());

					// Subtract face subentity set from element subentity set, such that only internal and ghost subentities are left
					faceNeighborSubentityIndex[iTmp][iCodim] = Dune::VectorHelper::sortedSetComplement(faceNeighborSubentityIndex[iTmp][iCodim], faceSubentityIndex[iCodim]);

					// 2) Add entities to the map unless they are already added.
					// **********************************************************************
					Local2LocalMap & thisLocalMap = gridbase_.selectCommMap(iCodim, tmpTypes[iTmp]);

					for (int iEntity = 0; iEntity < faceNeighborSubentityIndex[iTmp][iCodim].size(); iEntity++ )
					{
						LocalIndexType thisEntityLocalIndex = faceNeighborSubentityIndex[iTmp][iCodim][iEntity];
						StructuralType thisEntityType = gridbase_.entityStructuralType(iCodim, thisEntityLocalIndex);

						std::stringstream log_str;
						log_str << "CurvilinearPostConstructor: ---- Iterating over iTmp=" << gridstorage_.PartitonTypeName[tmpTypes[iTmp]];
						log_str << " codim=" << iCodim;
						log_str << " subentityNo=" << iEntity;
						log_str << " localindex =" << thisEntityLocalIndex;
						log_str << " gives structural type =" << gridstorage_.PartitonTypeName[thisEntityType];
						Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_str.str());


						// If this is a PB Entity, that happens to not be on the face, it is possible
						// that it is on a face of a different process, then it is a PB->G link
						// Otherwise, this entity is either internal or ghost, and it contributest to one
						// of the new maps
						if (thisEntityType == ProcessBoundaryType)
						{
							std::vector<int> & thisEntityPBNeighbors = gridbase_.commEntityNeighborRankSet(iCodim, thisEntityLocalIndex, ProcessBoundaryType, ProcessBoundaryType);
							bool isNewRank = !Dune::VectorHelper::isInside(thisEntityPBNeighbors, thisFaceNeighborRank);
							//std::cout << "testing isInside vector=(" << Dune::VectorHelper::vector2string(thisEntityPBNeighbors) << ") elem=" << thisFaceNeighborRank << " isNew=" << isNewRank << std::endl;

							// If this rank is not already in PB-PB of this entity, then it must be in PB-G
							if (isNewRank)
							{
								LocalIndexType thisEntityPBIndex = gridstorage_.processBoundaryIndexMap_[iCodim][thisEntityLocalIndex];
								gridstorage_.PB2GNeighborRank_[iCodim][thisEntityPBIndex].push_back(thisFaceNeighborRank);
							}
						}
						else
						{
							// Check if the type of the entity is expected
					    	bool entityTypeExpected =
					    		((iTmp == 0)&&((thisEntityType == InternalType)||(thisEntityType == DomainBoundaryType))) ||
					    		((iTmp == 1)&&(thisEntityType == GhostType));

					    	if (!entityTypeExpected) {
					    		GlobalIndexType thisEntityGlobalIndex;

					    		if (!gridbase_.findEntityGlobalIndex(iCodim, thisEntityLocalIndex, thisEntityGlobalIndex)) {
					    			std::cout << "Global index not found " << std::endl;
					    		}

					    		std::cout << "error in codim " << iCodim << " entity localindex=" << thisEntityLocalIndex << " globalIndex=" << thisEntityGlobalIndex << std::endl;

					    		std::string expectedTypeName = gridstorage_.PartitonTypeName[tmpTypes[iTmp]];
					    		std::string receivedTypeName = gridstorage_.PartitonTypeName[thisEntityType];

					    		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Unexpected type name expected=" + expectedTypeName + ", received="+receivedTypeName);
					    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected type name");
					    	}


							Local2LocalIterator thisIter = thisLocalMap.find(thisEntityLocalIndex);
							LocalIndexType thisEntitySubsetIndex;

							if (thisIter == thisLocalMap.end())
							{
								thisEntitySubsetIndex = thisLocalMap.size();
								thisLocalMap[thisEntityLocalIndex] = thisEntitySubsetIndex;

								// Enlarge the neighbour rank vector such that it has enough entries for all entities of subset
								if (iTmp == 0)  { gridstorage_.BI2GNeighborRank_[iCodim].resize(thisEntitySubsetIndex + 1); }

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

		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Finished generating communication maps");
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
    void communicateCommunicationEntityNeighborRanks()
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();

    	for (int iCodim = 0; iCodim <= cdim; iCodim++)
    	{
    		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Started communicating entity ranks for codim=" + std::to_string(iCodim));


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

    	    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Finished communicating entity ranks for codim=" + std::to_string(iCodim));
    	}

    }




protected:





    // ************************************************************************************
	// Auxiliary methods for Iterator list generation
	// ************************************************************************************

    // Loop over all corners of the mesh (including ghost) and add them to the interator set
    // Do this by iterating over all corners of all elements
    void fillPartitionIteratorCorner()
    {
    	// Loop over all elements
    	for (LocalIndexType iElem = 0; iElem < gridstorage_.element_.size(); iElem++)
    	{
        	// Get LocalCornerIndices
    		std::vector<LocalIndexType> thisCornerIndex =  gridbase_.entityCornerLocalIndex(ELEMENT_CODIM, iElem);

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




    // ************************************************************************************
	// Auxiliary methods for communicating communication entity neighbor ranks
	// ************************************************************************************


    // [TODO] Possibly set complement unnecessary, since already only adding neighbor if it is not a neighbor already
    void communicatePBG(int codim, MPI_Comm comm)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: -- Started ProcessBoundary-Ghost communication construction");

    	Local2LocalIterator iterB = gridstorage_.processBoundaryIndexMap_[codim].begin();
    	Local2LocalIterator iterE = gridstorage_.processBoundaryIndexMap_[codim].end();

		//1) Iterate over all PB entities, find entities with non-zero PB-G candidates
		//   Mark these entities to send to all their PB neighbors
		// ************************************************************************************

    	std::vector<int> nPBGEntitySend(size_);
    	std::vector<int> nPBGEntityRecv(size_);
    	std::vector<int> nPBGRankPerProcessSend(size_);
    	std::vector<int> nPBGRankPerProcessRecv(size_);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			// PB2G may contain repeating entities, and is not sorted, need to compactify
			Dune::VectorHelper::compactify(gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex]);

			//1.1) divide provisional PB->G set by PB->PB set to see which provisional PB->G are new
			gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex] = Dune::VectorHelper::sortedSetComplement(
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

		// Communicate entity number and candidate number per process
		MPI_Alltoall(nPBGEntitySend.data(), 1, MPI_INT, reinterpret_cast<int*>(nPBGEntityRecv.data()), 1, MPI_INT, comm);
		MPI_Alltoall(nPBGRankPerProcessSend.data(), 1, MPI_INT, reinterpret_cast<int*>(nPBGRankPerProcessRecv.data()), 1, MPI_INT, comm);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of ProcessBoundaries with candidates");


		//2) For each PB-G entity, communicate its global index
		// as well as communicate how many candidates it is going to send
		// ************************************************************************************

		// Construct displacements
		std::vector<int> displSend(size_);
		std::vector<int> displRecv(size_);
		std::vector<int> displSendTmp(size_);
		int sendSize = 0;
		int recvSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSend[i] = (i == 0) ? 0 : displSend[i-1] + nPBGEntitySend[i-1];
			displRecv[i] = (i == 0) ? 0 : displRecv[i-1] + nPBGEntityRecv[i-1];
			sendSize += nPBGEntitySend[i];
			recvSize += nPBGEntityRecv[i];

		}
		displSendTmp = displSend;


		// Fill candidate number send arrays
    	std::vector<int> nPBGRankPerEntitySend(sendSize);
    	std::vector<int> nPBGRankPerEntityRecv(recvSize);
    	std::vector<int> PBGEntityGlobalIndexSend(sendSize);
    	std::vector<int> PBGEntityGlobalIndexRecv(recvSize);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalIndex = (*iter).first;
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			// 1.2) Mark number of real PB-G candidates for each process
			int PBPBSize = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int realPBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();

			if (realPBGSize > 0)
			{
    			for (int i = 0; i < PBPBSize; i++)
    			{
    				int thisPBNeighborRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][i];
    				int thisTmpIndex = displSendTmp[thisPBNeighborRank]++;

    				if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, PBGEntityGlobalIndexSend[thisTmpIndex]))
    				{
    					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found by its local index");
    				}
    				nPBGRankPerEntitySend[thisTmpIndex] = realPBGSize;
    			}
			}
		}


		// Communicate PBG entity global indices
		MPI_Alltoallv (
			                       PBGEntityGlobalIndexSend.data(),  nPBGEntitySend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(PBGEntityGlobalIndexRecv.data()), nPBGEntityRecv.data(), displRecv.data(), MPI_INT,
			comm
		);

		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated global indices");


		// Communicate candidate rank numbers
		MPI_Alltoallv (
				                   nPBGRankPerEntitySend.data(),  nPBGEntitySend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(nPBGRankPerEntityRecv.data()), nPBGEntityRecv.data(), displRecv.data(), MPI_INT,
			comm
		);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of ranks per entity");



		//3) Communicate PB-G candidate ranks
		// ************************************************************************************

		sendSize = 0;
		recvSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSend[i] = (i == 0) ? 0 : displSend[i-1] + nPBGRankPerProcessSend[i-1];
			displRecv[i] = (i == 0) ? 0 : displRecv[i-1] + nPBGRankPerProcessRecv[i-1];
			sendSize += nPBGRankPerProcessSend[i];
			recvSize += nPBGRankPerProcessRecv[i];
		}
		displSendTmp = displSend;

		std::cout << "process_" << rank_ << " a with size = " << sendSize << std::endl;

		// Fill in communication arrays
		std::vector<int> neighborPBGRankSetSend(sendSize);
		std::vector<int> neighborPBGRankSetRecv(recvSize);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			// Send all provisional PB-G ranks to all PB neighbors of this entity
			int PBPBSize = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int realPBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();

			for (int i = 0; i < PBPBSize; i++)
			{
    			for (int j = 0; j < realPBGSize; j++)
    			{
    				int thisPBPBRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][i];
    				int thisPBGRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
    				int tmpIndex = displSendTmp[thisPBPBRank]++;

    				std::cout << "process_" << rank_ << " e2 tmpindex=" << tmpIndex << std::endl;

    				neighborPBGRankSetSend[tmpIndex] = thisPBGRank;
    			}
			}
		}

		MPI_Alltoallv (
			                       neighborPBGRankSetSend.data(),  nPBGRankPerProcessSend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(neighborPBGRankSetRecv.data()), nPBGRankPerProcessRecv.data(), displRecv.data(), MPI_INT,
			comm
		);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated ranks");


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
				StructuralType thisEntityType = gridbase_.entityStructuralType(codim, thisEntityLocalIndex);
				assert(thisEntityType == ProcessBoundaryType);

				LocalIndexType thisEntityLocalPBIndex = gridstorage_.processBoundaryIndexMap_[codim][thisEntityLocalIndex];
				for (int k = 0; k < nRankPerEntity; k++)
				{
					int thisNeighborRank = neighborPBGRankSetRecv[iRankData++];
					gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].push_back(thisNeighborRank);

				}

			}
		}
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Filled in received data");


		//5) Compactify PB-G arrays (sort and eliminate repeating)
		// ************************************************************************************

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			// Compactify the neighbor ranks
			Dune::VectorHelper::compactify(gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex]);

			// divide PB->G set by PB->PB set
			// This is just to make sure that no PB-PB links were picked up in the process
			gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex] = Dune::VectorHelper::sortedSetComplement(
				gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex],
				gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex]
			);
		}

		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: -- Finished ProcessBoundary-Ghost communication construction");
    }


    /** \brief For each PB entity that has non-zero PB-G, communicate its globalIndex and rank of self
     *
     * Algorithm:
     * 1) Communicate number of G-PB entities to send to each process
     * 2) Communicate global indices for each G-PB entity
     * 3) On receiving end, mark sender's rank on all received G-PB
     *
     *
     *
     * */
    void communicateGPB(int codim, MPI_Comm comm)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: -- Started Ghost-ProcessBoundary communication construction");

    	Local2LocalIterator iterB = gridstorage_.processBoundaryIndexMap_[codim].begin();
    	Local2LocalIterator iterE = gridstorage_.processBoundaryIndexMap_[codim].end();

		// 1) Communicate number of G-PB entities to send to each process
		// ********************************************************************
    	std::vector<int> nGPBEntitySend(size_);
    	std::vector<int> nGPBEntityRecv(size_);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			for (int i = 0; i < PBGSize; i++)
			{
				int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
				nGPBEntitySend[thisGNeighborRank]++;
			}
		}
		MPI_Alltoall(nGPBEntitySend.data(), 1, MPI_INT, reinterpret_cast<int*>(nGPBEntityRecv.data()), 1, MPI_INT, comm);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of Ghost candidates per process");


		// 2) Communicate global indices for each G-PB entity
		// ********************************************************************

		// Create displacement arrays
		std::vector<int> displSend(size_);
		std::vector<int> displRecv(size_);
		std::vector<int> displSendTmp(size_);
		int sendSize = 0;
		int recvSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSend[i] = (i == 0) ? 0 : displSend[i-1] + nGPBEntitySend[i-1];
			displRecv[i] = (i == 0) ? 0 : displRecv[i-1] + nGPBEntityRecv[i-1];
			sendSize += nGPBEntitySend[i];
			recvSize += nGPBEntityRecv[i];

		}
		displSendTmp = displSend;

		// fill global index array
    	std::vector<int> GPBEntityGlobalIndexSend(sendSize);
    	std::vector<int> GPBEntityGlobalIndexRecv(recvSize);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalIndex = (*iter).first;
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			int PBGSize = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			for (int i = 0; i < PBGSize; i++)
			{
				int thisGNeighborRank = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex][i];
				int thisTmpIndex = displSendTmp[thisGNeighborRank]++;

				// Check if element with this local index exists at all, otherwise bug in the map
				if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, GPBEntityGlobalIndexSend[thisTmpIndex])) {
					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found by its local index");
				}
			}
		}

		// Communicate
		MPI_Alltoallv (
				                   GPBEntityGlobalIndexSend.data(),  nGPBEntitySend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(GPBEntityGlobalIndexRecv.data()), nGPBEntityRecv.data(), displRecv.data(), MPI_INT,
			comm
		);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Sent own rank to all ghosts");


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
				StructuralType thisEntityType = gridbase_.entityStructuralType(codim, thisEntityLocalIndex);
				if (thisEntityType != GhostType)
				{
					Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   codim" + std::to_string(codim) + "entity globalIndex=" + std::to_string(thisEntityGlobalIndex) + " expected type=" + gridstorage_.PartitonTypeName[GhostType] + " received=" + gridstorage_.PartitonTypeName[thisEntityType]);
					assert(thisEntityType == GhostType);
				}




				LocalIndexType thisEntityLocalGhostIndex = gridstorage_.ghostIndexMap_[codim][thisEntityLocalIndex];
				gridstorage_.G2BIPBNeighborRank_[codim][thisEntityLocalGhostIndex].push_back(i);
			}
		}
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   filled received data");


		// 4) Compactify G-PB arrays (sort and eliminate repeating)
		// ********************************************************************
		Local2LocalIterator iterGB = gridstorage_.ghostIndexMap_[codim].begin();
		Local2LocalIterator iterEB = gridstorage_.ghostIndexMap_[codim].end();

		for (Local2LocalIterator iter = iterGB; iter != iterEB; iter++)
		{
			LocalIndexType thisEntityGhostLocalIndex = (*iter).second;

			//! Compactify the neighbor ranks
			//! \note no need to set-divide here, since only BI and PB were communicated
			Dune::VectorHelper::compactify(gridstorage_.G2BIPBNeighborRank_[codim][thisEntityGhostLocalIndex]);
		}

		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: -- Finished Ghost-ProcessBoundary communication construction");
    }


    /** \brief For each PB entity that has non-zero PB-G, and is lowest rank among all its neighbor PB,
     *   communicate to each its neighbour G all other G on the list
     *
     * [TODO]  Seems like this procedure could be compactified in 1 loop since there is no feedback
     *
     * */
    void communicateGG(int codim, MPI_Comm comm)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: -- Started Ghost-Ghost communication construction");

    	Local2LocalIterator iterB = gridstorage_.processBoundaryIndexMap_[codim].begin();
    	Local2LocalIterator iterE = gridstorage_.processBoundaryIndexMap_[codim].end();


    	// 1) Communicate number of Ghost entities this process will communicate to each
    	// It will communicate to each ghost the rest of the ghosts neighboring this entity
    	// It will only communicate if this process owns this entity
    	// ********************************************************************************

    	std::vector<int> nGGEntitySend(size_);
    	std::vector<int> nGGEntityRecv(size_);
    	std::vector<int> nGGRanksPerProcessSend(size_);
    	std::vector<int> nGGRanksPerProcessRecv(size_);

    	// Loop over all PB
		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

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
		MPI_Alltoall(nGGEntitySend.data(), 1, MPI_INT, reinterpret_cast<int*>(nGGEntityRecv.data()), 1, MPI_INT, comm);
		MPI_Alltoall(nGGRanksPerProcessSend.data(), 1, MPI_INT, reinterpret_cast<int*>(nGGRanksPerProcessRecv.data()), 1, MPI_INT, comm);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated number of ghost entities");


    	// 2) Communicate global indices of elements, as well as
		// number of ranks to communicate for each element
    	// ********************************************************************************

		// Construct displacements
		std::vector<int> displSend(size_);
		std::vector<int> displRecv(size_);
		std::vector<int> displSendTmp(size_);
		int sendSize = 0;
		int recvSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSend[i] = (i == 0) ? 0 : displSend[i-1] + nGGEntitySend[i-1];
			displRecv[i] = (i == 0) ? 0 : displRecv[i-1] + nGGEntityRecv[i-1];
			sendSize += nGGEntitySend[i];
			recvSize += nGGEntityRecv[i];

		}
		displSendTmp = displSend;


		// Fill candidate number send arrays
    	std::vector<int> nGGRankPerEntitySend(sendSize);
    	std::vector<int> nGGRankPerEntityRecv(recvSize);
    	std::vector<int> GGEntityGlobalIndexSend(sendSize);
    	std::vector<int> GGEntityGlobalIndexRecv(recvSize);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalIndex = (*iter).first;
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

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
    				int thisTmpIndex = displSendTmp[thisGNeighborRank]++;

    				// Check if element with this local index exists at all, otherwise bug in the map
    				if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, GGEntityGlobalIndexSend[thisTmpIndex])) {
    					DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Element not found by its local index");
    				}

    				nGGRankPerEntitySend[thisTmpIndex] = nGhostNeighbors - 1;
    			}
			}
		}


		// Communicate PBG entity global indices
		MPI_Alltoallv (
				                   GGEntityGlobalIndexSend.data(),  nGGEntitySend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(GGEntityGlobalIndexRecv.data()), nGGEntityRecv.data(), displRecv.data(), MPI_INT,
			comm
		);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated global indices");


		// Communicate candidate rank numbers
		MPI_Alltoallv (
				                   nGGRankPerEntitySend.data(),  nGGEntitySend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(nGGRankPerEntityRecv.data()), nGGEntityRecv.data(), displRecv.data(), MPI_INT,
			comm
		);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated data sizes");



		//3) Communicate G-G candidate ranks
		// ************************************************************************************

		sendSize = 0;
		recvSize = 0;

		for (int i = 0; i < size_; i++)
		{
			displSend[i] = (i == 0) ? 0 : displSend[i-1] + nGGRanksPerProcessSend[i-1];
			displRecv[i] = (i == 0) ? 0 : displRecv[i-1] + nGGRanksPerProcessRecv[i-1];
			sendSize += nGGRanksPerProcessSend[i];
			recvSize += nGGRanksPerProcessRecv[i];
		}
		displSendTmp = displSend;


		std::cout << "process_" << rank_ << "aa" << std::endl;

		// Fill in communication arrays
		std::vector<int> neighborGGRankSetSend(sendSize);
		std::vector<int> neighborGGRankSetRecv(recvSize);

		for (Local2LocalIterator iter = iterB; iter != iterE; iter++)
		{
			LocalIndexType thisEntityLocalPBIndex = (*iter).second;

			int nGhostNeighbors = gridstorage_.PB2GNeighborRank_[codim][thisEntityLocalPBIndex].size();
			int candidateRank = gridstorage_.PB2PBNeighborRank_[codim][thisEntityLocalPBIndex][0];

			std::cout << "process_" << rank_ << "bb" << std::endl;

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
							std::cout << "process_" << rank_ << "cc" << std::endl;

							neighborGGRankSetSend[displSendTmp[thisGNeighborRank]++] = thisGNeighborRankSend;
						}

					}
				}
			}
		}

		std::cout << "process_" << rank_ << " ggg procsetsend=(" << Dune::VectorHelper::vector2string(nGGRanksPerProcessSend);
		std::cout << " sendlength=" << neighborGGRankSetSend.size();
		std::cout << " recvlength=" << neighborGGRankSetRecv.size() << std::endl;

		MPI_Alltoallv (
				                   neighborGGRankSetSend.data(),  nGGRanksPerProcessSend.data(), displSend.data(), MPI_INT,
			reinterpret_cast<int*>(neighborGGRankSetRecv.data()), nGGRanksPerProcessRecv.data(), displRecv.data(), MPI_INT,
			comm
		);
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Communicated ranks");


		//4) Fill in
		// ************************************************************************************
		int iRankData = 0;
		int iEntityData = 0;

		for (int i = 0; i < size_; i++)
		{
			for (int j = 0; j < nGGEntityRecv[i]; j++)
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
				StructuralType thisEntityType = gridbase_.entityStructuralType(codim, thisEntityLocalIndex);
				assert(thisEntityType == GhostType);

				LocalIndexType thisEntityGhostLocalIndex = gridstorage_.ghostIndexMap_[codim][thisEntityLocalIndex];

				for (int k = 0; k < nRankPerEntity; k++)
				{
					int thisNeighborRank = neighborGGRankSetRecv[iRankData++];
					gridstorage_.G2GNeighborRank_[codim][thisEntityGhostLocalIndex].push_back(thisNeighborRank);
				}
			}
		}
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: --   Filled in the received data");


		//5) Compactify G-G arrays (sort and eliminate repeating)
		// ************************************************************************************

		Local2LocalIterator iterGB = gridstorage_.ghostIndexMap_[codim].begin();
		Local2LocalIterator iterGE = gridstorage_.ghostIndexMap_[codim].end();


		std::cout << "process_" << rank_ << " has total G2G=" << gridstorage_.G2GNeighborRank_[codim].size() << std::endl;

		for (Local2LocalIterator iter = iterGB; iter != iterGE; iter++)
		{
			LocalIndexType thisEntityGhostLocalIndex = (*iter).second;

			// Compactify the neighbor ranks
			Dune::VectorHelper::compactify(gridstorage_.G2GNeighborRank_[codim][thisEntityGhostLocalIndex]);

			//! No need to perform division, as it is assumed that only ghost entities were communicated,
			//! if all the previous steps were done correctly
		}

		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: -- Finished Ghost-Ghost communication construction");
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
