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

#ifndef DUNE_CURVILINEARPERIODICCONSTRUCTOR_HH
#define DUNE_CURVILINEARPERIODICCONSTRUCTOR_HH

#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <numeric>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/utility/allcommunication.hh>
#include <dune/curvilineargrid/utility/globalcommmap.hh>



namespace Dune {

namespace CurvGrid {


// Forward declaration
//template<class ct, int cdim, bool isCached>
//class CurvilinearGridStorage;


template <class GridBase>
class CurvilinearPeriodicConstructor {
public:

    /* public types */
	static const int dimension = GridBase::dimension;
	typedef typename GridBase::ctype                  ctype;
	typedef typename GridBase::GridStorageType                  GridStorageType;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
    typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    // Containers
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

    // Coordinates
	typedef typename GridStorageType::GlobalCoordinate  GlobalCoordinate;
	typedef typename GridStorageType::LocalCoordinate3D  LocalCoordinate3D;
	typedef typename GridStorageType::LocalCoordinate2D  LocalCoordinate2D;
	typedef typename GridStorageType::LocalCoordinate1D  LocalCoordinate1D;

	////////////////////////////////////////////
	// Auxiliary gather types
	////////////////////////////////////////////

	// [TODO] This method is hard-coded for triangular faces. When implementing other elements (e.g. tetrahedra), generalize
	// [TODO] Communication of CoM is excessive. Compute it after the communication
	// Note: avoid using vectors and other dynamic arrays. This type must be POD - have compile-time const size
	struct PeriodicGatherData {
		int ownerRank_;						// Rank of the process that stores this DB
		GlobalIndexType gind_;			// Global Index of this DB
		GlobalCoordinate com_;			// Center of Mass of this DB
		GlobalCoordinate corner_[3];	// Corner coordinates of this  of this DB

		PeriodicGatherData() {}
		PeriodicGatherData(int ownerRank, GlobalIndexType gind, GlobalCoordinate corner[3], GlobalCoordinate com) :
			ownerRank_(ownerRank),  gind_(gind), com_(com)
		{
			corner_[0] = corner[0];
			corner_[1] = corner[1];
			corner_[2] = corner[2];
		}

		bool operator<(const PeriodicGatherData & other) {
			GlobalCoordinate diff = com_ - other.com_;

			if			(fabs(diff[0]) > NUMERICS_RELATIVE_TOLERANCE)	{ return diff[0] < 0; }
			else if	(fabs(diff[1]) > NUMERICS_RELATIVE_TOLERANCE)	{ return diff[1] < 0; }
			else 																						{ return diff[2] < 0; }
		}
	};

	////////////////////////////////////////////
	// Auxiliary scatter types
	////////////////////////////////////////////
	struct PeriodicScatterData {
		GlobalIndexType gind_;
		int ownerRank_;
		unsigned int thisPermutationIndex_;

		PeriodicScatterData() {}
		PeriodicScatterData(GlobalIndexType gind, int ownerRank, unsigned int thisPermutationIndex) :
			gind_(gind),
			ownerRank_(ownerRank),
			thisPermutationIndex_(thisPermutationIndex)
		{}
	};

	typedef std::pair<GlobalIndexType, PeriodicScatterData> GlobalIndexPeriodicScatterDataPair;
	typedef std::map<GlobalIndexType, PeriodicScatterData> GlobalIndexPeriodicScatterDataMap;

	typedef GlobalCommMap<GlobalIndexType, PeriodicScatterData> GlobalIndexPeriodicScatterDataMapGenerator;


public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
	CurvilinearPeriodicConstructor(GridStorageType & gridstorage, GridBase & gridbase, MPIHelper &mpihelper) :
        gridstorage_(gridstorage),
		gridbase_(gridbase),
        mpihelper_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

    	// Generate Unit vectors
        // [TODO] Move to some reasonable helper class
        eCartesian_[0] = GlobalCoordinate(0.0);  eCartesian_[0][0] = 1.0;
        eCartesian_[1] = GlobalCoordinate(0.0);  eCartesian_[1][1] = 1.0;
        eCartesian_[2] = GlobalCoordinate(0.0);  eCartesian_[2][2] = 1.0;

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "Initialized CurvilinearPeriodicConstructor");
    }


    /** Communicates the Periodic Boundaries and Periodic Ghost elements
     *
     *  Prerequisites:
     *    * Exist all entities and global index
     *    * Domain boundary is a cuboid with conforming faces on the opposite directions
     *
     *  Algorithm:
     *    1) For each own DB
     *      * Find its normal, check if it is L,R,U,D,F,B.
     *      * If this dim is periodic, push to array GInd and CoM of face
     *    2) Communicate all L,R,U,D,F,B to master process
     *    3) On master process, find all pairs l-R, U-D, F-B
     *    4) Construct GCommMap for each direction
     *    5) On each process find the neighboring DB, store reduced map
     *    6) For all neighbors that are not on the same process, communicate entire elements and PB subentity index, and mark them as periodic ghosts
     *    7) For each periodic DB store local index of the neighboring local (if both on same process) or periodic Ghost entity
     *
     *    Note: It is allowed to ask a Ghost element for subentities, but not ok to ask subentities of ghost for neighbor elements
     *
     *    [TODO] Periodic gather map has excessive storage
     *    [TODO] Periodic scatter map
     *
     * */
    void generate()
    {
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Determining Periodic Boundaries");
    	std::vector<PeriodicGatherData> periodicData[6];
    	determinePeriodicDBNormalDirections(periodicData);
    	int nPeriodicFaceTotThisProcess = 0;
    	for (int i = 0; i < 6; i++) { nPeriodicFaceTotThisProcess += periodicData[i].size(); }

    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Gather PRB on Master");
    	std::vector<int> nPeriodicFaceProc[6];
    	std::vector<PeriodicGatherData> periodicDataTot[6];
    	communicatePeriodicBoundaryInfo(periodicData, periodicDataTot, nPeriodicFaceProc);

    	if (rank_ == MPI_MASTER_RANK) {
			std::stringstream logperiodic;
			logperiodic << "...Found periodic boundaries with normals LeftRightFrontBackDownUp={";
			for (int iNormDim = 0; iNormDim < 6; iNormDim++) { logperiodic << periodicDataTot[iNormDim].size() << " "; }
			LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>( __FILE__, __LINE__, logperiodic.str() + "}");
    	}

    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Determining neighbors");
    	std::vector<GlobalIndexPeriodicScatterDataPair> periodicConnectorTot;
    	determinePeriodicNeighbors(periodicDataTot, periodicConnectorTot);


    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Scatter pairs from Master");
    	GlobalIndexPeriodicScatterDataMapGenerator periodicMapGlobal;
    	periodicMapGlobal.init(mpihelper_, periodicConnectorTot);


    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assemble PRB gind map");
    	assemblePeriodicGlobalIndexLocalMap(periodicMapGlobal, nPeriodicFaceTotThisProcess);


    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Mark all local periodic neighbors into the grid");
    	markLocalPeriodicNeighbors();
    }

    const GlobalIndexPeriodicScatterDataMap & map() const { return gind2periodicDataMap_; }
    const GlobalIndexPeriodicScatterDataMap & mapInverse() const { return gind2periodicDataMapInverse_; }

    void clear()  {
    	gind2periodicDataMap_.clear();
    	gind2periodicDataMapInverse_.clear();
    }



protected:

    // Find out if two vectors are parallel, antiparallel or orthogonal, and throw error if neither is the case
    int orthogonalAssert(auto & v1, auto & v2, auto & faceGeom) {
    	ctype dot = v1 * v2;

    			if (fabs(dot) < NUMERICS_RELATIVE_TOLERANCE) { return 0; }
    	else if (fabs(dot - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { return 1; }
    	else if (fabs(dot + 1.0) < NUMERICS_RELATIVE_TOLERANCE) { return -1; }
    	else {
    		std::stringstream lerr;
    		lerr << "Dot product between unit vector " << v1 << " and the surface normal " << v2 << " is " << dot << std::endl;
    		lerr << "Errors for 1, -1 were: " << fabs(dot - 1.0) << " and " << fabs(dot + 1.0) << std::endl;
    		lerr << "Vertices of the face geometry were: " << VectorHelper::vector2string(faceGeom.vertexSet());

    		DUNE_THROW(Dune::IOError, lerr.str());
    	}
    }


    // Matching points must be translations along the periodicity axis
    bool matchPeriodicPoints(
    		const GlobalCoordinate & p1,
			const GlobalCoordinate & p2,
			int dim, bool assert)
    {
		GlobalCoordinate diff = p1 - p2;
		ctype diffMag = diff.two_norm();

		ctype errMag = fabs(diffMag - gridstorage_.periodicCuboidLength_[dim]);
		ctype errDir = ((eCartesian_[dim] * diff) / diffMag - 1.0);

		bool match = (errMag < NUMERICS_RELATIVE_TOLERANCE) && (errDir < NUMERICS_RELATIVE_TOLERANCE);

		// If match is expected, report and throw error
		if (assert && !match) {
			std::stringstream serr;
			serr << "For dim " << dim << " supposedly conforming periodic vertices=" << p1 << " and " << p2;
			serr << " have errMag=" << errMag << " and errDir=" << errDir;
			std::cerr << serr.str() << std::endl;
			DUNE_THROW(Dune::IOError, "Suspect non-conformal periodic mesh");
		}

		return match;
    }


    // Assert that two faces are translation-conformal w.r.t provided normal
    // Find the permutation of v2 corners such that it conforms to v1
    // Return the permutation index
    std::pair<unsigned int, unsigned int> matchPeriodicFace(
    		const GlobalCoordinate v1[],
			const GlobalCoordinate v2[],
			int dim)
    {
    	const unsigned int nFacePoint = 3;  // [TODO] Make me nice
    	const unsigned int NON_INIT = 10000;  // Set some unrealistic size to check if was init

    	std::vector<unsigned int> permForwards(nFacePoint, NON_INIT);			// how to permute v1 to become v2
    	std::vector<unsigned int> permBackwards(nFacePoint, NON_INIT);		// how to permute v2 to become v1

    	for (unsigned int i = 0; i < nFacePoint; i++) {
    		unsigned int nMatchThisCorner = 0;
        	for (unsigned int j = 0; j < nFacePoint; j++) {
        		if (matchPeriodicPoints(v1[i], v2[j], dim, false)) {
        			nMatchThisCorner++;
        			assert(permForwards[i] == NON_INIT);  // check surjectivity
        			assert(permBackwards[j] == NON_INIT);  // check surjectivity
        			permForwards[i] = j;
        			permBackwards[j] = i;
        		}
        	}
        	assert(nMatchThisCorner == 1);  // check injectivity
    	}

    	return std::pair<unsigned int, unsigned int> (
    			Dune::CurvilinearGeometryHelper::permutation2index(permForwards),
				Dune::CurvilinearGeometryHelper::permutation2index(permBackwards)
    	);
    }


    // Determine which dimension each of the local DB corresponds to
    // Optimization: Can skip putting into global map the periodic faces which are both local
    void determinePeriodicDBNormalDirections(std::vector<PeriodicGatherData> periodicData[6])
    {
        for (auto && faceIter : gridstorage_.boundarySegmentIndexMap_)
        {
            LocalIndexType thisFaceLocalIndex = faceIter.first;
            LocalIndexType thisFaceLocalDBIndex = faceIter.second;

            FaceStorage & thisFace_ = gridstorage_.face_[thisFaceLocalIndex];
            InternalIndexType thisFaceSubIndex = thisFace_.element1SubentityIndex;

            auto parentGeom = gridbase_.template entityGeometry<ELEMENT_CODIM>(thisFace_.element1Index);
            auto faceGeom = parentGeom.template subentityGeometry<dimension-1>(thisFaceSubIndex);
            LocalCoordinate2D faceCenterLocal2D = faceGeom.refElement().position( 0, 0 );
            LocalCoordinate3D faceCenterLocal3D = CurvilinearGeometryHelper::coordinateInParent<ctype, dimension-1, dimension>(parentGeom.type(), thisFaceSubIndex, faceCenterLocal2D);
            GlobalCoordinate faceUnitOuterNormal = parentGeom.subentityUnitNormal(thisFace_.element1SubentityIndex, faceCenterLocal3D);
            GlobalCoordinate faceCenterGlobal = faceGeom.global(faceCenterLocal2D);

            int ortX = orthogonalAssert(eCartesian_[0], faceUnitOuterNormal, faceGeom);
            int ortY = orthogonalAssert(eCartesian_[1], faceUnitOuterNormal, faceGeom);
            int ortZ = orthogonalAssert(eCartesian_[2], faceUnitOuterNormal, faceGeom);

            int ortDim;
            if			((ortX == -1)&&(ortY == 0)&&(ortZ == 0)) { ortDim = 0; }
            else if	((ortX == 1)&&(ortY == 0)&&(ortZ == 0)) { ortDim = 1; }
            else if	((ortX == 0)&&(ortY == -1)&&(ortZ == 0)) { ortDim = 2; }
            else if	((ortX == 0)&&(ortY == 1)&&(ortZ == 0)) { ortDim = 3; }
            else if	((ortX == 0)&&(ortY == 0)&&(ortZ == -1)) { ortDim = 4; }
            else if	((ortX == 0)&&(ortY == 0)&&(ortZ == 1)) { ortDim = 5; }
            else { DUNE_THROW(Dune::IOError, "Unexpected orthogonality pattern " + std::to_string(ortX) + " " + std::to_string(ortY) + " " + std::to_string(ortZ)); }

            // If it is decided that this face is a periodic face
            if (gridstorage_.periodicCuboidDimensions_[ortDim / 2]) {
				// Check that the face is originally a domain boundary, and set its structural type to periodic
				assert(gridstorage_.face_[thisFaceLocalIndex].boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary);
				gridstorage_.face_[thisFaceLocalIndex].boundaryType = GridStorageType::FaceBoundaryType::PeriodicBoundary;

				// Add its data to the gather data array
            	periodicData[ortDim].push_back(PeriodicGatherData(rank_, thisFace_.globalIndex, faceGeom.cornerSet().data(), faceCenterGlobal));
            }
        }
    }


    // Gather on Master the GInd and CoM of all periodic faces, subdivided by by normal arrays
    void communicatePeriodicBoundaryInfo(
        	std::vector<PeriodicGatherData> periodicData[6],
        	std::vector<PeriodicGatherData> periodicDataTot[6],
        	std::vector<int> nPeriodicFaceProc[6])
    {
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	Dune::CurvGrid::AllCommunication allcomm(mpihelper_);

    	for (int iNormDim = 0; iNormDim < 6; iNormDim++) {
    		if (gridstorage_.periodicCuboidDimensions_[iNormDim / 2]) {

        		allcomm.gatherv(periodicData[iNormDim], periodicDataTot[iNormDim], nPeriodicFaceProc[iNormDim]);

        		// Sanity self-check
            	if (rank_ == MPI_MASTER_RANK)	{ assert(nPeriodicFaceProc[iNormDim].size() == size_); }
            	else													{ assert(nPeriodicFaceProc[iNormDim].size() == 0); }
    		}
    	}
    }


    void determinePeriodicNeighbors(
    		std::vector<PeriodicGatherData> periodicDataTot[6],
			std::vector<GlobalIndexPeriodicScatterDataPair> & periodicConnectorTot)
    {
    	if (rank_ == MPI_MASTER_RANK) {
			ctype posPlane[6];
			for (int iDim = 0; iDim < dimension; iDim++) {
				if (gridstorage_.periodicCuboidDimensions_[iDim]) {
					for (int iDir = 0; iDir < 2; iDir++) {
						// Find plane coordinates for each direction, verify that all faces are coplanar
						int iNormDim = iDim * 2 + iDir;
						posPlane[iNormDim] = periodicDataTot[iNormDim][0].com_[iDim];
						for (auto && gc : periodicDataTot[iNormDim]) { assert(fabs(posPlane[iNormDim] - gc.com_[iDim]) < NUMERICS_RELATIVE_TOLERANCE); }

						// Sort all faces in X > Y > Z order of CoM
						std::sort(periodicDataTot[iNormDim].begin(), periodicDataTot[iNormDim].end());
					}

					// Define directions
					int iNormDimP = iDim * 2 + 1;
					int iNormDimM = iDim * 2;
					assert(periodicDataTot[iNormDimM].size() == periodicDataTot[iNormDimP].size());

					// Determine the distance between periodic planes
					gridstorage_.periodicCuboidLength_[iDim] = posPlane[iNormDimP] - posPlane[iNormDimM];

					// Verify that the sorted faces exactly correspond to each other by comparing the CoM distance
					for (int iFace = 0; iFace < periodicDataTot[iNormDimM].size(); iFace++) {

						// Note: Since the routine performs the assert itself, the return value is irrelevant
						matchPeriodicPoints(
								periodicDataTot[iNormDimP][iFace].com_,
								periodicDataTot[iNormDimM][iFace].com_,
								iDim, true);

						// Match the corresponding corners, thus get the permutation index of (+) face, when keeping the (-) face orientation fixed
						// VERY IMPORTANT: THIS STRATEGY RELIES ON LOCAL COORDINATE SYSTEM OF EACH ELEMENT STAYING FIXED
						// IF GRID DECIDES TO PERMUTE THE LOCAL COORDINATES OF AN ELEMENT, THIS METHOD WILL BECOME INVALID.
						std::pair<unsigned int, unsigned int> faceCoordPermutation = matchPeriodicFace(
								periodicDataTot[iNormDimM][iFace].corner_,
								periodicDataTot[iNormDimP][iFace].corner_,
								iDim
						);
						unsigned int permIndexM = faceCoordPermutation.first;
						unsigned int permIndexP = faceCoordPermutation.second;

						// Store the map in both directions
						// NOTE:The permutation index corresponds to permuting arg face, not val face
						// [TODO] Communicating 2-directional map is unnecessary. Can construct a 1-directional map, communicate, then construct the opposite direction
						GlobalIndexType gindm = periodicDataTot[iNormDimM][iFace].gind_;
						GlobalIndexType gindp = periodicDataTot[iNormDimP][iFace].gind_;
						periodicConnectorTot.push_back(GlobalIndexPeriodicScatterDataPair(gindm, PeriodicScatterData(gindp, periodicDataTot[iNormDimP][iFace].ownerRank_, permIndexM)));
						periodicConnectorTot.push_back(GlobalIndexPeriodicScatterDataPair(gindp, PeriodicScatterData(gindm, periodicDataTot[iNormDimM][iFace].ownerRank_, permIndexP)));
					}
				}
			}
		}
    }


    void assemblePeriodicGlobalIndexLocalMap(const GlobalIndexPeriodicScatterDataMapGenerator & periodicMapGlobal, int nPeriodicFaceTotThisProcess)
    {
    	const GlobalIndexPeriodicScatterDataMap & gind2gindrankmap = periodicMapGlobal.map();

    	int countFoundPeriodicFaces = 0;
    	for (auto && faceIter : gridstorage_.boundarySegmentIndexMap_) {
			LocalIndexType thisFaceLocalIndex = faceIter.first;
			const FaceStorage & thisFace = gridstorage_.face_[thisFaceLocalIndex];
			LocalIndexType thisFaceGlobalIndex = thisFace.globalIndex;

			auto g2grmInnerIter = gind2gindrankmap.find(thisFaceGlobalIndex);
			if (g2grmInnerIter != gind2gindrankmap.end()) {
				// Check that this boundary is already marked as a periodic boundary
				assert(gridstorage_.face_[thisFaceLocalIndex].boundaryType == GridStorageType::FaceBoundaryType::PeriodicBoundary);

				// Count the periodic face
				countFoundPeriodicFaces++;

				// Filter out the periodic boundary connection map only with faces local to this process
				PeriodicScatterData neighborData = g2grmInnerIter->second;
				gind2periodicDataMap_[thisFaceGlobalIndex] = neighborData;
				gind2periodicDataMapInverse_[neighborData.gind_] = PeriodicScatterData(thisFaceGlobalIndex, rank_, 0); // [FIXME] Permutation index set to whatever number since it is not used

				// Generate the periodic boundary local index
				assert(gridstorage_.periodicBoundaryIndexMap_.find(thisFaceLocalIndex) == gridstorage_.periodicBoundaryIndexMap_.end());
				LocalIndexType periodicLocalIndex = gridstorage_.periodicBoundaryIndexMap_.size();
				gridstorage_.periodicBoundaryIndexMap_[thisFaceLocalIndex] = periodicLocalIndex;
				gridstorage_.PERB2PERBNeighborRank_[FACE_CODIM].push_back(std::vector<int>{neighborData.ownerRank_});

				// Store the periodic permutations
				auto g2grmOuterIter = gind2gindrankmap.find(g2grmInnerIter->second.gind_);
				assert(g2grmOuterIter != gind2gindrankmap.end());
				assert(g2grmOuterIter->second.gind_ == thisFaceGlobalIndex);

				gridstorage_.periodicFaceMatchPermutationIndexInner_.push_back(g2grmInnerIter->second.thisPermutationIndex_);
				gridstorage_.periodicFaceMatchPermutationIndexOuter_.push_back(g2grmOuterIter->second.thisPermutationIndex_);
			}
		}

    	assert(gind2periodicDataMap_.size() == gind2periodicDataMapInverse_.size());

    	// Check that all periodic faces were present in the map
    	assert(countFoundPeriodicFaces == nPeriodicFaceTotThisProcess);
    }


    // [TODO] Adapt for the case with no ghost elements whatsoever - then the 2nd neighbor of face should be another periodic face
    void markLocalPeriodicNeighbors() {
        for (auto const & faceIter: gind2periodicDataMap_)
		{
			// Get global and local face index
			GlobalIndexType thisFaceGlobalIndex = faceIter.first;
			GlobalIndexType thisNeighborFaceGInd = faceIter.second.gind_;
			auto thisFaceIndexIter = gridstorage_.entityIndexMap_[FACE_CODIM].find(thisFaceGlobalIndex);
			assert(thisFaceIndexIter != gridstorage_.entityIndexMap_[FACE_CODIM].end());
			LocalIndexType thisFaceLocalIndex = thisFaceIndexIter->second;

			// If periodic neighbor exists on this process, do not communicate it, but mark it straight away
			// Optimization: The assignment can be done in pairs
			auto thisNeighborFaceIndexIter = gridstorage_.entityIndexMap_[FACE_CODIM].find(thisNeighborFaceGInd);
			if (thisNeighborFaceIndexIter != gridstorage_.entityIndexMap_[FACE_CODIM].end()) {
				LocalIndexType thisNeighborFaceLInd = thisNeighborFaceIndexIter->second;
				gridstorage_.face_[thisFaceLocalIndex].element2Index = gridstorage_.face_[thisNeighborFaceLInd].element1Index;
				gridstorage_.face_[thisFaceLocalIndex].element2SubentityIndex = gridstorage_.face_[thisNeighborFaceLInd].element1SubentityIndex;

				assert(gridstorage_.face_[thisNeighborFaceLInd].element1Index >= 0);
				assert(gridstorage_.face_[thisNeighborFaceLInd].element1SubentityIndex >= 0);
			}
		}
    }




private: // Private members

    // Cartesian Unit Vectors
    // [TODO] Move to some reasonable helper class
    GlobalCoordinate eCartesian_[dimension];

    GridStorageType & gridstorage_;
    GridBase & gridbase_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;

    // Map
    GlobalIndexPeriodicScatterDataMap gind2periodicDataMap_;				// From this process Periodic face Global index to Neighbor process periodic face data
    GlobalIndexPeriodicScatterDataMap gind2periodicDataMapInverse_;	// From neighbor process Periodic face Global index to this process periodic face data
};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARPERIODICCONSTRUCTOR_HH
