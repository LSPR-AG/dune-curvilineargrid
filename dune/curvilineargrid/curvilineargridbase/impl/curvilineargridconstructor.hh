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
    description          : Upgraded the mesh to curvilinear grid,
      implemented global index, ghost elements, dense and sparse communication, interior and periodic boundaries
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEARGRIDCONSTRUCTOR_HH
#define DUNE_CURVILINEARGRIDCONSTRUCTOR_HH

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
#include <dune/common/parallel/mpicommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/loggingtimer.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/impl/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearglobalindexconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearghostconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearperiodicconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearpostconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearlooseoctree.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>


namespace Dune {

namespace CurvGrid {


// Forwards-declatation of the base class
// **********************************************
//template <class ct, int cdim, bool isCached>
//class CurvilinearGridBase;


// Constructor class
// **********************************************
template <class GridBase>
class CurvilinearGridConstructor {
public:

    /* public types */
	typedef typename GridBase::ctype   ctype;
	static const int dimension = GridBase::dimension;
	typedef typename Dune::ReferenceElements<ctype, dimension>  ReferenceElements;

	typedef          GridBase                                   GridBaseType;
	typedef typename GridBase::GridStorageType                  GridStorageType;
    typedef typename GridBase::LoggingTimer                     LoggingTimer;

    typedef CurvilinearGhostConstructor<GridBase>         GridGhostConstructor;
    typedef CurvilinearPostConstructor<GridBase>		    GridPostConstructor;
    typedef CurvilinearGlobalIndexConstructor<GridBase>   GridGlobalIndexConstructor;
    typedef CurvilinearPeriodicConstructor<GridBase>  GridPeriodicConstructor;

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


	typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;
    static const unsigned int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridConstructor(GridBaseType & gridbase, MPIHelper &mpihelper) :
        gridstorage_(gridbase.gridstorage()),
        gridbase_(gridbase),
        mpihelper_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "Initialized CurvilinearGridConstructor");
    }


public:

    /* ***************************************************************************
     * Section: Loading the mesh
     * ***************************************************************************/


    /** \brief Add a new vertex to the mesh
     * \param[in] p                coordinate of this vertex
     * \param[in] globalIndex      global index of this vertex
     * \param[in] structtype       (optional) partition type of this vertex. User SHOULD use default value
     * */
    void insertVertex(GlobalCoordinate p, GlobalIndexType globalIndex, Dune::PartitionType ptype = Dune::PartitionType::InteriorEntity)
    {
        VertexStorage point;
        point.coord = p;
        point.globalIndex = globalIndex;
        point.ptype = ptype;

        gridstorage_.entityIndexMap_[VERTEX_CODIM][globalIndex] = gridstorage_.point_.size();
        gridstorage_.point_.push_back(point);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted vertex LocalIndex=" << gridstorage_.point_.size()-1 << " GlobalIndex=" << globalIndex;
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
    }


    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] elementIndex     global index of the element provided by the GMSH reader
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * Note: Even though we pass the globalId as a parameter from GMSH, it is a globalIndex for the set Elements+BoundaryFaces,
     * therefore obtaining globalIndex for elements from it is not possible, it will have to be communicated
     *
     * */
    void insertElement(
    	Dune::GeometryType gt,
    	const std::vector<LocalIndexType> & vertexIndexSet,
    	GlobalIndexType elementIndex,
    	InterpolatoryOrderType order,
    	PhysicalTagType physicalTag)
    {
        if (!gt.isTetrahedron() || (vertexIndexSet.size() != CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
        	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertElement() unexpected element type or number of interpolatory points");
            DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertElement() unexpected element type or number of interpolatory points");
        }

        EntityStorage thisElement;

        thisElement.geometryType = gt;
        thisElement.globalIndex = elementIndex;        // At this stage the globalIndex might not be defined if GMSH does not provide it
        thisElement.ptype = Dune::PartitionType::InteriorEntity;
        thisElement.interpOrder = order;
        thisElement.physicalTag = physicalTag;
        thisElement.vertexIndexSet = vertexIndexSet;

        LocalIndexType thisLocalIndex = gridstorage_.element_.size();
        gridstorage_.element_.push_back(thisElement);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted Element Type=" << CurvilinearGeometryHelper::geometryName(gt);
        log_stream << " LocalIndex=" << thisLocalIndex;
        log_stream << " Order=" << order;
        log_stream << " PhysicalTag=" << physicalTag;
        log_stream << " VertexIndices=(" << VectorHelper::vector2string(vertexIndexSet) << ")";
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
    }


    /** Insert a boundary segment into the mesh
     *
     *     Note: It is expected that all domain boundaries are inserted by the factory before finalising
     *     Note: Only domain boundary faces have initial globalId given by GMSH. Therefore, we ignore it, and generate our own
     *     globalId for all faces at a later stage.
     *
     *     Note: Support for interior boundaries is optional. There is no restriction on the number of inserted interior boundaries
     *     The final effect of this operation is that the faces corresponding to the inserted boundary segments will be marked
     *     with the provided physicalTag, to distinguish them from the other faces
     *
     *  \param[in] gt                       geometry type of the face (should be a triangle)
     *  \param[in] vertexIndexSet           local indices of the interpolatory vertices of this face
     *  \param[in] order                    interpolatory order of the face
     *  \param[in] physicalTag              physical tag of the element (material property)
     *  \param[in] isDomainBoundary              determines whether the inserted boundary segment belongs to the domain or interior boundary
     *
     * */
    void insertBoundarySegment(
    	Dune::GeometryType gt,
    	//LocalIndexType associatedElementIndex,
    	const std::vector<LocalIndexType> & vertexIndexSet,
    	InterpolatoryOrderType order,
    	PhysicalTagType physicalTag,
		bool isDomainBoundary)
    {
        if (!gt.isTriangle() || (vertexIndexSet.size() != CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
        	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
            DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
        }


        // Get corners of this face
        // **********************************************************************************
        std::vector<LocalIndexType> faceCorners = CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 2>(gt, vertexIndexSet, order);
        FaceKey thisFaceKey;
        thisFaceKey.node0 = gridstorage_.point_[faceCorners[0]].globalIndex;
        thisFaceKey.node1 = gridstorage_.point_[faceCorners[1]].globalIndex;
        thisFaceKey.node2 = gridstorage_.point_[faceCorners[2]].globalIndex;

        // Sort in ascending order
        thisFaceKey.sort();

        // Map the face key to the associated physical tag
        if (isDomainBoundary)	{ domainBoundaryFaceKey2TagMap_[thisFaceKey] = physicalTag; }
        else									{ interiorBoundaryFaceKey2TagMap_[thisFaceKey] = physicalTag; }
    }


    /** \brief Compulsory: insert the total number of vertices in the mesh before constructing the grid */
    void insertNVertexTotal(int nVertexTotal)  { gridstorage_.nEntityTotal_[VERTEX_CODIM] = nVertexTotal; }

    /** \brief Compulsory: insert the total number of elements in the mesh before constructing the grid */
    void insertNElementTotal(int nElementTotal)  { gridstorage_.nEntityTotal_[ELEMENT_CODIM] = nElementTotal; }


    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    /** Calls the subroutines for transforming the inserted data into a functional mesh.
     *  Note: It is expected, that all necessary data (vertices, elements and boundary segments) have been added before this function is called.
     * */

    // 0) Compute process bounding box - coordinates of a box - compact, but sufficient to fit in the elements of this process
    // 1) Generate all edges
    // 2) Generate all faces - process boundary and internal (as domain boundary already inserted)
    // 3) Generate process boundary corner set - from process boundary faces
    // 4) Generate process boundary edge set  - from process boundary faces
    // 5) Generate Global Indices
    //        - Communicate neighbor ranks for all process boundary corners
    //      - Compute neighbor ranks for all edges from set intersection. Mark as owned if your rank is top
    //      - Compute neighbor ranks for all faces from set intersection. Mark as owned if your rank is top
    //      - [Nothing changes here] Enumerate owned
    //      - Communicate non-owned. No such thing as complex vertex at this stage - just communicate to all other processes who own it
    // 6) Generate Ghost Elements
    // 7) Construct OCTree

    void generateMesh() {
    	LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "[[CurvilinearGridBase: Generating curvilinear mesh...");
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Initializing mesh");

        // Construct missing parts of the mesh
        // ************************************************************
        gridstorage_.nEntityTotal_[EDGE_CODIM] = 0;  // Will be updated later
        gridstorage_.nEntityTotal_[FACE_CODIM] = 0;  // Will be updated later

        LoggingTimer::time("CurvilinearGridConstructor: EntityGeneration");
        generateEdges();
        generateFaces();
        markBorderVertex();
        markBorderEdge();
        LoggingTimer::time("CurvilinearGridConstructor: EntityGeneration");


		// Generate Global Index
		// **********************************************************
        if (size_ > 1)  // Parallel case
		{
			LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Parallel Grid");

			LoggingTimer::time("CurvilinearGridConstructor: Global Index Generation");

			GridGlobalIndexConstructor giConstructor(
				gridstorage_,
				gridbase_,
				mpihelper_,

				edgeKey2LocalIndexMap_,
				internalFaceKey2LocalIndexMap_,
				domainBoundaryFaceKey2LocalIndexMap_,
				processBoundaryFaceKey2LocalIndexMap_);

			giConstructor.generateEdgeGlobalIndex();
			giConstructor.generateFaceGlobalIndex();
			if (!gridstorage_.withElementGlobalIndex_) { giConstructor.generateElementGlobalIndex(); }   // Generate element global index unless it was explicitly provided by the factory

			// Fill in Global2Local maps
			for (unsigned int iEdge = 0; iEdge < gridstorage_.edge_.size(); iEdge++)    { gridstorage_.entityIndexMap_[EDGE_CODIM][gridstorage_.edge_[iEdge].globalIndex] = iEdge; }
			for (unsigned int iFace = 0; iFace < gridstorage_.face_.size(); iFace++)    { gridstorage_.entityIndexMap_[FACE_CODIM][gridstorage_.face_[iFace].globalIndex] = iFace; }
			for (unsigned int iElem = 0; iElem < gridstorage_.element_.size(); iElem++) {	gridstorage_.entityIndexMap_[ELEMENT_CODIM][gridstorage_.element_[iElem].globalIndex] = iElem; }

			LoggingTimer::time("CurvilinearGridConstructor: Global Index Generation");
		}
		else  // Serial case
		{
			LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Serial Grid");

			// Serial case:
			// * Boundary Neighbors not necessary, since all boundaries are domain boundaries
			// * No ghost elements, even if requested by user
			// * Fake globalIndex by making it equal to localIndex

			LoggingTimer::time("CurvilinearGridConstructor: Global Index Generation");

			gridstorage_.withGhostElements_ = false;
			for (unsigned int i = 0; i < gridstorage_.edge_.size();    i++)  { gridstorage_.edge_[i].globalIndex = i;     gridstorage_.entityIndexMap_[EDGE_CODIM][i] = i; }
			for (unsigned int i = 0; i < gridstorage_.face_.size();    i++)  { gridstorage_.face_[i].globalIndex = i;     gridstorage_.entityIndexMap_[FACE_CODIM][i] = i; }
			for (unsigned int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = i;  gridstorage_.entityIndexMap_[ELEMENT_CODIM][i] = i;  gridstorage_.entityInternalIndexSet_[ELEMENT_CODIM].insert(i); }

			for (FaceMapIterator faceIter = internalFaceKey2LocalIndexMap_.begin(); faceIter != internalFaceKey2LocalIndexMap_.end(); faceIter++)              { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.entityInternalIndexSet_[FACE_CODIM].insert(localIndex); }
			for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)  { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.faceDomainBoundaryIndexSet_.insert(localIndex); }

			gridstorage_.nEntityTotal_[EDGE_CODIM] = gridstorage_.edge_.size();
			gridstorage_.nEntityTotal_[FACE_CODIM] = gridstorage_.face_.size();

			LoggingTimer::time("CurvilinearGridConstructor: Global Index Generation");
		}


        // Deletes all temporary memory
        // ************************************************************
        edgeKey2LocalIndexMap_.clear();
        internalFaceKey2LocalIndexMap_.clear();
        domainBoundaryFaceKey2LocalIndexMap_.clear();
        processBoundaryFaceKey2LocalIndexMap_.clear();


        // Construct periodic boundary
        // ************************************************************
        GridPeriodicConstructor periodicConstructor(gridstorage_, gridbase_, mpihelper_);
        if (gridstorage_.periodicCuboidDimensions_.size() > 0) { periodicConstructor.generate(); }


    	// Generate Ghost elements
    	// **********************************************************
        if (size_ > 1)  // Parallel case
        {
            if (gridstorage_.withGhostElements_)
            {
            	LoggingTimer::time("CurvilinearGridConstructor: Ghost Element Constructor");

            	if (gridstorage_.periodicCuboidDimensions_.size() > 0)		{
            		GridGhostConstructor ghostConstructor(gridstorage_, mpihelper_, &periodicConstructor);
            		ghostConstructor.generate();
            	} else {
            		GridGhostConstructor ghostConstructor(gridstorage_, mpihelper_);
            		ghostConstructor.generate();
            	}

            	LoggingTimer::time("CurvilinearGridConstructor: Ghost Element Constructor");
            }
        }


        // Free up memory by deleting the periodic constructor
        // ************************************************************
        periodicConstructor.clear();


        // Create sets that will be used for iteration over the map
        // Create maps of all entity subsets which can be communicated over
        // Find neighbor ranks for all entities that can be communicated over
        // ************************************************************
        LoggingTimer::time("CurvilinearGridConstructor: Index Set Generation");
        GridPostConstructor postConstructor(gridstorage_, gridbase_, mpihelper_);
        postConstructor.generateCornerIndex();
        postConstructor.generateIteratorSets();
        LoggingTimer::time("CurvilinearGridConstructor: Index Set Generation");



        // The PB-PB communication interface is available by default. The below procedures enable the communication interfaces
        // involving ghost entities, and require existence of ghost entities
        if ((size_ > 1)&& gridstorage_.withGhostElements_)
        {
#if HAVE_MPI
        	LoggingTimer::time("CurvilinearGridConstructor: Generation of communication maps");
        	postConstructor.generateCommunicationMaps();
        	LoggingTimer::time("CurvilinearGridConstructor: Generation of communication maps");

        	// Does not exist any more
//        	LoggingTimer::time("CurvilinearGridConstructor: Communication of neighbor ranks");
//            postConstructor.communicateCommunicationEntityNeighborRanks();
//            LoggingTimer::time("CurvilinearGridConstructor: Communication of neighbor ranks");
#endif
        }


        // ************************************************************
        // Construct OCTree
        // ************************************************************

        //LoggingTimer::time("CurvilinearGridConstructor: Constructing OCTree");
        computeProcessBoundingBox();
        //constructOctree();
        //LoggingTimer::time("CurvilinearGridConstructor: Constructing OCTree");


        // ************************************************************
        // Diagnostics output
        // ************************************************************
        std::stringstream log_stream;
        log_stream << "CurvilinearGridBase: Constructed Mesh ";
        log_stream << " nVertexPerMesh="             << gridbase_.property().nEntityTotal(VERTEX_CODIM);
        log_stream << " nEdgePerMesh="               << gridbase_.property().nEntityTotal(EDGE_CODIM);
        log_stream << " nFacePerMesh="               << gridbase_.property().nEntityTotal(FACE_CODIM);
        log_stream << " nElementPerMesh="            << gridbase_.property().nEntityTotal(ELEMENT_CODIM);

        log_stream << std::endl; "    *** ";
        log_stream << " nCorner="                    << gridbase_.property().nEntity(VERTEX_CODIM);
        log_stream << " nCornerInterior="            << gridbase_.property().nEntity(VERTEX_CODIM, PartitionType::InteriorEntity);
        log_stream << " nCornerBorder="              << gridbase_.property().nEntity(VERTEX_CODIM, PartitionType::BorderEntity);
        log_stream << " nCornerGhost="               << gridbase_.property().nEntity(VERTEX_CODIM, PartitionType::GhostEntity);

        log_stream << std::endl; "    *** ";
        log_stream << " nEdge="                      << gridbase_.property().nEntity(EDGE_CODIM);
        log_stream << " nEdgeInterior="              << gridbase_.property().nEntity(EDGE_CODIM, PartitionType::InteriorEntity);
        log_stream << " nEdgeBorder="                << gridbase_.property().nEntity(EDGE_CODIM, PartitionType::BorderEntity);
        log_stream << " nEdgeGhost="                 << gridbase_.property().nEntity(EDGE_CODIM, PartitionType::GhostEntity);

        log_stream << std::endl; "    *** ";
        log_stream << " nFace="                      << gridbase_.property().nEntity(FACE_CODIM);
        log_stream << " nFaceInterior="              << gridbase_.property().nEntity(FACE_CODIM, PartitionType::InteriorEntity);
        log_stream << " nFaceBoundarySegment="       << gridbase_.property().nEntity(FACE_CODIM, PartitionType::InteriorEntity, DOMAIN_BOUNDARY_TYPE);
        log_stream << " nFaceBorder="                << gridbase_.property().nEntity(FACE_CODIM, PartitionType::BorderEntity);
        log_stream << " nFaceGhost="                 << gridbase_.property().nEntity(FACE_CODIM, PartitionType::GhostEntity);

        log_stream << std::endl; "    *** ";
        log_stream << " nElement="                   << gridbase_.property().nEntity(ELEMENT_CODIM);
        log_stream << " nElementInterior="           << gridbase_.property().nEntity(ELEMENT_CODIM, PartitionType::InteriorEntity);
        log_stream << " nElementGhost="              << gridbase_.property().nEntity(ELEMENT_CODIM, PartitionType::GhostEntity);

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "...CurvilinearGridBase: Finished generating curvilinear mesh]]");
    }




protected:



    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    /** Compute center and extent (halved) of the bounding box of the whole mesh.
     *
     * FIXME: Current implementation underestimates bounding box due to curvilinear effects.
     * Suggested correction:
     * 1) Min and max over all interpolation points [oh, this is already done]
     * 2) Correction factor for Excess curvature
     * 2.1) For example, compute average element linear size, and take half of that, and enlarge extent by it
     *
     * */
    void computeProcessBoundingBox()
    {
        GlobalCoordinate min = gridstorage_.point_[0].coord;
        GlobalCoordinate max = min;

        for (unsigned int i = 1; i < gridstorage_.point_.size(); i ++) {
        	LoggingMessage::writePatience(" Computing process bounding box...", i, gridstorage_.point_.size());

            min[0] = std::min(min[0], gridstorage_.point_[i].coord[0]);
            min[1] = std::min(min[1], gridstorage_.point_[i].coord[1]);
            min[2] = std::min(min[2], gridstorage_.point_[i].coord[2]);

            max[0] = std::max(max[0], gridstorage_.point_[i].coord[0]);
            max[1] = std::max(max[1], gridstorage_.point_[i].coord[1]);
            max[2] = std::max(max[2], gridstorage_.point_[i].coord[2]);
        }
        gridstorage_.boundingBoxCenter_ = min + max;  gridstorage_.boundingBoxCenter_ *= 0.5;
        gridstorage_.boundingBoxExtent_ = max - min;  gridstorage_.boundingBoxExtent_ *= 0.5;
    }


    /** \brief Generates all edges
     *
     * Algorithm:
     * 1) Loop over all elements
     * 2) For each element, construct EdgeKeys of all its edges, based on the corners of the element
     * 3) Add each EdgeKey to the map. If the edge did not exist before, give it a local index
     * 4) Note for each edge the local index of the parent that created it
     * 5) Mark edge local index as a subentity of containing element
     *
     * \note Orientation of elements ensured via Dune::ReferenceElement::subEntity()
     * [TODO]  Mark domain boundary edges correctly
     * [TODO]  Use more generic functions when extending to any mesh other than tetrahedral
     *
     * */
    void generateEdges()
    {
    	// Init the subentity index vector
    	int nEdgePerTetrahedron = 6;  // [FIXME] Use ref elem subentity size
    	int nElem = gridstorage_.element_.size();
    	gridstorage_.elementSubentityCodim2_ = std::vector<std::vector<int> > (nElem, std::vector<int>(nEdgePerTetrahedron));

        // Loop over all elements and their edges
        for (int iElem = 0; iElem < nElem; iElem++)
        {
        	LoggingMessage::writePatience("Generating edges...", iElem, nElem);

        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<LocalIndexType> elementCornerLocalIndexSet = CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

            for (int iEdge = 0; iEdge < nEdgePerTetrahedron; iEdge++)
            {

                // Get internal indices of the corners of this face wrt its associated element
                InternalIndexType node0SubIndex = ReferenceElements::general(thisElem.geometryType).subEntity(iEdge, EDGE_CODIM, 0, VERTEX_CODIM);
                InternalIndexType node1SubIndex = ReferenceElements::general(thisElem.geometryType).subEntity(iEdge, EDGE_CODIM, 1, VERTEX_CODIM);

                // Define (key = sorted globalIndices of corners)
                EdgeKey thisKey;
                thisKey.node0 = gridstorage_.point_[elementCornerLocalIndexSet[node0SubIndex]].globalIndex;
                thisKey.node1 = gridstorage_.point_[elementCornerLocalIndexSet[node1SubIndex]].globalIndex;

                // Sort in ascending order
                thisKey.sort();
                //std::cout << "process_" << rank_ << "GenerateEdgeKey=(" << thisKey.node0 << "," << thisKey.node1 << ")" << std::endl;


                EdgeMapIterator edgeIter = edgeKey2LocalIndexMap_.find(thisKey);

                // If this edge has not been added already, add it to the map
                // Find local index of this edge and note it as this element subentity
                if (edgeIter == edgeKey2LocalIndexMap_.end())
                {
                    // Store map (key -> edgeIndex)
                	LocalIndexType localEdgeIndex = gridstorage_.edge_.size();
                    edgeKey2LocalIndexMap_[thisKey] = localEdgeIndex;

                    // Store vector (edgeId -> elemId + edgeElemIndex)
                    // Note: Edges do not have physical tag at all so we do not even store it
                    EdgeStorage thisEdge;
                    thisEdge.globalIndex = 0;        // GlobalId for edge determined later using global communication
                    thisEdge.ptype = Dune::PartitionType::InteriorEntity;   // For process boundaries will be redefined later
                    thisEdge.elementIndex = iElem;
                    thisEdge.subentityIndex = iEdge;

                    // Log output
                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: Added Edge";
                    log_stream << " LocalIndex=" << localEdgeIndex;
                    log_stream << " AssociatedElementIndex=" << iElem;
                    log_stream << " InternalSubentityIndex=" << iEdge;
                    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

                    gridstorage_.edge_.push_back(thisEdge);

                    gridstorage_.elementSubentityCodim2_[iElem][iEdge] = localEdgeIndex;
                } else {
                	gridstorage_.elementSubentityCodim2_[iElem][iEdge] = (*edgeIter).second;
                }
            }
        }
    }


    /** Generates Internal and ProcessBoundary Faces. (!!!) Assumes that all Domain Boundary Faces have been added.
     *
     * Algorithm:
     * 1) Loop over all elements
     * 2) For each element, construct FaceKeys of all its edges, based on the corners of the element
     * 3) Add each FaceKey to a temporary map, unless it is already in the DomainBoundary map. In the latter case, assign its local index as the containing element subentity
     * 4) For each FaceKey in the temporary map, check if it has been added once or twice, this determines if it is a process boundary or an internal face.
     * 5) Add each face to the total face map and to the corresponding structure map (internal/processboundary). Mark face local index as a subentity of containing element
     * 6) Note for each face the local index of the parent that created it
     * 7) For each process boundary face, create map from localFaceIndex to dummy creation index
     * 8) Resize the future neighbor rank array with the number of PB faces
     *
     * \note Orientation of elements ensured via Dune::ReferenceElement::subEntity()
     * [TODO]  If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
     *
     * */
    void generateFaces()
    {
    	// Init the subentity index vector
    	// [FIXME] Use ref elem subentity size
        int nFacePerTetrahedron = 4;
    	int nElem = gridstorage_.element_.size();
    	gridstorage_.elementSubentityCodim1_ = std::vector<std::vector<LocalIndexType> > (nElem, std::vector<LocalIndexType>(nFacePerTetrahedron));

        typedef std::map<FaceKey, std::vector<int>> tmpFace2InfoMap;
        typedef typename tmpFace2InfoMap::iterator  tmpMapIterator;
        tmpFace2InfoMap tmpFaceMap;

        // Loop over all elements and their faces
        for (int iElem = 0; iElem < nElem; iElem++)
        {
        	LoggingMessage::writePatience("Generating faces...", iElem, nElem);

        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<int> elementCornerLocalIndexSet = CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

            // Store info for all faces except of domain boundaries
            // Store it in a map, not to store internal faces twice
            for (int iFace = 0; iFace < nFacePerTetrahedron; iFace++)
            {

                // Get internal indices of the corners of this face wrt its associated element
                InternalIndexType node0SubIndex = ReferenceElements::general(thisElem.geometryType).subEntity(iFace, FACE_CODIM, 0, VERTEX_CODIM);
                InternalIndexType node1SubIndex = ReferenceElements::general(thisElem.geometryType).subEntity(iFace, FACE_CODIM, 1, VERTEX_CODIM);
                InternalIndexType node2SubIndex = ReferenceElements::general(thisElem.geometryType).subEntity(iFace, FACE_CODIM, 2, VERTEX_CODIM);

                // Define (key = sorted globalIndices of corners)
                FaceKey thisKey;
                thisKey.node0 = gridstorage_.point_[elementCornerLocalIndexSet[node0SubIndex]].globalIndex;
                thisKey.node1 = gridstorage_.point_[elementCornerLocalIndexSet[node1SubIndex]].globalIndex;
                thisKey.node2 = gridstorage_.point_[elementCornerLocalIndexSet[node2SubIndex]].globalIndex;

                // Sort in ascending order
                thisKey.sort();

                std::vector<int> connectedFaceInfo;
                tmpMapIterator iter = tmpFaceMap.find(thisKey);

                // If the face already exists, update it by adding extra related elements
                if (iter != tmpFaceMap.end()) { connectedFaceInfo = std::vector<int> ( (*iter).second ); }
                connectedFaceInfo.push_back(iElem);
                connectedFaceInfo.push_back(iFace);
                tmpFaceMap[thisKey] = connectedFaceInfo;

                std::stringstream log_stream;
                log_stream << "CurvilinearGridConstructor: Adding FaceKey=(" << thisKey.node0 << ", " << thisKey.node1 << ", " << thisKey.node2 << ") attached to total of " << connectedFaceInfo.size() / 2 << " elements";
                LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());


                // FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.find(thisKey);

                // Mark this face for creation if it is not an already existing Domain Boundary
                // Otherwise note its local index
                //if (faceIter == domainBoundaryFaceKey2LocalIndexMap_.end())
                //{

                //} else {
                	//LocalIndexType localFaceIndex = (*faceIter).second;
                	//gridstorage_.elementSubentityCodim1_[iElem][iFace] = localFaceIndex;
                //}
            }
        }


        // Add internal and process boundary faces to the mesh
        int dbCount = 0;  // Count DB and IB to later check if all have been mapped
        int ibCount = 0;
        for (tmpMapIterator iter = tmpFaceMap.begin();  iter != tmpFaceMap.end();  iter++)
        {
            FaceStorage thisFace;
            LocalIndexType localFaceIndex = gridstorage_.face_.size();
            std::vector<int> connectedFaceInfo = (*iter).second;

            // Each face should have 1 or 2 neighbor elements
            assert((connectedFaceInfo.size() == 2) || (connectedFaceInfo.size() == 4));
            LocalIndexType    thisAssociatedElementIndex = connectedFaceInfo[0];
            InternalIndexType thisFaceSubentityIndex = connectedFaceInfo[1];

            // Store the face local index as element subentity
            gridstorage_.elementSubentityCodim1_[thisAssociatedElementIndex][thisFaceSubentityIndex] = localFaceIndex;

            // Recover parental information for this face
            thisFace.globalIndex = 0;       // GlobalId is defined at a later stage
            thisFace.element1Index = thisAssociatedElementIndex;
            thisFace.element1SubentityIndex = thisFaceSubentityIndex;
            thisFace.physicalTag = -1;    // At the moment physicalTag of an internal face is not defined as it could be inbetween two different elements

            // Find the geometry type of the face from its parent face
            Dune::GeometryType parentGeometry = gridstorage_.element_[thisFace.element1Index].geometryType;
            thisFace.geometryType = ReferenceElements::general(parentGeometry).type(thisFace.element1SubentityIndex, 1);

            FaceMapIterator faceDBIter = domainBoundaryFaceKey2TagMap_.find(iter->first);
            FaceMapIterator faceIBIter = interiorBoundaryFaceKey2TagMap_.find(iter->first);
            bool foundDB = (faceDBIter != domainBoundaryFaceKey2TagMap_.end());
            bool foundIB = (faceIBIter != interiorBoundaryFaceKey2TagMap_.end());

            assert(!(foundDB && foundIB));  // The same boundary can not be an interior and domain boundary at the same time
            assert((!foundDB) || (connectedFaceInfo.size() == 2));  // If this face is identified as DB, it should be connected to only one element

            std::string boundaryName  = "";
            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: Added Face";
            log_stream << " LocalIndex=" << localFaceIndex;
            log_stream << " AssociatedElementIndex=" << thisFace.element1Index;
            log_stream << " InternalSubentityIndex=" << thisFace.element1SubentityIndex;
            //log_stream << " Order=" << order;
            log_stream << " PhysicalTag=" << faceDBIter->second;

            // Treat domain boundaries
            // NOTE: This will also include periodic boundaries, which will be specialised from domain boundaries later
            if (foundDB) {
            	dbCount++;
            	thisFace.physicalTag = faceDBIter->second;  // Here physical tag is very important as it need not match the tag of the element
            	thisFace.boundaryType = GridStorageType::FaceBoundaryType::DomainBoundary;  // !! When periodic and internal boundaries are introduced this line will change
            	thisFace.element2Index = -1;                  // Boundary Segments do not have a 2nd neighbor
            	thisFace.element2SubentityIndex = - 1; // Boundary Segments do not have a 2nd neighbor
            	thisFace.ptype = Dune::PartitionType::InteriorEntity; // Note: By Dune classification DB are interior entities (see Dune-Grid-Howto manual)

                // Store face in a Map (key -> faceIndex) for constructor purposes
                // Also create a domain boundary index for future indexing
                LocalIndexType localFaceDBIndex = gridstorage_.boundarySegmentIndexMap_.size();
                domainBoundaryFaceKey2LocalIndexMap_[faceDBIter->first] = localFaceIndex;
                gridstorage_.boundarySegmentIndexMap_[localFaceIndex] = localFaceDBIndex;
                gridstorage_.boundarySegment2LocalIndexMap_[localFaceDBIndex] = localFaceIndex;

            	boundaryName += " BoundaryType=domainBoundary";

            } else { // Treat regular interior faces, interior boundaries, and periodic boundaries

             	/* Note: We do not actually store IB as boundary segments at the moment, because there is no particular need for it
             	 * At the moment the IB are in no way different to regular faces, other than having a special physicalTag */
            	 if (foundIB) {
                 	ibCount++;
                 	thisFace.physicalTag = faceIBIter->second;
                 	thisFace.boundaryType = GridStorageType::FaceBoundaryType::InteriorBoundary;  // !! When periodic and internal boundaries are introduced this line will change
                     boundaryName += " BoundaryType=interiorBoundary";
            	 } else {
            		 // If the face is not specifically marked as DBSegment or IBSegment, it is a regular face - interior or PB
            		 thisFace.boundaryType = GridStorageType::FaceBoundaryType::None;
            	 }


                 if (connectedFaceInfo.size() == 2)		{
                     processBoundaryFaceKey2LocalIndexMap_[(*iter).first] = localFaceIndex;  // Store Map (key -> faceIndex)
                     thisFace.ptype = Dune::PartitionType::BorderEntity;
                     thisFace.element2Index = -1;                 // Eventually this will be the Ghost Element Index
                     thisFace.element2SubentityIndex = - 1; //Eventually this will be the subentity index in Ghost

                     // Add this face to the process boundary map
                     LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].size();
                     gridstorage_.processBoundaryIndexMap_[FACE_CODIM][localFaceIndex] = thisFaceLocalPBIndex;

                     boundaryName += " PartitionType=processBoundary";
                 } else {
             		// If it exists, the 2nd neighbouring element needs to be mapped
                 	LocalIndexType    thisAssociatedElement2Index = connectedFaceInfo[2];
                 	InternalIndexType thisFaceSubentityIndex2 = connectedFaceInfo[3];
                 	gridstorage_.elementSubentityCodim1_[thisAssociatedElement2Index][thisFaceSubentityIndex2] = localFaceIndex;
                 	thisFace.element2Index = thisAssociatedElement2Index;              // This is the 2nd neighbor of this internal face
                 	thisFace.element2SubentityIndex = thisFaceSubentityIndex2;
                 	thisFace.ptype = Dune::PartitionType::InteriorEntity;

                 	log_stream << " AssociatedElement2Index=" << thisAssociatedElement2Index;
                 	log_stream << " InternalSubentityIndex2=" << thisFaceSubentityIndex2;

                     // Add this face to the internal map
                     internalFaceKey2LocalIndexMap_[(*iter).first] = localFaceIndex;    // Store Map (key -> faceIndex)
                 	 boundaryName += " PartitionType=interior";
                 }
            }
            log_stream << boundaryName;

            // Add face to the mesh
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
            gridstorage_.face_.push_back(thisFace);

            // Update neighbor index storage size
            int nPBFace = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].size();
            gridstorage_.PB2PBNeighborRank_[FACE_CODIM].resize(nPBFace);
        }


        // Check that all boundaries have been exactly mapped
        assert(dbCount == domainBoundaryFaceKey2TagMap_.size());
        assert(ibCount == interiorBoundaryFaceKey2TagMap_.size());

        // Save space
        domainBoundaryFaceKey2TagMap_.clear();
        interiorBoundaryFaceKey2TagMap_.clear();

        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished generating faces");
    }


    /** Mark correct structural type for all boundary vertices, by making it equal to the
     * structural type of the containing face. When sharing several structural types, the
     * priority order is ProcessBoundary > DomainBoundary > Internal
     *
     * Then, give all process boundary corners a unique index by storing them in the
     * processBoundaryIndexMap_[VERTEX_CODIM]. These are later used for communication
     *
     * NOTE: VERTICES CAN NO LONGER BE DOMAIN BOUNDARIES. ONLY PROCESS BOUNDARIES, INTERNALS, AND GHOSTS
     *
     * [TODO] To make work with arbitrary geometries, must replace numbers, as well as make FaceKeys more abstract
     *
     * */
    void markBorderVertex()
    {
    	// Mark process boundary vertices - overwrite domain boundary vertices where necessary
        // ********************************************************
        int faceCount = 0;
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	LoggingMessage::writePatience("Marking process boundary vertices...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

        	EntityStorage thisFace = gridbase_.entity().data(1, (*faceIter).second);

            for (unsigned int i = 0; i < thisFace.vertexIndexSet.size(); i++)
            {
            	LocalIndexType thisVertexLocalIndex = thisFace.vertexIndexSet[i];
            	gridstorage_.point_[thisVertexLocalIndex].ptype = Dune::PartitionType::BorderEntity;
            }
        }

        // Construct the set of process boundary corners - corners necessary to make process boundary faces on this process
        // ********************************************************
        faceCount = 0;
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	LoggingMessage::writePatience("Generating Boundary Corners...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;

            GlobalIndexType thisVertexKey[3] = {thisFaceKey.node0, thisFaceKey.node1, thisFaceKey.node2};

            // [FIXME] Use ref elem subentity size
            for (int i = 0; i < 3; i++)
            {
            	LocalIndexType thisCornerLocalIndex = gridstorage_.entityIndexMap_[VERTEX_CODIM][thisVertexKey[i]];

            	// For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].find(thisCornerLocalIndex) ==
                	gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].end())
                {
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    LocalIndexType processBoundaryCornerIndex = gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].size();
                	gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM][thisCornerLocalIndex] = processBoundaryCornerIndex;

                    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Adding boundary corner GlobalIndex=" + std::to_string(thisVertexKey[i]));
                }
            }
        }

        // Resize the neighbor rank vector such that it can store for each process boundary corner
        gridstorage_.PB2PBNeighborRank_[VERTEX_CODIM].resize(gridstorage_.processBoundaryIndexMap_[VERTEX_CODIM].size());
    }


    /** Mark correct structural type for all boundary edges, by making it equal to the
     * structural type of the containing face. When sharing several structural types, the
     * priority order is ProcessBoundary > DomainBoundary > Internal
     *
     * Then, give all process boundary edges a unique index by storing them in the
     * processBoundaryIndexMap_[EDGE_CODIM]. These are later used for communication
     *
     * NOTE: EDGES CAN NO LONGER BE DOMAIN BOUNDARIES. ONLY PROCESS BOUNDARIES, INTERNALS, AND GHOSTS
     *
     * [TODO] To make work with arbitrary geometries, must replace numbers, as well as make FaceKeys more abstract
     *
     * */
    void markBorderEdge()
    {
        // Construct the set of EdgeKeys corresponding to edges of processBoundaries
        // ********************************************************
        int faceCount = 0;
        for (FaceMapIterator faceIter = processBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != processBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)
        {
        	LoggingMessage::writePatience("Marking process boundary edges...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

            // Get info of this face wrt associated element
            FaceKey thisFaceKey = (*faceIter).first;
            LocalIndexType thisFaceLocalIndex = (*faceIter).second;
            LocalIndexType assocElementIndex = gridstorage_.face_[thisFaceLocalIndex].element1Index;
            InternalIndexType thisFaceSubIndex = gridstorage_.face_[thisFaceLocalIndex].element1SubentityIndex;
            Dune::GeometryType assocElemGT = gridstorage_.element_[assocElementIndex].geometryType;

            // Get global indices of the associated vertices from the map
            EdgeKey thisEdgeKey[3];
            thisEdgeKey[0].node0 = thisFaceKey.node0;  thisEdgeKey[0].node1 = thisFaceKey.node1;
            thisEdgeKey[1].node0 = thisFaceKey.node0;  thisEdgeKey[1].node1 = thisFaceKey.node2;
            thisEdgeKey[2].node0 = thisFaceKey.node1;  thisEdgeKey[2].node1 = thisFaceKey.node2;

            // [FIXME] Use ref elem subentity size
            for (int i = 0; i < 3; i++)
            {
            	LocalIndexType thisEdgeLocalIndex = edgeKey2LocalIndexMap_[thisEdgeKey[i]];

            	// Mark this edge as subentity of associated boundary internal element
            	// It is strongly suspected that this operation is unnecessary, because it is already done at generateEdges()
            	//InternalIndexType thisEdgeSubIndex = ReferenceElements::general(assocElemGT).subEntity(thisFaceSubIndex, FACE_CODIM, i, EDGE_CODIM);
            	//gridstorage_.elementSubentityCodim2_[assocElementIndex][thisEdgeSubIndex] = thisEdgeLocalIndex;

                // For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].find(thisEdgeLocalIndex) ==
                	gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end())
                {
                	// Change the structural type of this edge to be ProcessBoundary
                	gridstorage_.edge_[thisEdgeLocalIndex].ptype = Dune::PartitionType::BorderEntity;

                	// Insert this edge into the local2local map
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    LocalIndexType processBoundaryEdgeIndex = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size();
                	gridstorage_.processBoundaryIndexMap_[EDGE_CODIM][thisEdgeLocalIndex] = processBoundaryEdgeIndex;

                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: -- From face index= " << thisFaceLocalIndex << " marking process boundary edge index=" << thisEdgeLocalIndex << " EdgeKey= (" << thisEdgeKey[i].node0 << ", " << thisEdgeKey[i].node1 << ")";
                    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
                }
            }
        }

        // Resize the neighbor rank vector such that it can store for each process boundary corner
        gridstorage_.PB2PBNeighborRank_[EDGE_CODIM].resize(gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].size());
    }


    /** Construct OCTree for locating tetrahedrons in mesh
     *
     * Algorithm:
     * 0) Uses OctreeNodes, structures that point to the original grid and contain the localIndex to the associated element.
     * 1) For each element of the mesh, create an OCTreeNode and add it to the octree
     * 1.1) When OctreeNode is created, it computes the rectangular box in which the element fits.
     * 1.2) When OctreeNode is added, it is recursively sinks into a tree starting with the root Octant.
     * 1.3) The root Octant contains the whole mesh of this process. Each octant contains 8 children by spitting itself into 8 smaller boxes.
     * 1.4) The element is fit into a given child octant if its bounding box intersects with that of the child. Thus the element may and will likely be contained inside several octants
     *
     * [FIXME] Revisit OCTree. See how the depth is regulated. See whether all candidates for this element from all octants are returned upon request
     *
     * [TODO]  Use standard logging message
     * [TODO]  Original octree has diagnostics output under #if 0, can append at later stage
     * [TODO]  Replace OCTree pointer by just an instance
     *
     * */
    void constructOctree() {
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started OCTree construction");

        // bounding box of whole mesh
        GlobalCoordinate center = gridstorage_.boundingBoxCenter_;
        GlobalCoordinate extent = gridstorage_.boundingBoxExtent_;

        // octree length is the largest component of extent
        double length = extent[0];
        if (extent[1] > length)  { length = extent[1]; }
        if (extent[2] > length)  { length = extent[2]; }

        // construct LooseOctree with large max depth
        gridstorage_.octree_ = new CurvilinearLooseOctree(center, length, 100, mpihelper_);

        // loop over all tets and insert them in the octree
        for (unsigned int iElem = 0; iElem < gridstorage_.element_.size(); iElem++)
        {
        	LoggingMessage::writePatience("Filling OCTree with elements...", iElem, gridstorage_.element_.size());

            NodeType* thisNode = new NodeType(gridbase_, iElem);
            gridstorage_.octree_->addNode(thisNode);
        }


        int maxDepth, nOctant, nNode;
        double avgNodeDepth;
        gridstorage_.octree_->statistics(maxDepth, avgNodeDepth, nOctant, nNode);

        std::stringstream outputString;
        outputString << "CurvilinearGridConstructor: Constructed OCTree MaxDepth=" << maxDepth;
        outputString << ", #octants=" << nOctant;
        outputString << ", #nodes=" << nNode;
        outputString << ", avg. node depth=" << avgNodeDepth;
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, outputString.str());
    }



private: // Private members

    // Temporary maps necessary to locate and communicate entities during grid base construction
    EdgeKey2EdgeIndexMap edgeKey2LocalIndexMap_;                    // (global edgeKey -> edge_ index)
    FaceKey2FaceIndexMap internalFaceKey2LocalIndexMap_;            // (global faceKey -> gridstorage_.face_ index)
    FaceKey2FaceIndexMap domainBoundaryFaceKey2LocalIndexMap_;      // (global faceKey -> gridstorage_.face_ index)
    FaceKey2FaceIndexMap processBoundaryFaceKey2LocalIndexMap_;     // (global faceKey -> gridstorage_.face_ index)

    FaceKey2FaceIndexMap domainBoundaryFaceKey2TagMap_;      // (global faceKey -> face physical tag)
    FaceKey2FaceIndexMap interiorBoundaryFaceKey2TagMap_;      // (global faceKey -> face physical tag

    // Curvilinear Grid Storage Class
    GridStorageType & gridstorage_;

    // Reference to Curvilinear Grid Base - necessary for OCTree construction
    GridBaseType & gridbase_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDCONSTRUCTOR_HH
