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
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/loggingtimer.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilinearglobalindexconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearghostconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearpostconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>


namespace Dune {




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

    typedef Dune::CurvilinearGhostConstructor<GridBase>         GridGhostConstructor;
    typedef Dune::CurvilinearPostConstructor<GridBase>		    GridPostConstructor;
    typedef Dune::CurvilinearGlobalIndexConstructor<GridBase>   GridGlobalIndexConstructor;


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


	typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridConstructor(
    		GridStorageType & gridstorage,
    		GridBaseType & gridbase,
    		MPIHelper &mpihelper) :
        gridstorage_(gridstorage),
        gridbase_(gridbase),
        mpihelper_(mpihelper)
    {
        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "Initialized CurvilinearGridConstructor");
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
    void insertVertex(Vertex p, GlobalIndexType globalIndex, Dune::PartitionType ptype = Dune::PartitionType::InteriorEntity)
    {
        VertexStorage point;
        point.coord = p;
        point.globalIndex = globalIndex;
        point.ptype = ptype;

        gridstorage_.entityIndexMap_[VERTEX_CODIM][globalIndex] = gridstorage_.point_.size();
        gridstorage_.point_.push_back(point);

        std::stringstream log_stream;
        log_stream << "CurvilinearGridConstructor: Inserted vertex LocalIndex=" << gridstorage_.point_.size()-1 << " GlobalIndex=" << globalIndex;
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
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
        if (!gt.isTetrahedron() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
        	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertElement() unexpected element type or number of interpolatory points");
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
        log_stream << "CurvilinearGridConstructor: Inserted Element Type=" << Dune::CurvilinearGeometryHelper::geometryName(gt);
        log_stream << " LocalIndex=" << thisLocalIndex;
        log_stream << " Order=" << order;
        log_stream << " PhysicalTag=" << physicalTag;
        log_stream << " VertexIndices=(" << Dune::VectorHelper::vector2string(vertexIndexSet) << ")";
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
    }

    /** Insert a boundary segment into the mesh
     *
     *     Note: It is expected that all faces - domain an process boundaries - are inserted by the factory before finalising
     *     Note: Only domain boundary faces have initial globalId given by GMSH. Therefore, we ignore it, and generate our own
     *     globalId for all faces at a later stage.
     *
     *  \param[in] gt                       geometry type of the face (should be a triangle)
     *  \param[in] associatedElementIndex   local index of the element this face is associated to
     *  \param[in] vertexIndexSet           local indices of the interpolatory vertices of this face
     *  \param[in] order                    interpolatory order of the face
     *  \param[in] physicalTag              physical tag of the element (material property)
     *
     * */

    void insertBoundarySegment(
    	Dune::GeometryType gt,
    	LocalIndexType associatedElementIndex,
    	const std::vector<LocalIndexType> & vertexIndexSet,
    	InterpolatoryOrderType order,
    	PhysicalTagType physicalTag)
    {
        if (!gt.isTriangle() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
        	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
            DUNE_THROW(Dune::IOError, "CurvilinearGridConstructor: insertBoundarySegment() unexpected number of interpolatory points");
        }


        // Get corners of this face
        // **********************************************************************************
        std::vector<LocalIndexType> faceCorners = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 2>(gt, vertexIndexSet, order);
        FaceKey thisFaceKey;
        thisFaceKey.node0 = gridstorage_.point_[faceCorners[0]].globalIndex;
        thisFaceKey.node1 = gridstorage_.point_[faceCorners[1]].globalIndex;
        thisFaceKey.node2 = gridstorage_.point_[faceCorners[2]].globalIndex;

        // Sort in ascending order
        thisFaceKey.sort();


        // Take associated element, get all its corners, get all keys, compare to face key
        // **********************************************************************************
        std::vector<LocalIndexType> elementCorners = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 3>(
        		gridstorage_.element_[associatedElementIndex].geometryType,
        		gridstorage_.element_[associatedElementIndex].vertexIndexSet,
        		gridstorage_.element_[associatedElementIndex].interpOrder
        );



        // Search for the face among subentities of the element
        // **********************************************************************************
        int j = 0;
        int nFacePerTetrahedron = 4; // [FIXME] Replace number by ref elem subentity size
        bool found_face = false;

        while (!found_face)
        {
            if (j == nFacePerTetrahedron)  {
            	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridConstructor: insertBoundarySegment() did not find the face in the associated element");
                DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertBoundarySegment() did not find the face in the associated element");
            }


            // Get internal indices of the corners of this face wrt its associated element
            Dune::GeometryType assocElementGeometryType = gridstorage_.element_[associatedElementIndex].geometryType;
            InternalIndexType node0SubIndex = ReferenceElements::general(assocElementGeometryType).subEntity(j, FACE_CODIM, 0, VERTEX_CODIM);
            InternalIndexType node1SubIndex = ReferenceElements::general(assocElementGeometryType).subEntity(j, FACE_CODIM, 1, VERTEX_CODIM);
            InternalIndexType node2SubIndex = ReferenceElements::general(assocElementGeometryType).subEntity(j, FACE_CODIM, 2, VERTEX_CODIM);

            // Define (key = sorted localIndices of corners)
            FaceKey thisKey;
            thisKey.node0 = gridstorage_.point_[elementCorners[node0SubIndex]].globalIndex;
            thisKey.node1 = gridstorage_.point_[elementCorners[node1SubIndex]].globalIndex;
            thisKey.node2 = gridstorage_.point_[elementCorners[node2SubIndex]].globalIndex;


            // Sort in ascending order
            thisKey.sort();

            // By comparison find internalIndex of this face
            if (thisKey == thisFaceKey)
            {
                found_face = true;

                // Store face in a Map (key -> faceIndex) for constructor purposes
                // Also create a domain boundary index for future indexing
                LocalIndexType localFaceIndex = gridstorage_.face_.size();
                LocalIndexType localFaceDBIndex = gridstorage_.boundarySegmentIndexMap_.size();
                domainBoundaryFaceKey2LocalIndexMap_[thisKey] = localFaceIndex;
                gridstorage_.boundarySegmentIndexMap_[localFaceIndex] = localFaceDBIndex;

                // Store Vector (faceId -> associated element)
                FaceStorage thisFaceAsSubentity;
                thisFaceAsSubentity.geometryType = gt;
                thisFaceAsSubentity.globalIndex  = 0;                  // At this stage the globalId is not known yet
                thisFaceAsSubentity.ptype = Dune::PartitionType::InteriorEntity;
                thisFaceAsSubentity.boundaryType = GridStorageType::FaceBoundaryType::DomainBoundary;  // !! When periodic and internal boundaries are introduced this line will change

                thisFaceAsSubentity.element1Index = associatedElementIndex;
                thisFaceAsSubentity.element2Index = -1;              // Boundary Segments do not have a 2nd neighbor
                thisFaceAsSubentity.element1SubentityIndex = j;
                thisFaceAsSubentity.physicalTag = physicalTag;    // Here physical tag is very important as it need not match the tag of the element

                gridstorage_.face_.push_back(thisFaceAsSubentity);


                std::stringstream log_stream;
                log_stream << "CurvilinearGridConstructor: Inserted BoundarySegment Type=" << Dune::CurvilinearGeometryHelper::geometryName(gt);
                log_stream << " LocalIndex=" << gridstorage_.face_.size()-1;
                log_stream << " Order=" << order;
                log_stream << " PhysicalTag=" << physicalTag;
                log_stream << " AssociatedElementIndex=" << associatedElementIndex;
                log_stream << " InternalSubentityIndex=" << j;
                LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
            }

            j++;
        }
    }


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

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Initializing mesh");

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

        if (size_ > 1)  // Parallel case
        {
#if HAVE_MPI
        	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Parallel Grid");

        	// Generate Global Index
        	// **********************************************************
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


        	// Generate Ghost elements
        	// **********************************************************

            if (gridstorage_.withGhostElements_)
            {
            	LoggingTimer::time("CurvilinearGridConstructor: Ghost Element Constructor");
            	GridGhostConstructor ghostConstructor(gridstorage_, mpihelper_);
            	ghostConstructor.generate();
            	LoggingTimer::time("CurvilinearGridConstructor: Ghost Element Constructor");
            }
#endif
        }
        else  // Serial case
        {
        	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Assembling Serial Grid");
            // Serial case:
            // * Boundary Neighbors not necessary, since all boundaries are domain boundaries
            // * No ghost elements, even if requested by user
            // * Fake globalIndex by making it equal to localIndex

            gridstorage_.withGhostElements_ = false;
            for (unsigned int i = 0; i < gridstorage_.edge_.size();    i++)  { gridstorage_.edge_[i].globalIndex = i;     gridstorage_.entityIndexMap_[EDGE_CODIM][i] = i; }
            for (unsigned int i = 0; i < gridstorage_.face_.size();    i++)  { gridstorage_.face_[i].globalIndex = i;     gridstorage_.entityIndexMap_[FACE_CODIM][i] = i; }
            for (unsigned int i = 0; i < gridstorage_.element_.size(); i++)  { gridstorage_.element_[i].globalIndex = i;  gridstorage_.entityIndexMap_[ELEMENT_CODIM][i] = i;  gridstorage_.entityInternalIndexSet_[ELEMENT_CODIM].insert(i); }

            for (FaceMapIterator faceIter = internalFaceKey2LocalIndexMap_.begin(); faceIter != internalFaceKey2LocalIndexMap_.end(); faceIter++)              { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.entityInternalIndexSet_[FACE_CODIM].insert(localIndex); }
            for (FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.begin(); faceIter != domainBoundaryFaceKey2LocalIndexMap_.end(); faceIter++)  { LocalIndexType localIndex = (*faceIter).second;        gridstorage_.faceDomainBoundaryIndexSet_.insert(localIndex); }

            gridstorage_.nEntityTotal_[EDGE_CODIM] = gridstorage_.edge_.size();
            gridstorage_.nEntityTotal_[FACE_CODIM] = gridstorage_.face_.size();
        }

        // Deletes all temporary memory
        // ************************************************************
        edgeKey2LocalIndexMap_.clear();
        internalFaceKey2LocalIndexMap_.clear();
        domainBoundaryFaceKey2LocalIndexMap_.clear();
        processBoundaryFaceKey2LocalIndexMap_.clear();


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

        	LoggingTimer::time("CurvilinearGridConstructor: Communication of neighbor ranks");
            postConstructor.communicateCommunicationEntityNeighborRanks();
            LoggingTimer::time("CurvilinearGridConstructor: Communication of neighbor ranks");
#endif
        }


        // Construct OCTree
        // ************************************************************

        //LoggingTimer::time("CurvilinearGridConstructor: Constructing OCTree");
        computeProcessBoundingBox();
        //constructOctree();
        //LoggingTimer::time("CurvilinearGridConstructor: Constructing OCTree");
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
        Vertex min = gridstorage_.point_[0].coord;
        Vertex max = min;

        for (unsigned int i = 1; i < gridstorage_.point_.size(); i ++) {
        	Dune::LoggingMessage::writePatience(" Computing process bounding box...", i, gridstorage_.point_.size());

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
        	Dune::LoggingMessage::writePatience("Generating edges...", iElem, nElem);

        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<LocalIndexType> elementCornerLocalIndexSet = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

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
                    LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());

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
        	Dune::LoggingMessage::writePatience("Generating faces...", iElem, nElem);

        	EntityStorage & thisElem = gridstorage_.element_[iElem];
            std::vector<int> elementCornerLocalIndexSet = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ctype, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

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


                FaceMapIterator faceIter = domainBoundaryFaceKey2LocalIndexMap_.find(thisKey);

                // Mark this face for creation if it is not an already existing Domain Boundary
                // Otherwise note its local index
                if (faceIter == domainBoundaryFaceKey2LocalIndexMap_.end())
                {
                    std::vector<int> connectedFaceInfo;
                    tmpMapIterator iter = tmpFaceMap.find(thisKey);

                    if (iter != tmpFaceMap.end()) { connectedFaceInfo = std::vector<int> ( (*iter).second ); }

                    connectedFaceInfo.push_back(iElem);
                    connectedFaceInfo.push_back(iFace);

                    tmpFaceMap[thisKey] = connectedFaceInfo;


                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridConstructor: Adding FaceKey=(" << thisKey.node0 << ", " << thisKey.node1 << ", " << thisKey.node2 << ") attached to total of " << connectedFaceInfo.size() / 2 << " elements";
                    LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
                } else {
                	LocalIndexType localFaceIndex = (*faceIter).second;
                	gridstorage_.elementSubentityCodim1_[iElem][iFace] = localFaceIndex;
                }
            }
        }


        // Add internal and process boundary faces to the mesh
        for (tmpMapIterator iter = tmpFaceMap.begin();  iter != tmpFaceMap.end();  iter++)
        {
            FaceStorage thisFace;
            LocalIndexType localFaceIndex = gridstorage_.face_.size();
            std::vector<int> connectedFaceInfo = (*iter).second;

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


            std::stringstream log_stream;
            log_stream << "CurvilinearGridConstructor: Added Face";
            log_stream << " LocalIndex=" << localFaceIndex;
            log_stream << " AssociatedElementIndex=" << thisFace.element1Index;
            log_stream << " InternalSubentityIndex=" << thisFace.element1SubentityIndex;

            // Store internal, domain and process boundaries separately for faster iterators
            if (connectedFaceInfo.size() == 2)
            {
                thisFace.ptype = Dune::PartitionType::BorderEntity;
                thisFace.boundaryType = GridStorageType::FaceBoundaryType::None;

                processBoundaryFaceKey2LocalIndexMap_[(*iter).first] = localFaceIndex;  // Store Map (key -> faceIndex)
                thisFace.element2Index = 0;                                             // Eventually this will be the Ghost Element Index

                // Add this face to the process boundary map
                LocalIndexType thisFaceLocalPBIndex = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].size();
                gridstorage_.processBoundaryIndexMap_[FACE_CODIM][localFaceIndex] = thisFaceLocalPBIndex;

                log_stream << " StructuralType=processBoundary";
            }
            else
            {
            	// In this case 2nd neighbouring element needs to be mapped
                LocalIndexType    thisAssociatedElement2Index = connectedFaceInfo[2];
                InternalIndexType thisFaceSubentityIndex2 = connectedFaceInfo[3];
                gridstorage_.elementSubentityCodim1_[thisAssociatedElement2Index][thisFaceSubentityIndex2] = localFaceIndex;
                log_stream << " AssociatedElement2Index=" << thisAssociatedElement2Index;
                log_stream << " InternalSubentityIndex2=" << thisFaceSubentityIndex2;

                // Add this face to the internal map
                thisFace.ptype = Dune::PartitionType::InteriorEntity;
                thisFace.boundaryType = GridStorageType::FaceBoundaryType::None;

                internalFaceKey2LocalIndexMap_[(*iter).first] = localFaceIndex;    // Store Map (key -> faceIndex)
                thisFace.element2Index = thisAssociatedElement2Index;              // This is the 2nd neighbor of this internal face
                log_stream << " StructuralType=internal";
            }

            // Add face to the mesh
            LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
            gridstorage_.face_.push_back(thisFace);

            // Update neighbor index storage size
            int nPBFace = gridstorage_.processBoundaryIndexMap_[FACE_CODIM].size();
            gridstorage_.PB2PBNeighborRank_[FACE_CODIM].resize(nPBFace);
        }

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Finished generating faces");
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
        	Dune::LoggingMessage::writePatience("Marking process boundary vertices...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

        	EntityStorage thisFace = gridbase_.entityData(1, (*faceIter).second);

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
        	Dune::LoggingMessage::writePatience("Generating Boundary Corners...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

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

                    LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: -- Adding boundary corner GlobalIndex=" + std::to_string(thisVertexKey[i]));
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
        	Dune::LoggingMessage::writePatience("Marking process boundary edges...", faceCount++, processBoundaryFaceKey2LocalIndexMap_.size());

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
                    LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
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
    	LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, "CurvilinearGridConstructor: Started OCTree construction");

        // bounding box of whole mesh
        Vertex center = gridstorage_.boundingBoxCenter_;
        Vertex extent = gridstorage_.boundingBoxExtent_;

        // octree length is the largest component of extent
        double length = extent[0];
        if (extent[1] > length)  { length = extent[1]; }
        if (extent[2] > length)  { length = extent[2]; }

        // construct LooseOctree with large max depth
        gridstorage_.octree_ = new CurvilinearLooseOctree(center, length, 100, mpihelper_);

        // loop over all tets and insert them in the octree
        for (unsigned int iElem = 0; iElem < gridstorage_.element_.size(); iElem++)
        {
        	Dune::LoggingMessage::writePatience("Filling OCTree with elements...", iElem, gridstorage_.element_.size());

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
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, outputString.str());
    }



private: // Private members

    // Temporary maps necessary to locate and communicate entities during grid base construction
    EdgeKey2EdgeIndexMap edgeKey2LocalIndexMap_;                    // (global edgeKey -> edge_ index)
    FaceKey2FaceIndexMap internalFaceKey2LocalIndexMap_;            // (global faceKey -> gridstorage_.face_ index)
    FaceKey2FaceIndexMap domainBoundaryFaceKey2LocalIndexMap_;      // (global faceKey -> gridstorage_.face_ index)
    FaceKey2FaceIndexMap processBoundaryFaceKey2LocalIndexMap_;     // (global faceKey -> gridstorage_.face_ index)

    // Curvilinear Grid Storage Class
    GridStorageType & gridstorage_;

    // Reference to Curvilinear Grid Base - necessary for OCTree construction
    GridBaseType & gridbase_;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDCONSTRUCTOR_HH
