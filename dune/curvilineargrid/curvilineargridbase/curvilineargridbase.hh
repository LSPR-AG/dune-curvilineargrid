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

#ifndef DUNE_CURVILINEARGRIDBASE_HH
#define DUNE_CURVILINEARGRIDBASE_HH

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

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridconstructor.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>


/* ***************************************************************************
 * Specifications: Parallel Curvilinear Mesh Manager
 *
 * Existent Functionality:
 *     - All vertices, edges, faces, elements are accessible through their globalIndex and as subentities
 *     - All elements and faces possess physicalTag - integer associated to their material property
 *     - InterProcessBoundaries automatically generate
 *     - GhostElements automatically communicated and generated (optional)
 *     - GlobalIndex automatically generated for edges, faces and elements
 *     - OCTree functionality allows to find element and its local coordinate corresponding to a given global coordinate
 *     - Supports non-uniformly p-refinemed meshes (not tested)
 *
 *
 * Missing functionality:
 *  - [DESIGN] reader must provide globalIds for vertices and elements, they are not generated automatically.
 *  - [DESIGN] reader must provide all BoundarySegments (process boundaries), they are not communicated automatically
 *  - [TODO] Does NOT support refinement - it is not possible to dynamically refine or coarsen the mesh at the moment
 *  - [TODO] Does NOT support hanging nodes - it is not possible to load non-uniform h-refined mesh at the moment
 *  - [TODO] Does NOT support overlapping elements at the moment
 *  - [TODO] Does NOT support non-tetrahedral meshes. Generalization to arbitrary and mixed geometry meshes is possible but will be cumbersome
 *
 *
 *
 *
 * Development log
 *  - [FIXME] Global indices for vertices must start at 0. GMSH returns them at 1. Must edit GMSHReader to lower globalIndex of vertices and assoc elem/belem
 *
 *  - [FIXME] Need to add normal and outerNormal
 *  - [FIXME] Need to wrap for Dune
 *  - [FIXME] Need to match Dune's internal subentity id convention
 *  - [FIXME] When returning Ghost elements, must check if they are defined, and throw error if not
 *
 * Usage:
 *  - [TODO] Disable all the vertex2string output for multiprocessor case - too much output
 *
 *  ***************************************************************************/




namespace Dune {



template <class ct>
class CurvilinearGridBase {
public:

    /* Stage Structure
     * *******************************************************************/
	struct Stage
	{
		enum
		{
			GRID_CONSTRUCTION,
			GRID_OPERATION
		};
	};



    /* public types
     * *******************************************************************/
    typedef Dune::CurvilinearGridBase<ct>          GridBaseType;
    typedef Dune::CurvilinearGridStorage<ct>       GridStorageType;
    typedef Dune::CurvilinearGridConstructor<ct>   GridConstructorType;

    typedef typename GridStorageType::Vertex                 Vertex;
    typedef typename GridStorageType::VertexStorage          VertexStorage;
    typedef typename GridStorageType::EdgeStorage            EdgeStorage;
    typedef typename GridStorageType::FaceStorage            FaceStorage;
    typedef typename GridStorageType::EntityStorage          EntityStorage;

    typedef typename GridStorageType::EdgeKey                EdgeKey;
    typedef typename GridStorageType::FaceKey                FaceKey;

    typedef typename GridStorageType::EdgeKey2EdgeIdMap      EdgeKey2EdgeIdMap;
    typedef typename GridStorageType::FaceKey2FaceIdMap      FaceKey2FaceIdMap;
    typedef typename GridStorageType::Index2IndexMap         Index2IndexMap;

    typedef typename EdgeKey2EdgeIdMap::iterator             EdgeMapIterator;
    typedef typename FaceKey2FaceIdMap::iterator             FaceMapIterator;
    typedef typename Index2IndexMap::iterator                IndexMapIterator;

    typedef typename GridStorageType::EdgeGeometry           EdgeGeometry;
    typedef typename GridStorageType::FaceGeometry           FaceGeometry;
    typedef typename GridStorageType::ElementGeometry        ElementGeometry;

    typedef Dune::CurvilinearOctreeNode<ct>                       NodeType;
    typedef Dune::CurvilinearLooseOctree<ct, 3, NodeType>         CurvilinearLooseOctree;

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridBase(bool withGhostElements, bool verbose, bool processVerbose, MPIHelper &mpihelper ) :
        withGhostElements_(withGhostElements),
        verbose_(verbose),
        processVerbose_(processVerbose),
        gridstage_(0),
        mpihelper_(mpihelper),
        gridstorage_(),
        gridconstructor_(withGhostElements, verbose, processVerbose, gridstorage_, *this, mpihelper)
    {
        //assert(instance_ == 0);
        //instance_ = this;

        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        std::string log_string = "Initialized CurvilinearGridBase withGhostElements=" + std::to_string(withGhostElements);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
    }

private:
    /** Copy constructor: private, undefined: disallow copy */
    CurvilinearGridBase(const CurvilinearGridBase&);

    /** Assignment operator: private, undefined: disallow assignment */
    CurvilinearGridBase& operator=(const CurvilinearGridBase&);

public:

    /* ***************************************************************************
     * Section: Loading the mesh
     * ***************************************************************************/


    /** \brief Add a new vertex to the mesh
     * \param[in] globalIndex      global index of this vertex
     * \param[in] p                coordinate of this vertex
     * */
    void insertVertex(Vertex p, int globalIndex)
    {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridconstructor_.insertVertex(p, globalIndex);
    }

    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] globalId         the global index of this element
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * */
    void insertElement(Dune::GeometryType gt, int globalId, const std::vector<int> & vertexIndexSet, int order, int physicalTag)
    {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridconstructor_.insertElement(gt, globalId, vertexIndexSet, order, physicalTag);
    }

    /** Insert a boundary segment into the mesh
     *
     *     Note: It is expected that all faces - domain an process boundaries - are inserted by the factory before finalising
     *     Note: Only domain boundary faces have initial globalId given by GMSH. Therefore, we ignore it, and generate our own
     *     globalId for all faces at a later stage.
     *
     *  \param[in] gt                       geometry type of the face (should be a triangle)
     *  \param[in] globalId                 id of the element this face is associated to
     *  \param[in] associatedElementIndex   local index of the element this face is associated to
     *  \param[in] vertexIndexSet           local indices of the interpolatory vertices of this face
     *  \param[in] order                    interpolatory order of the face
     *  \param[in] physicalTag              physical tag of the element (material property)
     *
     * */

    void insertBoundarySegment(Dune::GeometryType gt, int globalId, int associatedElementIndex, const std::vector<int> & vertexIndexSet, int order, int physicalTag)
    {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridconstructor_.insertBoundarySegment(gt, globalId, associatedElementIndex, vertexIndexSet, order, physicalTag);
    }


    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    /** Calls the subroutines for transforming the inserted data into a functional mesh.
     *  Note: It is expected, that all necessary data (vertices, elements and boundary segments) have been added before this function is called.
     * */

    void generateMesh(int nVertexTotalMesh, int nElementTotalMesh) {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridstage_ = Stage::GRID_OPERATION;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Initializing mesh");

        gridconstructor_.generateMesh(nVertexTotalMesh, nElementTotalMesh);

        // Diagnostics output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridBase: Constructed Mesh ";
        log_stream << " nVertexPerMesh="             << nVertexTotal();
        log_stream << " nEdgePerMesh="               << nEdgeTotal();
        log_stream << " nFacePerMesh="               << nFaceTotal();
        log_stream << " nElementPerMesh="            << nElementTotal();
        log_stream << " nVertex="                    << nVertex();
        log_stream << " nEdge="                      << nEdge();
        log_stream << " nFace="                      << nFace();
        log_stream << " nFaceInternal="              << nFaceInternal();
        log_stream << " nFaceProcessBoundary="       << nFaceProcessBoundary();
        log_stream << " nFaceDomainBoundary="        << nFaceDomainBoundary();
        log_stream << " nElement="                   << nElement();
        log_stream << " nGhostElement="              << nGhost();

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
    }

    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    /** Get total number of entities in a mesh  */
    int nVertexTotal() const   { return gridstorage_.nVertexTotal_; }
    int nEdgeTotal() const     { return gridstorage_.nEdgeTotal_; }
    int nFaceTotal() const     { return gridstorage_.nFaceTotal_; }
    int nElementTotal() const  { return gridstorage_.nElementTotal_; }

    /** Get total number of entities on this process  */
    int nVertex() const   { return gridstorage_.point_.size(); }
    int nEdge() const     { return gridstorage_.edge_.size(); }
    int nFace() const     { return gridstorage_.face_.size(); }
    int nElement() const  { return gridstorage_.element_.size(); }
    int nGhost() const    { return gridstorage_.ghostElement_.size(); }

    int nFaceInternal() const        { return gridstorage_.internalFaceIndex_.size(); }
    int nFaceProcessBoundary() const { return gridstorage_.processBoundaryFaceIndex_.size(); }
    int nFaceDomainBoundary() const  { return gridstorage_.domainBoundaryFaceIndex_.size(); }


    /** Vertex coordinate
     *  \param[in] localIndex            local vertex index (insertion index)
     * */
    Vertex vertex(int localIndex) const { return gridstorage_.point_[localIndex].coord; }

    int facePhysicalTag(int localIndex) const           { return gridstorage_.face_[localIndex].physicalTag;}
    int elementPhysicalTag(int localIndex) const        { return gridstorage_.element_[localIndex].physicalTag;}
    int ghostElementPhysicalTag(int localIndex) const   { return gridstorage_.ghostElement_[localIndex].physicalTag;}
    int faceStructuralType(int localIndex) const        { return gridstorage_.face_[localIndex].structuralType; }


    /** Returns the local index of a subentity of a given entity
     *  \param[in] entityIndex              local index of the entity
     *  \param[in] codim                    codimension of the entity
     *  \param[in] subcodim                 codimension of the subentity
     *  \param[in] subentityInternalIndex   subentity internal index wrt entity
     *
     *
     *  Algorithm:
     *  1) Find the parent tetrahedron index from (entityIndex, codim)
     *  2) Find the requested subentity internal index wrt parent tetrahedron
     *  3) Access tetrahedral subentity index storage to return stored answer
     *
     *  Note: To comply with DUNE regulations, if index of a vertex are requested,
     *  that index will be over the corners of the element as if it was a linear element.
     *  Curvilinear interpolation vertices will not be accessible by this method
     *
     * */
    int subentityIndex (int entityIndex, int codim, int subcodim, int subentityInternalIndex)
    {
    	// Stage 1) Find this subentity as an element subentity
    	// **************************************************************************
    	Dune::GeometryType tetrahedronGeometry;
    	tetrahedronGeometry.makeTetrahedron();
    	Dune::ReferenceElement<ct,3> & thisRefElement = Dune::ReferenceElements<ct,3>::general(tetrahedronGeometry);

    	if (subcodim >= codim) { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair"); }

    	int elementLocalId;
    	int elementSubentityInternalIndex1;

    	switch (codim)
    	{
    		case 0 :
    		{
    			elementLocalId = entityIndex;
    			elementSubentityInternalIndex1 = subentityInternalIndex;
    		} break;
    		case 1 :
    		{
    			int elementSubentityInternalIndex2 = gridstorage_.face_[entityIndex].element1SubentityIndex;

    			elementLocalId = gridstorage_.face_[entityIndex].element1Index;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    		case 2 :
    		{
    			int elementSubentityInternalIndex2 = gridstorage_.edge[entityIndex].subentityIndex;

    			elementLocalId = gridstorage_.edge_[entityIndex].element1Index;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    	}


    	// Stage 2) Find index of the element subentity
    	// **************************************************************************

    	switch (subcodim)
    	{
    		// Face
    		case 1 :  return gridstorage_.elementSubentityCodim1_[elementLocalId][elementSubentityInternalIndex1];  break;
    		// Edge
    		case 2 :  return gridstorage_.elementSubentityCodim2_[elementLocalId][elementSubentityInternalIndex1];  break;
    		// Corner
    		case 3 :
    		{
    			int interpolationOrder = gridstorage_.element_[elementLocalId].interpOrder;
    			int cornerInternalIndex = Dune::CurvilinearGeometryHelper::cornerID(tetrahedronGeometry, interpolationOrder, elementSubentityInternalIndex1);
    			return gridstorage_.element_[elementLocalId].vertexIndexSet[cornerInternalIndex];
    		} break;
    	}

    }


    /** \brief local index of the element that is neighbor to this face
     *
     * \param[in]    localIndex               local index of this face
     * \param[in]    internalNeighborIndex    {0,1} - determines which of the two neighbors to address
     *
     * \return local index of the element that is neighbor to this face.
     * NOTE: If the neighbor is a GHOST ELEMENT, the returned int will be for the GhostElementLocalIndex
     * and not the ElementIndex
     *
     *
     *
     * Conventions of internalNeighborIndex for face types
     * * For Domain Boundary there is only one neighbor
     * * For Process Boundary 2nd neighbor is always the Ghost Element
     * * For Internal Face there is no convention on order of the neighbors
     * FIXME: Check if this contradicts Dune convention in any way
     *
     *  */
    int faceNeighbor(int localIndex, int internalNeighborIndex) const
    {
        int rez;

        switch(internalNeighborIndex)
        {
        case 0 : rez = gridstorage_.face_[localIndex].element1Index;  break;
        case 1 : rez = gridstorage_.face_[localIndex].element2Index;  break;
        default: DUNE_THROW(Dune::IOError, "CurvilinearGrid: faceNeighbor() unexpected neighbor index");  break;
        }

        return rez;
    }


    /** Storage data related to this edge, except of explicit vertex coordinates
     *  \param[in] localIndex            local edge index
     *
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    EntityStorage edgeData(int localIndex) const
    {
    	EdgeStorage & thisEdgeData =    gridstorage_.edge_[localIndex];
        EntityStorage & assocElement =  gridstorage_.element_[thisEdgeData.elementIndex];

        EntityStorage thisEdge;
        thisEdge.geometryType.makeLine();
        thisEdge.globalIndex = thisEdgeData.globalIndex;
        thisEdge.interpOrder = assocElement.interpOrder;
        thisEdge.physicalTag = -1;        // Note: Edges do not have a physical tag

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<int> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, 3>(assocElement.geometryType, thisEdge.interpOrder, 2, thisEdgeData.subentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisEdge;

    }


    /** Storage data related to this face, except of explicit vertex coordinates
     *  \param[in] localIndex            local face index
     *
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    EntityStorage faceData(int localIndex) const
    {
    	FaceStorage & thisFaceData = gridstorage_.face_[localIndex];
        EntityStorage & assocElement = gridstorage_.element_[thisFaceData.element1Index];

        EntityStorage thisFace;
        thisFace.geometryType.makeTriangle();
        thisFace.globalIndex = thisFaceData.globalIndex;
        thisFace.interpOrder = assocElement.interpOrder;
        thisFace.physicalTag = thisFaceData.physicalTag;

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<int> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, 3>(assocElement.geometryType, thisFace.interpOrder, 1, thisFaceData.element1SubentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisFace;
    }


    /** Storage data related to this element, except of explicit vertex coordinates
     *  \param[in] localIndex            local element index
     *
     * */
    EntityStorage elementData(int localIndex) const { return gridstorage_.element_[localIndex]; }


    /** Storage data related to this element, except of explicit vertex coordinates
     *  \param[in] localIndex            local ghost element index
     *
     * */
    EntityStorage ghostElementData(int localIndex) const { return gridstorage_.ghostElement_[localIndex];}

    // Construct and retrieve geometry class associated to this edge
    EdgeGeometry edgeGeometry(int localIndex) const
    {
        EntityStorage edgeData = edgeData(localIndex);
        return entityGeometry<1>(edgeData.geometryType, edgeData.vertexIndexSet, edgeData.interpOrder);
    }

    // Construct and retrieve geometry class associated to this face
    FaceGeometry faceGeometry(int localIndex) const
    {
        EntityStorage faceData = faceData(localIndex);
        return entityGeometry<2>(faceData.geometryType, faceData.vertexIndexSet, faceData.interpOrder);
    }

    // Construct and retrieve geometry class associated to this element
    ElementGeometry elementGeometry(int localIndex) const
    {
        EntityStorage thisElementData = elementData(localIndex);
        return entityGeometry<3>(thisElementData.geometryType, thisElementData.vertexIndexSet, thisElementData.interpOrder);
    }

    /** Return pointer to Octree or 0 if it is not constructed. */
    const CurvilinearLooseOctree & octree() const { return *gridstorage_.octree_; }

    void processBoundingBox(Vertex & center, Vertex & extent) const
    {
        center = gridstorage_.boundingBoxCenter_;
        extent = gridstorage_.boundingBoxExtent_;
    }


    /** Searches for elements containing this global coordinate
     *
     * \param[in]    globalC         the global coordinate of the point looked for
     * \param[in]    containerIndex  returns all element indices in which the point was found.
     * In principle there can be more than one if the point is close to the boundary between elements.
     * \param[in]    localC          returns local coordinates associated with this point within all elements it is found in
     *
     * \returns      Whether the point found at all on this process. If not, return arrays should be empty
     *
     *
     * Algorithm:
     *
     * 1) Check if this point is inside the bounding box
     * 2) If yes, check if this point is inside OCTree
     * 3) OCTree returns a list of candidate elements - those located in the lowest level octant the point is found in
     * 4) Loop over candidate elements, use CurvilinearGeometry search functionality
     * 5) If the point is found inside by CurvilinearGeometry, it also returns the associated local coordinate at no extra cost
     * 6) All the elements in which the point is found, as well as the associated local coordinates are returned
     *
     * TODO: PBE file allows femaxx to count time. Use alternative in Dune?
     * */

    bool locateCoordinate(const Vertex & globalC, std::vector<int> & containerIndex, std::vector<Vertex> & localC) const {

        // If the point is not even in the process bounding box, then the coordinate is definitely not on this process
        if (isInsideProcessBoundingBoxGracious(globalC))
        {
            //pbe_start(132, "Octree traversal");

            // Get list of indices of elements which may contain the coordinate
            // This corresponds to the elements at the lowest level of the OCTree
            // Note that no internal element search is performed by the OCTree itself
            int nNodeVisited = 0;
            std::vector<int> elementIndices;
            gridstorage_.octree_->findNode(globalC, elementIndices, nNodeVisited, &isInsideProcessBoundingBoxGracious);



            // Loop over candidate elements and search for local coordinate using internal Geometry mechanism
            for (int i = 0; i < elementIndices.size(); i++) {
                ElementGeometry thisGeometry = elementGeometry(elementIndices[i]);

                Vertex thisLocalC;
                bool isInside = thisGeometry.local( globalC, thisLocalC );

                // If the point is found inside the geometry, add the element index and the found local coordinate to the output
                if (isInside)
                {
                    containerIndex.push_back(elementIndices[i]);
                    localC.push_back(thisLocalC);
                }
            }

            // pbe_stop(132);
            // rDebug("find_tets_by_point: nof_found=%d, nNodeVisited=%d", static_cast<int>(nodes.size()), nNodeVisited);

            return (containerIndex.size() > 0);
        } else {
            return false;
        }
    }

    // Iterators over local indices of the mesh
    // NOTE: There are no iterators over entities because there is no entity object in the core mesh
    // There will be generic entity in the wrapper because the wrapper will define an entity object
    IndexMapIterator elementIndexBegin()  { return gridstorage_.elementGlobal2LocalMap_.begin(); }
    IndexMapIterator faceIndexBegin()     { return gridstorage_.faceGlobal2LocalMap_.begin(); }
    IndexMapIterator edgeIndexBegin()     { return gridstorage_.edgeGlobal2LocalMap_.begin(); }
    IndexMapIterator vertexIndexBegin()   { return gridstorage_.vertexGlobal2LocalMap_.begin(); }

    IndexMapIterator elementIndexEnd()  { return gridstorage_.elementGlobal2LocalMap_.end(); }
    IndexMapIterator faceIndexEnd()     { return gridstorage_.faceGlobal2LocalMap_.end(); }
    IndexMapIterator edgeIndexEnd()     { return gridstorage_.edgeGlobal2LocalMap_.end(); }
    IndexMapIterator vertexIndexEnd()   { return gridstorage_.vertexGlobal2LocalMap_.end(); }

    // This construction allows fast iteration over faces of different type
    std::vector<int>::iterator faceInternalBegin()         { return gridstorage_.internalFaceIndex_.begin(); }
    std::vector<int>::iterator faceDomainBoundaryBegin()   { return gridstorage_.domainBoundaryFaceIndex_.begin(); }
    std::vector<int>::iterator faceProcessBoundaryBegin()  { return gridstorage_.processBoundaryFaceIndex_.begin(); }

    std::vector<int>::iterator faceInternalEnd()         { return gridstorage_.internalFaceIndex_.end(); }
    std::vector<int>::iterator faceDomainBoundaryEnd()   { return gridstorage_.domainBoundaryFaceIndex_.end(); }
    std::vector<int>::iterator faceProcessBoundaryEnd()  { return gridstorage_.processBoundaryFaceIndex_.end(); }



    /** Return pointer to single CurvilinearGridBase instance. */
    //static CurvilinearGridBase* get_instance() { return instance_; }


protected:


    /* ***************************************************************************
     * Section: Auxiliary Methods
     * ***************************************************************************/

    void assertStage(int expectedStage)
    {
    	if ((gridstage_ == Stage::GRID_OPERATION) && (expectedStage == Stage::GRID_CONSTRUCTION)) { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Attempted to insert entities into grid after construction"); }
    }

    // Get curved geometry of an entity
    // TODO: assert mydim == element geometry type dim
    template<int mydim>
    Dune::CurvilinearGeometry<ct, mydim, 3> entityGeometry(
            Dune::GeometryType gt,
            std::vector<int> & vertexIndexSet,
            int order) const
    {
        std::vector<Vertex> entityVertices;
        for (int i = 0; i < vertexIndexSet.size(); i++) { entityVertices.push_back(gridstorage_.point_[vertexIndexSet[i]].coord); }

        return Dune::CurvilinearGeometry<ct, mydim, 3> (gt, entityVertices, order);
    }

    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    // Checks if given point fits into the bounding box of the mesh
    bool isInsideProcessBoundingBoxGracious(const Vertex & point) const {
        const double grace_tolerance = 1e-13;
        Vertex center, extent;

        bool isInside = true;

        isInside &= fabs(gridstorage_.boundingBoxCenter_[0] - point[0]) <= gridstorage_.boundingBoxExtent_[0] * (1.0 + grace_tolerance);
        isInside &= fabs(gridstorage_.boundingBoxCenter_[1] - point[1]) <= gridstorage_.boundingBoxExtent_[1] * (1.0 + grace_tolerance);
        isInside &= fabs(gridstorage_.boundingBoxCenter_[2] - point[2]) <= gridstorage_.boundingBoxExtent_[2] * (1.0 + grace_tolerance);

        return isInside;
    }

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::stringstream tmp_stream;

        int nEntry = V.size();
        if (nEntry == 0)  { tmp_stream << "Null"; }
        for (int i = 0; i < nEntry; i++) {
        	tmp_stream << V[i];
        	if (i != nEntry - 1) { tmp_stream << " "; }
        }
        return tmp_stream.str();
    }

private: // Private members

    bool verbose_;
    bool processVerbose_;
    bool withGhostElements_;


    // The stage of the grid determines if the grid has already been assembled or not
    int gridstage_;

    // Curvilinear Grid Constructor Class
    GridConstructorType gridconstructor_;

    // Curvilinear Grid Storage Class
    GridStorageType gridstorage_;
    
    /** Pointer to single CurvilinearGridBase instance (Singleton) */
    //static CurvilinearGridBase* instance_ = 0;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDBASE_HH
