/***************************************************************************
                          tetmesh.h  -  description
                             -------------------
    begin                : Mon Dec 15 2003
    copyright            : (C) 2003 by Roman Geus
    email                : roman.geus@psi.ch
    edited by            : Hua Guo, Sep 2010
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TETMESH_H
#define TETMESH_H

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>
#include "looseoctree.h"


#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/exceptions.hh>


/* ***************************************************************************
 * Specifications: Parallel Curvilinear Mesh Manager
 *
 * Existent Functionality:
 * 	- All vertices, edges, faces, elements are accessible through their globalIndex and as subentities
 * 	- All elements and faces possess physicalTag - integer associated to their material property
 * 	- InterProcessBoundaries automatically generate
 * 	- GhostElements automatically communicated and generated (optional)
 * 	- GlobalIndex automatically generated for edges and faces
 * 	- OCTree functionality allows to find element and its local coordinate corresponding to a given global coordinate
 * 	- Supports non-uniformly p-refinemed meshes (not tested)
 *
 *
 * Missing functionality:
 *  - [DESIGN] reader must provide globalIds for vertices and elements, they are not generated automatically.
 *  - [DESIGN] reader must provide all BoundarySegments (process boundaries), they are not communicated automatically
 *  - [TODO] Does NOT support refinement - it is not possible to dynamically refine or coarsen the mesh at the moment
 *  - [TODO] Does NOT support hanging nodes - it is not possible to load non-uniform h-refined mesh at the moment
 *  - [TODO] Does NOT support overlapping elements at the moment
 *
 *
 *
 *
 * Development log
 *  - [FIXME] Need to wrap for Dune
 *  - [FIXME] Need to match Dune's internal subentity id convention
 *  - [FIXME] Need to rewrite OCTree so it fits new convention
 *  - [FIXME] Need to introduce parallel logging through verbose as in other classes
 *  - [FIXME] Need to add normal and outerNormal
 *  - [FIXME] Introduce some convention between int and id_t
 *
 *
 *
 *  ***************************************************************************/




namespace Dune {







template <class ct>
struct FemaxxVertexStorage
{
	id_t global_id;
	Dune::FieldVector<ct, 3> coord;
};

struct FemaxxEdgeStorage
{
	id_t global_id;
	id_t element_id;
	int subentityNo;
};

// Face stores indices to 2 intersecting elements, and subentity index for the first one
// Note: element1_id is always an internal element index
// Note: element2_id is an internal element index only for internal faces, it is a ghost element index for process boundaries, and it is -1 for domain boundaries.
struct FemaxxFaceStorage
{
	id_t global_id;
	id_t element1_id;
	id_t element2_id;
	int element1_subentityNo;
	int physicalTag;
};

struct FemaxxEntityStorage
{
	Dune::GeometryType geometryType;
	id_t global_id;
	std::vector<int> vertexIds;
	int interpOrder;
	int physicalTag;
};



// This is a minimal info necessary to recognize an edge among all processes
// The sorted globalId's of vertices of an edge
struct FemaxxEdgeMapKey
{
	int node0;
	int node1;

	void sort() {
		if (node0 > node1)  { std::swap(node0, node1); }
	}
};

// This is a minimal info necessary to recognize a face among all processes
// The sorted globalId's of vertices of an edge
struct FemaxxFaceMapKey
{
	int node0;
	int node1;
	int node2;

	void sort()
	{
		if (node0 > node1)  { std::swap(node0, node1); }
		if (node1 > node2)  { std::swap(node1, node2); }
		if (node0 > node1)  { std::swap(node0, node1); }
	}
};

enum {
	FemaxxFaceInternal = 0,
	FemaxxFaceDomainBoundary = 1,
	FemaxxFaceProcessBoundary = 2
};

template <class ct>
class TetMesh {
public:

	/* public types */
	typedef Dune::FieldVector<ct, 3> Vertex;
	typedef std::vector<Vertex> VertexVector;
	typedef FemaxxVertexStorage<ct> FemaxxPointStorage;

    typedef std::map<FemaxxEdgeMapKey, id_t> EdgeKey2EdgeIdMap;
    typedef std::map<FemaxxFaceMapKey, id_t> FaceKey2FaceIdMap;

    typedef Dune::LagrangeGeometry<ct, 1, 3> EdgeGeometry;
    typedef Dune::LagrangeGeometry<ct, 2, 3> FaceGeometry;
    typedef Dune::LagrangeGeometry<ct, 3, 3> ElementGeometry;

    typedef std::map<id_t, id_t> Index2IndexMap;



public: /* public methods */

    /** Default constructor - DO NOT USE*/
    TetMesh() { }

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    TetMesh(int totalVerticesPerMesh, int totalElementsPerMesh, bool withGhostElements, bool verbose,  MPIHelper &mpihelper ) :
    	_totalVerticesPerMesh(totalVerticesPerMesh),
    	_totalEdgesPerMesh(0),
    	_totalFacesPerMesh(0),
    	_totalElementsPerMesh(totalElementsPerMesh),
    	_withGhostElements(withGhostElements),
    	_verbose(verbose),
    	_octree(0),
    	 _mpihelper(mpihelper)
    {
    	assert(instance_ == 0);
    	instance_ = this;
    }

    /** Destructor */
    ~TetMesh() {
        // delete octree
        if (_octree)  { delete _octree; }
        instance_ = 0;
    }

private:
    /** Copy constructor: private, undefined: disallow copy */
    TetMesh(const TetMesh&);

    /** Assignment operator: private, undefined: disallow assignment */
    TetMesh& operator=(const TetMesh&);

public:

    /* ***************************************************************************
     * Section: Loading the mesh
     * ***************************************************************************/


    /** Reserve space for nofPoint mesh points */
    void init_points(int nofPoint)
    {
        _points.reserve(nofPoint);
    }

    void init_elements(int nofElement)
    {
        _elements.reserve(nofElement);
    }


    /** Add a new point to the mesh */
    void insertVertex(int id, Vertex p)
    {
    	FemaxxPointStorage point;
    	point.coord = p;
    	point.global_id = id;

    	_vertexGlobal2LocalMap[id] = _points.size();
    	_points.push_back(point);

    }

    /** Insert an element into the mesh */
    void insertElement(Dune::GeometryType gt, id_t global_id, std::vector<int> vertexIds, int order, int physicalTag)
    {
        if (!gt.isTetrahedron() || (vertexIds.size() != Dune::CurvilinearElementInterpolator<ct, 2, 3>::dofPerOrder(gt, order)))  {
        	DUNE_THROW(Dune::IOError, "FemaxxGrid: insertElement() unexpected element type or number of interpolatory points");
        }

    	FemaxxEntityStorage thisElement;

    	thisElement.geometryType = gt;
    	thisElement.global_id = global_id;
    	thisElement.interpOrder = order;
    	thisElement.physicalTag = physicalTag;
    	thisElement.vertexIds = vertexIds;

    	_elementGlobal2LocalMap[global_id] = _elements.size();
    	_elements.push_back(thisElement);
    }

    /** Insert an boundary segment into the mesh
     *
     * 	Note: It is expected that all faces - domain an process boundaries - are inserted by the factory before finalising
     * 	Note: Only domain boundary faces have initial globalId given by GMSH. Therefore, we ignore it, and generate our own
     * 	globalId for all faces at a later stage.
     *
     *  \param[in] gt				   	geometry type of the face (should be a triangle)
     *  \param[in] assoc_element_id		id of the element this face is associated to
     *  \param[in] vertexIds			ids of the interpolatory vertices of this face
     *  \param[in] order				interpolatory order of the face
     *  \param[in] physicalTag			physical tag of the element (material property)
     *
     * */
    void insertBoundarySegment(Dune::GeometryType gt, id_t assoc_element_id, std::vector<int> vertexIds, int order, int physicalTag)
    {
        if (!gt.isTriangle() || (vertexIds.size() != Dune::CurvilinearElementInterpolator<ct, 2, 3>::dofPerOrder(gt, order)))  {
        	DUNE_THROW(Dune::IOError, "FemaxxGrid: insertBoundarySegment() unexpected number of interpolatory points");
        }


    	// Get corners of this face
    	// **********************************************************************************
    	std::vector<id_t> faceCorners = getEntityCorners(gt, vertexIds, order);
    	FemaxxFaceMapKey thisFaceKey;
    	thisFaceKey.node0 = _points[faceCorners[0]].global_id;
    	thisFaceKey.node1 = _points[faceCorners[1]].global_id;
    	thisFaceKey.node2 = _points[faceCorners[2]].global_id;

		// Sort in ascending order
    	thisFaceKey.sort();


		// Take associated element, get all its corners, get all keys, compare to face key
		// **********************************************************************************
		std::vector<id_t> elementCorners = getEntityCorners(
				_elements[assoc_element_id].geometryType,
				_elements[assoc_element_id].vertexIds,
				_elements[assoc_element_id].interpOrder
		);
		std::vector<std::vector<int> > elementInternalIndicesFaces = internalFacesOfTetrahedron();



		// Search for the face among subentities of the element
		// **********************************************************************************
		int j = 0;
		bool found_face = false;

		while (!found_face)
		{
	        if (j == elementInternalIndicesFaces.size())  {
	        	DUNE_THROW(Dune::IOError, "FemaxxGrid: insertBoundarySegment() did not find the face in the associated element");
	        }

			// Define (key = sorted localIds of corners)
			FemaxxFaceMapKey thisKey;
			thisKey.node0 = _points[elementCorners[elementInternalIndicesFaces[j][0]]].global_id;
			thisKey.node1 = _points[elementCorners[elementInternalIndicesFaces[j][1]]].global_id;
			thisKey.node2 = _points[elementCorners[elementInternalIndicesFaces[j][2]]].global_id;


			// Sort in ascending order
			thisKey.sort();

			// By comparison find internalIndex of this face
			if (thisKey == thisFaceKey)
			{
				found_face = true;

		    	// Store domain and process boundaries separately for faster iterators
				// Store Map (key -> domainBoundaryIndex), Vector (domainBoundaryIndex -> faceIndex)

				_boundaryInternalMap[thisKey] = _domainBoundaryFaceIndex.size();
				_domainBoundaryFaceIndex.push_back(_faces.size());

				// Store Vector (faceId -> associated element)
				FemaxxFaceStorage thisFaceAsSubentity;
				thisFaceAsSubentity.global_id  = 0;				  // At this stage the globalId is not known yet
				thisFaceAsSubentity.element1_id = assoc_element_id;
				thisFaceAsSubentity.element2_id = -1;			  // Boundary Segments do not have a 2nd neighbor
				thisFaceAsSubentity.element1_subentityNo = j;
				thisFaceAsSubentity.physicalTag = physicalTag;    // Here physical tag is very important as it need not match the tag of the element

				_faces.push_back(thisFaceAsSubentity);
			}

			j++;
		}
    }


    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    /*  */
    void generate_mesh() {
    	typedef FaceKey2FaceIdMap::iterator IterType;

    	// Construct missing parts of the mesh
        generateEdges();
        generateFaces();
        generateProcessBoundaryNeighbors();
        generateGhostElements();
        generateGlobalIndices();
        construct_octree();


    	// Deletes all temporary memory
        _edgemap.clear();

        _internalFaceIndex.reserve(_internalInternalMap.size());
        _domainBoundaryFaceIndex.reserve(_boundaryInternalMap.size());
        _processBoundaryFaceIndex.reserve(_processInternalMap.size());

        for (IterType it = _internalInternalMap.begin(); it < _internalInternalMap.end();  it++) { _internalFaceIndex.push_back((*it).second); }
        for (IterType it = _boundaryInternalMap.begin(); it < _boundaryInternalMap.end();  it++) { _domainBoundaryFaceIndex.push_back((*it).second); }
        for (IterType it = _processInternalMap.begin();  it < _processInternalMap.end();   it++) { _processBoundaryFaceIndex.push_back((*it).second); }

        _internalInternalMap.clear();
        _boundaryInternalMap.clear();
        _processInternalMap.clear();
    }

    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    /* Get total number of entities in a mesh  */
    int getTotalVerticesPerMesh()  { return _totalVerticesPerMesh; }
    int getTotalEdgesPerMesh()     { return _totalEdgesPerMesh; }
    int getTotalFacesPerMesh()     { return _totalFacesPerMesh; }
    int getTotalElementsPerMesh()  { return _totalElementsPerMesh; }

    /* Get total number of entities on this process  */
    int getVertexNumber()  { return _points.size(); }
    int getEdgeNumber()    { return _edges.size(); }
    int getFaceNumber()    { return _faces.size(); }
    int getElementNumber()  { return _elements.size(); }


    // Returns vertex by its localId
    Vertex getVertex(id_t id) const { return _points[id]; }


    // Returns the storage data related to this edge, except of explicit vertex coordinates
    // This data is stored only partially - vertexIds are extracted from associated element
    FemaxxEntityStorage getEdgeData(id_t id)
    {
    	FemaxxEntityStorage & assocElement = _elements[_edges[id].element_id];
    	FemaxxEntityStorage thisEdge;
    	thisEdge.geometryType.makeLine();
    	thisEdge.global_id = _edges[id].global_id;
    	thisEdge.interpOrder = assocElement.interpOrder;
    	thisEdge.physicalTag = -1;		// Note: Edges do not have a physical tag

    	// Get the internal element vertex indices associated with this face as a subentity
    	std::vector<int> subentityVertexIndices = Dune::CurvilinearElementInterpolator<ct, 1, 3>::subentityInternalCoordinates(assocElement.geometryType, thisEdge.interpOrder, 2, _edges[id].subentityNo );

    	// Calculate the localId's of vertices of this face by extracting them from the element vertex Ids
    	for(int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIds.push_back(assocElement.vertexIds[subentityVertexIndices[i]]); }

    	return thisEdge;

    }


    // Returns the storage data related to this edge, except of explicit vertex coordinates
    // This data is stored only partially - vertexIds are extracted from associated element
    FemaxxEntityStorage getFaceData(id_t id)
    {
    	FemaxxEntityStorage assocElement = _elements[_faces[id].element1_id];
    	FemaxxEntityStorage thisFace;
    	thisFace.geometryType.makeTriangle();
    	thisFace.global_id = _faces[id].global_id;
    	thisFace.interpOrder = assocElement.interpOrder;
    	thisFace.physicalTag = _faces[id].physicalTag;

    	// Get the internal element vertex indices associated with this face as a subentity
    	std::vector<int> subentityVertexIndices = Dune::CurvilinearElementInterpolator<ct, 2, 3>::subentityInternalCoordinates(assocElement.geometryType, thisFace.interpOrder, 1, _faces[id].element1_subentityNo);

    	// Calculate the localId's of vertices of this face by extracting them from the element vertex Ids
    	for(int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIds.push_back(assocElement.vertexIds[subentityVertexIndices[i]]); }

    	return thisFace;
    }


    // Returns the storage data related to this edge, except of explicit vertex coordinates
    FemaxxEntityStorage getElementData(id_t id) { return _elements[id]; }


    // Construct and retrieve geometry class associated to this edge
    EdgeGeometry getEdgeGeometry(id_t id)
    {
    	FemaxxEntityStorage edgeData = getEdgeData(id);
    	return getEntityGeometry<1>(edgeData.geometryType, edgeData.vertexIds, edgeData.interpOrder);
    }

    // Construct and retrieve geometry class associated to this face
    FaceGeometry getFaceGeometry(id_t id)
    {
    	FemaxxEntityStorage faceData = getFaceData(id);
    	return getEntityGeometry<2>(faceData.geometryType, faceData.vertexIds, faceData.interpOrder);
    }

    // Construct and retrieve geometry class associated to this element
    ElementGeometry getElementGeometry(id_t id)
    {
    	FemaxxEntityStorage elementData = getElementData(id);
    	return getEntityGeometry<3>(elementData.geometryType, elementData.vertexIds, elementData.interpOrder);
    }

    /** Return pointer to Octree or 0 if it is not constructed. */
    LooseOctree<Tet>* get_octree() const { return _octree; }


    /** Compute center and extent (halved) of the bounding box of the whole mesh. */
    void get_bounding_box(Vertex & center, Vertex & extent) const
    {
    	Vertex min = _points[0];
    	Vertex max = min;

        for (id_t i = 1; i < _points.size(); i ++) {
            min[0] = std::min(min[0], _points[i][0]);
            min[1] = std::min(min[1], _points[i][1]);
            min[2] = std::min(min[2], _points[i][2]);

            max[0] = std::max(max[0], _points[i][0]);
            max[1] = std::max(max[1], _points[i][1]);
            max[2] = std::max(max[2], _points[i][2]);
        }
        center = min + max;  center *= 0.5;
        extent = max - min;  extent *= 0.5;
    }


    /** Return edge connecting node0 and node1. A runtime_error is
     *  thrown if the edge does not exist.
     */

    // FIXME: The maps will be deleted to save space. If the lookup functionality is not necessary it should just be deleted. Check which routines still use it
    id_t lookupEdge(id_t node0, id_t node1) {
    	FemaxxEdgeMapKey thisEdgeKey;
    	thisEdgeKey.node0 = node0;
    	thisEdgeKey.node1 = node1;

        if (thisEdgeKey.node0 > thisEdgeKey.node1)  { std::swap(thisEdgeKey.node0, thisEdgeKey.node1); }

        EdgeKey2EdgeIdMap::iterator iter = _edgemap.find(thisEdgeKey);
        if (iter == _edgemap.end())  { throw std::runtime_error("Accessed non-existing edge."); }

        return (*iter).second;
    }

    /** Return face connecting by node0, node1 and node2. A runtime_error is
     *  thrown if the face does not exist.
     */
    id_t lookupFace(id_t node0, id_t node1, id_t node2){
    	FemaxxFaceMapKey thisFaceKey;
    	thisFaceKey.node0 = node0;
    	thisFaceKey.node1 = node1;
    	thisFaceKey.node2 = node2;

        // sort node ids in ascending order
    	if (thisFaceKey.node0 > thisFaceKey.node1)  { std::swap(thisFaceKey.node0, thisFaceKey.node1); }
    	if (thisFaceKey.node1 > thisFaceKey.node2)  { std::swap(thisFaceKey.node1, thisFaceKey.node2); }
    	if (thisFaceKey.node0 > thisFaceKey.node1)  { std::swap(thisFaceKey.node0, thisFaceKey.node1); }

    	// Check if the face is present in any of the 3 face maps, otherwise throw error
    	FaceKey2FaceIdMap::iterator iterInt = _internalInternalMap.find(thisFaceKey);

        if (iterInt != _internalInternalMap.end()) { return (*iterInt).second; }
        else
        {
        	FaceKey2FaceIdMap::iterator iterProc = _processInternalMap.find(thisFaceKey);
        	if (iterProc != _processInternalMap.end()) { return (*iterProc).second; }
        	else
        	{
        		FaceKey2FaceIdMap::iterator iterBnd = _boundaryInternalMap.find(thisFaceKey);
        		if (iterBnd != _processInternalMap.end()) { return (*iterBnd).second; }
        		else { throw std::runtime_error("Accessed non-existing face."); }
        	}
        }
    }


    /** Return true if p is in inside the mesh, i.e. inside at least one tet of the mesh. */
    bool is_inside(const Vertex & p) const {
        std::vector<id_t> containerIds;
        find_tets_by_point(p, containerIds);
        return containerIds.size() > 0;
    }

    // Iterators over local indices of the mesh
    // NOTE: There are no iterators over entities because there is no entity object in the core mesh
    // There will be generic entity in the wrapper because the wrapper will define an entity object
    Index2IndexMap::iterator elementIndexBegin()  { return _elementGlobal2LocalMap.begin(); }
    Index2IndexMap::iterator faceIndexBegin()     { return _faceGlobal2LocalMap.begin(); }
    Index2IndexMap::iterator edgeIndexBegin()     { return _edgeGlobal2LocalMap.begin(); }
    Index2IndexMap::iterator vertexIndexBegin()   { return _vertexGlobal2LocalMap.begin(); }

    Index2IndexMap::iterator elementIndexEnd()  { return _elementGlobal2LocalMap.end(); }
    Index2IndexMap::iterator faceIndexEnd()     { return _faceGlobal2LocalMap.end(); }
    Index2IndexMap::iterator edgeIndexEnd()     { return _edgeGlobal2LocalMap.end(); }
    Index2IndexMap::iterator vertexIndexEnd()   { return _vertexGlobal2LocalMap.end(); }

    // This construction allows fast iteration over faces of different type
    std::vector<int>::iterator faceInternalBegin()         { return _internalFaceIndex.begin(); }
    std::vector<int>::iterator faceDomainBoundaryBegin()   { return _domainBoundaryFaceIndex.begin(); }
    std::vector<int>::iterator faceProcessBoundaryBegin()  { return _processBoundaryFaceIndex.begin(); }

    std::vector<int>::iterator faceInternalEnd()         { return _internalFaceIndex.end(); }
    std::vector<int>::iterator faceDomainBoundaryEnd()   { return _domainBoundaryFaceIndex.end(); }
    std::vector<int>::iterator faceProcessBoundaryEnd()  { return _processBoundaryFaceIndex.end(); }



    /** Return pointer to single TetMesh instance. */
    static TetMesh* get_instance() { return instance_; }


protected:


    /* ***************************************************************************
     * Section: Auxiliary Methods
     * ***************************************************************************/

    // Writes debug info to the command line
    // TODO: Use IFDEF to manipulate between no output, all output, or only master process output
    void print_debug(std::string s)
    {
        if (verbose) { std::cout << "Process_" << _mpihelper.rank() << ": " << s << std::endl; }
    }

    // Gets all sets local indices of a tetrahedron which correspond to edges
    std::vector<std::vector<int> > internalEdgesOfTriangle()
	{
    	std::vector<std::vector<int> > edges (3, std::vector<int>());

    	edges[0] = {0, 1};
    	edges[1] = {0, 2};
    	edges[2] = {1, 2};

    	return edges;
	}

    // Gets all sets local indices of a tetrahedron which correspond to edges
    std::vector<std::vector<int> > internalEdgesOfTetrahedron()
	{
    	std::vector<std::vector<int> > edges (6, std::vector<int>());

    	edges[0] = {0, 1};
    	edges[1] = {0, 2};
    	edges[2] = {0, 3};
    	edges[3] = {1, 2};
    	edges[4] = {1, 3};
    	edges[5] = {2, 3};

    	return edges;
	}

    // Gets all sets local indices of a tetrahedron which correspond to faces
    std::vector<std::vector<int> > internalFacesOfTetrahedron()
	{
    	std::vector<std::vector<int> > faces (4, std::vector<int>());

    	faces[0] = {0, 1, 2};
    	faces[1] = {0, 1, 3};
    	faces[2] = {0, 2, 3};
    	faces[3] = {1, 2, 3};

    	return faces;
	}

    // Get curved geometry of an entity
    // TODO: assert mydim == element geometry type dim
    template<int mydim>
    Dune::LagrangeGeometry<ct, mydim, 3> getEntityGeometry(
    		Dune::GeometryType gt,
    		std::vector<int> & vertexIds,
    		int order) const
    {
    	VertexVector entityVertices;
    	for (int i = 0; i < vertexIds.size(); i++) { entityVertices.push_back(_points[vertexIds[i]]); }

    	return Dune::LagrangeGeometry<ct, mydim, 3> (gt, entityVertices, order);
    }

    // Returns corner id's of this entity
    std::vector<id_t> getEntityCorners(
    		Dune::GeometryType gt,
    		std::vector<int> & vertexIds,
    		int order)
	{
    	std::vector<id_t> corners;

    	// Get corner number
    	int cornerNo = gt.dim() + 1;
    	//cornerNumber = ReferenceElements::general(geomType).size( geomType.dim() );

    	// Get corners
    	for (int j = 0; j < cornerNo; j++) {
    		int internalId = Dune::CurvilinearElementInterpolator::cornerID(gt, order, j );
    		corners.push_back(vertexIds[internalId]);
    	}

    	return corners;
	}



    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    // Generates all edges
    // FIXME: Original code seems to give edges some orientation
    // FIXME: Currently subentity orientation does not match the one of Dune
    void generateEdges()
    {
    	// Loop over all elements and their edges
    	for (int i = 0; i < _elements.size(); i++)
    	{
    		std::vector<id_t> elementCorners = getEntityCorners(_elements[i].geometryType, _elements[i].vertexIds, _elements[i].interpOrder);
    		std::vector<std::vector<int> > elementInternalIndicesEdges = internalEdgesOfTetrahedron();

    		for (int j = 0; j < elementInternalIndicesEdges.size(); j++)
    		{
    			// Define (key = sorted localIds of corners)
    			FemaxxEdgeMapKey thisKey;
    			thisKey.node0 = elementCorners[elementInternalIndicesEdges[j][0]];
    			thisKey.node1 = elementCorners[elementInternalIndicesEdges[j][1]];

    			// Sort in ascending order
    			thisKey.sort();

    			// If this edge has not been added already, add it
    			if (_edgemap.find(thisKey) == _edgemap.end())
    			{
    				// Store map (key -> edgeId)
    				_edgemap[thisKey] = _edges.size();

    				// Store vector (edgeId -> elemId + edgeElemIndex)
    				// Note: Edges do not have physical tag at all so we do not even store it
    				FemaxxEdgeStorage thisEdge;
    				thisEdge.global_id = 0;		// GlobalId for edge determined later using global communication
    				thisEdge.element_id = i;
    				thisEdge.subentityNo = j;
    				_edges.push_back(thisEdge);
    			}


    		}

    	}
    }

    // Generates Internal and ProcessBoundary Faces. (!!!) Assumes that all Domain Boundary Faces have been added.
    // FIXME: Original code seems to give edges some orientation
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
    void generateFaces()
    {
    	typedef std::map<FemaxxFaceMapKey, std::vector<int>> tmpFace2InfoMap;
    	tmpFace2InfoMap tmpFaceMap;

    	// Loop over all elements and their faces
    	for (int i = 0; i < _elements.size(); i++)
    	{
    		std::vector<id_t> elementCorners = getEntityCorners(_elements[i].geometryType, _elements[i].vertexIds, _elements[i].interpOrder);
    		std::vector<std::vector<int> > elementInternalIndicesFaces = internalFacesOfTetrahedron();

    		for (int j = 0; j < elementInternalIndicesFaces.size(); j++)
    		{
    			// Define (key = sorted localIds of corners)
    			FemaxxFaceMapKey thisKey;
    			thisKey.node0 = _points[elementCorners[elementInternalIndicesFaces[j][0]]].global_id;
    			thisKey.node1 = _points[elementCorners[elementInternalIndicesFaces[j][1]]].global_id;
    			thisKey.node2 = _points[elementCorners[elementInternalIndicesFaces[j][2]]].global_id;

    			// Sort in ascending order
    			thisKey.sort();

    			// Only do sth if this is not the Domain Boundary
    			if (_boundaryInternalMap.find(thisKey) == _boundaryInternalMap.end())
    			{
    				std::vector<int> connectedFaceInfo;
    				tmpFace2InfoMap::iterator iter = tmpFaceMap.find(thisKey);

    				if (iter == tmpFaceMap.end()) { connectedFaceInfo = (*iter).second; }

    				connectedFaceInfo.push_back(i);
    				connectedFaceInfo.push_back(j);

    				tmpFaceMap[thisKey] = connectedFaceInfo;
    			}
    		}

    		for (tmpFace2InfoMap::iterator iter = tmpFaceMap.begin();  iter != tmpFaceMap.end();  iter++)
    		{
    			FemaxxFaceStorage thisFace;
    			std::vector<int> connectedFaceInfo = (*iter).second;

				thisFace.global_id = 0;       // GlobalId is defined at a later stage
				thisFace.element1_id = connectedFaceInfo[0];
				thisFace.element1_subentityNo = connectedFaceInfo[1];
				thisFace.physicalTag = -1;    // At the moment physicalTag of an internal face is not defined as it could be inbetween two different elements

		    	// Store internal, domain and process boundaries separately for faster iterators
    			if (connectedFaceInfo.size() == 2)
    			{
    				// Store Map (key -> processBoundaryFaceIndex), Vector (processBoundaryFaceIndex -> faceIndex)
    				_processInternalMap[(*iter).first] = _processBoundaryFaceIndex.size();
    				_processBoundaryFaceIndex.push_back(_faces.size());

    				thisFace.element2_id = 0;	// Eventually this will be the Ghost Element Index
    			}
    			else
    			{
    				// Store Map (key -> internalFaceIndex), Vector (internalFaceIndex -> faceIndex)
    				_internalInternalMap[(*iter).first] = _internalFaceIndex.size();
    				_internalFaceIndex.push_back(_faces.size());

    				thisFace.element2_id = connectedFaceInfo[2];  // This is the 2nd neighbor of this internal face
    			}

    			_faces.push_back(thisFace);
    		}
    	}
    }


    /** Construct octree for locating tetrahedrons in mesh */
    void construct_octree() {
        rAssert(_octree == 0);

        // bounding box of whole mesh
        Vector3 center, extent;
        get_bounding_box(center, extent);

        // octree length is the largest component of extent
        double length = extent.x;
        if (extent.y > length)
            length = extent.y;
        if (extent.z > length)
            length = extent.z;

        // construct LooseOctree with large max depth
        _octree = new LooseOctree<Tet>(center, length, 100);

        // loop over all tets and insert them in the octree
        for (id_t t = 0; t < get_nof_tets(); t ++) {
            _octree->add_node(get_tet(t));
        }
        int max_depth, nof_octants, nof_nodes;
        double avg_node_depth;
        _octree->get_statistics(max_depth, avg_node_depth, nof_octants, nof_nodes);
        rInfo("Octree stats: max. depth=%d, #octants=%d, #nodes=%d, avg. #node per octant=%.2f, avg. node depth=%.2f",
              max_depth, nof_octants, nof_nodes, static_cast<double>(nof_nodes)/nof_octants, avg_node_depth);

    #if 0
        // Vector3 p (0.000001, 0.000001, 0.000001);
        // Vector3 p (0.000000, 0.000000, 0.000000);
        // Vector3 p (1.00000071059, 0.500001458653, 0.106184447056);
        Vector3 p (0.1, 0.0, 2.6);
        std::vector<OctreeNode *> nodes;
        int nof_visited;

        cout << "Finding tets containing " << p << "\n\n";

        tet_iterator tet_it;
        cout << "Exhausive search yields (is_inside_bounding_box_gracious)\n";
        for (tet_it = tet_begin(); tet_it != tet_end(); ++ tet_it) {
            Tet* tet = *tet_it;
            if (tet->is_inside_bounding_box_gracious(p)) {
                std::cout << "inside tet " << tet->get_id() << "\n";
                for (int i = 0; i < 4; i ++) {
                    cout << "   " << tet->get_corner(i)->get_coord() << "\n";
                }
            }
        }

        cout << "Exhausive search yields (is_inside_bounding_box)\n";
        for (tet_it = tet_begin(); tet_it != tet_end(); ++ tet_it) {
            Tet* tet = *tet_it;
            if (tet->is_inside_bounding_box(p)) {
                std::cout << "inside tet " << tet->get_id() << "\n";
                for (int i = 0; i < 4; i ++) {
                    cout << "   " << tet->get_corner(i)->get_coord() << "\n";
                }
            }
        }

        cout << "Exhausive search yields (is_inside)\n";
        for (tet_it = tet_begin(); tet_it != tet_end(); ++ tet_it) {
            Tet* tet = *tet_it;
            if (tet->is_inside(p)) {
                std::cout << "inside tet " << tet->get_id() << "\n";
                for (int i = 0; i < 4; i ++) {
                    cout << "   " << tet->get_corner(i)->get_coord() << "\n";
                }
            }
        }

        cout << "\nOctree: is_inside_bounding_box\n";
        nof_visited = 0; nodes.resize(0);
        _octree->find_nodes_by_point(p, nodes, nof_visited,
                                     (&OctreeNode::is_inside_bounding_box));
        std::cout << "nof. nodes visited=" << nof_visited << "\n";
        std::cout << "nof. nodes found=" << nodes.size() << "\n";
        std::vector<OctreeNode *>::iterator it;
        for (it = nodes.begin(); it != nodes.end(); it ++) {
            Tet *tet = dynamic_cast<Tet*>(*it);
            if (tet->is_inside(p)) {
                std::cout << "inside tet " << tet->get_id() << "\n";
                for (int i = 0; i < 4; i ++) {
                    cout << "   " << tet->get_corner(i)->get_coord() << "\n";
                }
            }
        }

        cout << "\nOctree: is_inside_bounding_box_gracious\n";
        nof_visited = 0; nodes.resize(0);
        _octree->find_nodes_by_point(p, nodes, nof_visited,
                                     (&OctreeNode::is_inside_bounding_box_gracious));
        std::cout << "nof. nodes visited=" << nof_visited << "\n";
        std::cout << "nof. nodes found=" << nodes.size() << "\n";
        for (it = nodes.begin(); it != nodes.end(); it ++) {
            Tet *tet = dynamic_cast<Tet*>(*it);
            if (tet->is_inside(p)) {
                std::cout << "inside tet " << tet->get_id() << "\n";
                for (int i = 0; i < 4; i ++) {
                    cout << "   " << tet->get_corner(i)->get_coord() << "\n";
                }
            }
        }

        cout << "\nOctree: Tet::is_inside\n";
        nof_visited = 0; nodes.resize(0);
        _octree->find_nodes_by_point(p, nodes, nof_visited,
                                     (LooseOctree::filter)(&Tet::is_inside));
        std::cout << "nof. nodes visited=" << nof_visited << "\n";
        std::cout << "nof. nodes found=" << nodes.size() << "\n";
        for (it = nodes.begin(); it != nodes.end(); it ++) {
            Tet *tet = dynamic_cast<Tet*>(*it);
            if (tet->is_inside(p)) {
                std::cout << "inside tet " << tet->get_id() << "\n";
                for (int i = 0; i < 4; i ++) {
                    cout << "   " << tet->get_corner(i)->get_coord() << "\n";
                }
            }
        }

        cout << "\nOctree, find_one: Tet::is_inside\n";
        nof_visited = 0;
        // this is a dangerous typecast: it's only valid since we know all nodes in the
        // octree are Tets.
        OctreeNode* node = _octree->find_one_node_by_point(p, nof_visited, (LooseOctree::filter)(&Tet::is_inside));
        std::cout << "nof. nodes visited=" << nof_visited << "\n";
        if (node) {
            std::cout << "nof. nodes found=" << 1 << "\n";
            Vector4 p_simplex;
            dynamic_cast<Tet *>(node)->cartesian_to_simplex(p, p_simplex);
            cout << "is_inside(p) is " << dynamic_cast<Tet *>(node)->is_inside(p) << "\n";
            std::cout << "inside tet " << dynamic_cast<Tet *>(node)->get_id() << "\n";
            for (int i = 0; i < 4; i ++) {
                cout << "   " << dynamic_cast<Tet*>(node)->get_corner(i)->get_coord() << "\n";
            }
            cout << "p_simplex=" << p_simplex << "\n";
        }

        {
            Vector3 center, extent;
            get_tet(6488)->get_bounding_box(center, extent);
            cout.precision(17);
            cout << "\n\ncenter: " << center << "\nextent: " << extent << "\n";
            center -= p;
            center.x = fabs(center.x);
            center.y = fabs(center.y);
            center.z = fabs(center.z);
            cout << "|c-p|: " << center
                 << "\nextend - |c-p|:" << extent - center
                 << "\nextend_ - |c-p|:" << extent*(1+1e-13) - center
                 << "\nbool: " << (center <= extent) << "\n";
            _octree->add_node(get_tet(6488))->dump(cout);
            cout << "eps_mach=" << std::numeric_limits<double>::epsilon() << "\n";
        }
    #endif
    }


	// [TODO] Possibly inefficient implementation:
	// 1) Excessive communication. Package sent to all can only be used by 1 other process
	// 2) Imbalance bottleneck. All processes have to wait until the ones with biggest number of procBoundaries communicate
    void generateProcessBoundaryNeighbors()
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = _mpihelper.getCollectiveCommunication();

    	// 0) Make space in the neighbor process rank vector
    	_processBoundaryNeighborProcess.reserve(_processInternalMap.size());

    	// 1) Each process must somehow learn which other process it shares each process boundary with
    	// 1.1) Every process tells everyone how many interprocessor boundaries he has. Calculate max = max(procbnd per process), loops 1 to max
    	int thisProcBoundarySize = _processInternalMap.size();
    	int maxProcBoundarySize = collective_comm.max(thisProcBoundarySize);


    	int thisProcFace = 0;
    	FaceKey2FaceIdMap::iterator iter = _processInternalMap.begin();

    	// 1.2) Per each loop allgather - send 1 face global key from each to all. If all faces already sent, send empty face
    	for (int thisProcFace = 0; thisProcFace < maxProcBoundarySize; thisProcFace++)
    	{
    		FemaxxFaceMapKey thisFaceKey;
    		if (thisProcFace < thisProcBoundarySize)  { thisFaceKey = (*iter).first; iter++; }
    		else {
    			// All process boundaries have been sent by this process. Send fake faces now.
    			// Do NOT increase iterator
    			thisFaceKey.node0 = -1;
    			thisFaceKey.node1 = -1;
    			thisFaceKey.node2 = -1;
    		}

    		std::vector<FemaxxFaceMapKey> keyset (_mpihelper.size());

    		_mpihelper.allgather(&thisFaceKey, 1, reinterpret_cast<FemaxxFaceMapKey*> (keyset.data()) );

    		// 1.2.2) Loop over all received faces. If a face is present, note its sender rank
    		for (int iProc = 0; iProc < keyset.size(); iProc++)
    		{
    			// Skip the vertex if it is fake, or if it is the one sent by this same process
    			if ((iProc != _mpihelper.rank()) && (keyset[iProc][0] >= 0))
    			{
    				FaceKey2FaceIdMap::iterator iter2 = _processInternalMap.find(keyset[iProc]);

    				// If this face is present, note its sender process
    				if (iter2 != _processInternalMap.end()) { _processBoundaryNeighborProcess[(*iter2).second] = iProc; }
    			}
    		}
    	}
    }


    // NOTE: Call after generating Ghost Elements
    void generateGlobalIndices()
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = _mpihelper.getCollectiveCommunication();


    	// 1) Get edges and faces on this process that are not owned by this process
    	// *************************************************************************
    	EdgeKey2EdgeIdMap edgesNonOwned;  globalComputeNonOwnedEdges(edgesNonOwned);
    	FaceKey2FaceIdMap facesNonOwned;  globalComputeNonOwnedFaces(facesNonOwned);

    	int edgesOwned = _edges.size() - edgesNonOwned.size();
    	int facesOwned = _faces.size() - facesNonOwned.size();


    	// 2) Communicate number of edges and faces owned by each process to all
    	// *************************************************************************
    	std::vector<int> edgesOnProcess(_mpihelper.size());   // owned edges [rank]
    	std::vector<int> facesOnProcess(_mpihelper.size());   // owned faces [rank]

    	collective_comm.allgather (&edgesOwned, 1, reinterpret_cast<int*> (edgesOnProcess.data()));
    	collective_comm.allgather (&facesOwned, 1, reinterpret_cast<int*> (facesOnProcess.data()));

    	int edgesBeforeMe = 0;   // Sum(edgesOwned : rank < thisRank)
    	int facesBeforeMe = 0;   // Sum(facesOwned : rank < thisRank)

    	for (int i = 0; i < _mpihelper.size(); i++)
    	{
    	    _totalEdgesPerMesh += edgesOnProcess[i];
    	    _totalFacesPerMesh += facesOnProcess[i];

    	    if (i < _mpihelper.rank())
    	    {
    	    	edgesBeforeMe += edgesOnProcess[i];
    	    	facesBeforeMe += facesOnProcess[i];
    	    }
    	}

    	// 3) Enumerate all edges and faces that you own
    	// *************************************************************************

    	int iEdgeGlobalId = edgesBeforeMe;
    	int iFaceGlobalId = facesBeforeMe;

    	// Faces that are not shared with other processes are automatically owned by this process
    	for (int i = 0; i < _internalFaceIndex.size(); i++)        { _faces[_internalFaceIndex[i]].global_id = iFaceGlobalId++; }
    	for (int i = 0; i < _domainBoundaryFaceIndex.size(); i++)  { _faces[_domainBoundaryFaceIndex[i]].global_id = iFaceGlobalId++; }

        for (FaceKey2FaceIdMap::iterator iter = _processInternalMap.begin(); iter != _processInternalMap.end(); iter++)
        {
        	// This process owns this face if it is in the processInternalMap and if it is not in the non-owned face map
        	if (facesNonOwned.find((*iter).first) == facesNonOwned.end())  { _faces[(*iter).second].global_id = iFaceGlobalId++; }
        }

    	for (EdgeKey2EdgeIdMap::iterator iter = _edgemap.begin(); iter != _edgemap.end(); iter++)
    	{
    		// This process owns this edge if it is in the edgemap and if it is not in the non-owned edge map
    		if (edgesNonOwned.find((*iter).first) == edgesNonOwned.end())  { _edges[(*iter).second].global_id = iEdgeGlobalId++; }
    	}


    	// 4) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
    	// *************************************************************************

    	globalDistributeMissingFaceGlobalId();
    	globalDistributeMissingEdgeGlobalId(edgesNonOwned);


    	// 1) Faces
    	// 1.1) allgather - Communicate to all how many faces you own
    	// Note: You own a processFace only if your rank is smaller than that of its other neighbor
    	// 1.2) Sum up all face numbers owned by processes rank less than you, enumerate faces you own starting from this position
    	// 1.3) MPI_alltoallv - send the numbers you enumerated for each enumerated process face to its other neighbor
    	// 1.4) Use received numbers to enumerate process faces you don't own

    	// 2) Edges
    	// 2.1) Go over all edges of all process boundaries. If neighbor of that face has higher rank than you, mark edge as not owned
    	// Note: If two process faces that share this edge report different processes as neighbors, mark this edge as complicated edge
    	// 2.2) allgather - send to all how many complicated edges you own
    	// 2.3) allgather - send to all list of keys of all complicated edges.
    	// 2.4) for each complicated edge, figure out if you own it or not by checking if your rank is above everyone else who owns it
    	// 2.5) allgather - Communicate to all how many edges you own
    	// 2.6) Sum up all edge numbers owned by processes rank less than you, enumerate edges you own starting from this position
    	// 2.7) MPI_alltoallv - send the numbers you enumerated for each enumerated process edge to its other neighbor
    	// 2.8) allgather - send to all the globalIds of complicated edges you enumerated
    	// 2.9) Use this data to enumerate edges you don't own
    }


    // 2) Communicate Ghost elements
    // Requires processBoundaries to know neighboring process rank
    // Requires existence of face global index
    // TODO: Only supports tetrahedral ghost elements at the moment
    // [FIXME] Check if everything makes sense in terms of self-to-self communication
    // [FIXME] Check if there is balance between send and receive
    void generateGhostElements()
    {

    	// 1) Compute how many neighbors have with each process and the exact interpolation orders of Ghost Elements
    	// *************************************************************************************
    	std::vector< std::vector<int> > thisProcessGhostElementIndices (_mpihelper.size(), std::vector<int>() );
    	std::vector< std::vector<int> > thisProcessGhostFaceGlobalIndices (_mpihelper.size(), std::vector<int>() );

    	for (int iPBFace = 0; iPBFace < _processBoundaryFaceIndex.size(); iPBFace++ )  {
    		int thisFaceIndex = _processBoundaryFaceIndex[iPBFace];

    		int thisGhostIndex = _faces[thisFaceIndex].element1_id;
    		int thisFaceGlobalIndex = _faces[thisFaceIndex].global_id;
    		int thisNeighborRank = _processBoundaryNeighborProcess[iPBFace];
    		thisProcessGhostElementIndices[thisNeighborRank].push_back(thisGhostIndex);
    		thisProcessGhostFaceGlobalIndices[thisNeighborRank].push_back(thisFaceGlobalIndex);
    	}


    	// 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
    	// *************************************************************************************
    	std::vector< std::vector<int> > neighborProcessGhostOrder (_mpihelper.size(), std::vector<int>() );
    	ghostDistributeInterpolationOrders(thisProcessGhostElementIndices, neighborProcessGhostOrder);


    	// 3) MPI_alltoallv - Package element globalIndex + elementPhysicalTag + all interpVertex globalIds
    	// *************************************************************************************
    	std::vector<int> packageGhostElementData;
        ghostDistributeGhostElements(thisProcessGhostElementIndices, thisProcessGhostFaceGlobalIndices, neighborProcessGhostOrder, packageGhostElementData );

    	thisProcessGhostElementIndices.clear();
    	thisProcessGhostFaceGlobalIndices.clear();


    	// 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
    	// *************************************************************************************
        std::vector<std::set<int> > missingVertices (_mpihelper.size(), std::vector<int>());
        ghostInsertGhostElements(neighborProcessGhostOrder, packageGhostElementData, missingVertices);

        neighborProcessGhostOrder.clear();
        packageGhostElementData.clear();


    	// 5) Communicates to each process the number of missing vertices out of the ones it had communicated
        // Then communicate the globalId's of all missing vertices
    	// *************************************************************************************
        std::vector<int> packageMissingVertexGlobalIndices;
        std::vector<int> verticesRequested;
        std::vector<int> verticesToSend;

        ghostDistributeMissingVertexGlobalIndices(missingVertices, packageMissingVertexGlobalIndices, verticesRequested, verticesToSend);


    	// 6) Distrubute vertex coordinates and add received coordinates to the mesh
    	// *************************************************************************************
        ghostDistributeMissingVertexCoordinates (missingVertices, packageMissingVertexGlobalIndices, verticesRequested, verticesToSend);
    }


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGlobalIndices()
     * ***************************************************************************/

    // Compute all edges that you do not own
    // Return pairs (edge local index -> owning process)
    void globalComputeNonOwnedEdges(EdgeKey2EdgeIdMap & edgesNonOwned)
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = _mpihelper.getCollectiveCommunication();

    	EdgeKey2EdgeIdMap             processBoundaryEdges;				// Edges that are shared by exactly 2 processes
    	EdgeKey2EdgeIdMap			  processBoundaryComplicatedEdges;	// Edges that are shared by more than 2 processes

    	// Mark each process boundary edge as either simple or complicated
    	// For simple edges immediately mark them as owned or non-owned
    	// *************************************************************************
    	for(FaceKey2FaceIdMap::iterator iter = _processInternalMap.begin(); iter != _processInternalMap.end(); iter++)
    	{
    		FemaxxFaceMapKey thisFace = (*iter).first;
    		int neighborRank = _processBoundaryNeighborProcess[(*iter).second];

    		FemaxxEdgeMapKey thisEdge[3];
    		thisEdge[0].node0 = thisFace.node0;  thisEdge[0].node1 = thisFace.node1;
    		thisEdge[1].node0 = thisFace.node0;  thisEdge[1].node1 = thisFace.node2;
    		thisEdge[2].node0 = thisFace.node1;  thisEdge[2].node1 = thisFace.node2;

    		for (int i = 0; i < 3; i++)
    		{
    			EdgeKey2EdgeIdMap::iterator tmpIter = processBoundaryEdges.find(thisEdge[i]);

    			// If this edge has not been visited, note its neighbor rank
    			// Otherwise, check if ranks match, if not, mark it as complicated
    			if (tmpIter == processBoundaryEdges.end()) { processBoundaryEdges[thisEdge[i]] = neighborRank; }
    			else
    			{
    				if ((*tmpIter).second != neighborRank)  {
    					processBoundaryEdges.erase(tmpIter);
    					processBoundaryComplicatedEdges[thisEdge[i]] = 1;
    				}
    				else if (neighborRank < _mpihelper.rank())  { edgesNonOwned[thisEdge[i]] = neighborRank; }
    			}
    		}
    	}


    	// Communicate number of own complicated edges to all processes
    	// *************************************************************************
    	int totalRecvKeys = 0;
    	std::vector<int> processComplicatedEdgeNumbers (_mpihelper.size());
    	std::vector<int> processComplicatedEdgeDispls (_mpihelper.size());
    	collective_comm.allgather(& processBoundaryComplicatedEdges.size(), 1, reinterpret_cast<int*> (processComplicatedEdgeNumbers.data()));

    	// Communicate complicate edge keys to all
    	// *************************************************************************
    	for (int i = 0; i < _mpihelper.size(); i++)
    	{
    		totalRecvKeys += processComplicatedEdgeNumbers[i];
    		processComplicatedEdgeDispls[i] = (i == 0) ? 0 : processComplicatedEdgeDispls[i - 1] + processComplicatedEdgeNumbers[i - 1];
    	}

    	std::vector<FemaxxEdgeMapKey> packageProcessComplicatedEdgeKeysSend(processBoundaryComplicatedEdges.begin(), processBoundaryComplicatedEdges.end());
    	std::vector<FemaxxEdgeMapKey> packageProcessComplicatedEdgeKeysRecv(totalRecvKeys);
    	collective_comm.allgatherv(
    			packageProcessComplicatedEdgeKeysSend.data(),
    			processBoundaryComplicatedEdges.size(),
    			reinterpret_cast<FemaxxEdgeMapKey*>(packageProcessComplicatedEdgeKeysRecv.data()),
    			processComplicatedEdgeNumbers.data(), processComplicatedEdgeDispls.data()
    			);


    	// Mark complicated edges as non-owned if they are on this process and also on a process with lower rank
    	// We only care about processes with lower rank when checking who owns the edge
    	// *************************************************************************
    	int iData = 0;
    	for (int i = 0; i < _mpihelper.rank(); i++)
    	{
    		for (int j = 0; j < processComplicatedEdgeNumbers[i]; j++)
    		{
    			FemaxxEdgeMapKey tmpEdgeKey = packageProcessComplicatedEdgeKeysRecv[iData++];
    			EdgeKey2EdgeIdMap::iterator tmpIter = processBoundaryComplicatedEdges.find(tmpEdgeKey);
    			if (tmpIter != processBoundaryComplicatedEdges.end())
    			{
    				edgesNonOwned[tmpEdgeKey] = i;
    				processBoundaryComplicatedEdges.erase(tmpIter);
    			}
    		}
    	}
    }


    // Compute total number of faces
    void globalComputeNonOwnedFaces(FaceKey2FaceIdMap & facesNonOwned)
    {
    	for (FaceKey2FaceIdMap::iterator iter = _processInternalMap.begin(); iter != _processInternalMap.end(); iter++)
    	{
    		int thisNeighborRank = _processBoundaryNeighborProcess[(*iter).second];
    		if (thisNeighborRank < _mpihelper.rank())  { facesNonOwned[(*iter).first] = thisNeighborRank; }
    	}
    }

    // Communicates all process boundary face global Id's to the neighbors if owned
    void globalDistributeMissingFaceGlobalId()
    {
    	typedef std::pair<FemaxxFaceMapKey, int>  FaceInfo;
    	std::vector< std::vector< FaceInfo > > facesToSend (_mpihelper.size());

    	int totalRecvSize = 0;
    	std::vector<int> sendbuf, sendcounts, sdispls;
    	std::vector<int> recvbuf, recvcounts, rdispls;

    	// 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
    	// ********************************************************************************************8
    	for (FaceKey2FaceIdMap::iterator iter = _processInternalMap.begin(); iter != _processInternalMap.end(); iter++)
    	{
    		int localIndex = (*iter).second;
    		int neighborRank = _processBoundaryNeighborProcess[localIndex];

    		// If the neighbor of this face has lower rank, then add it to recv, else to send
    		if (neighborRank < _mpihelper.rank())  { recvcounts[neighborRank]++; }
    		else
    		{
    			int thisGlobalIndex = _faces[localIndex].global_id;
    			facesToSend[neighborRank].push_back(FaceInfo ((*iter).first, thisGlobalIndex ));
    		}
    	}


    	// 2) Fill in communication arrays
    	// ********************************************************************************************8
    	for (int i = 0; i < _mpihelper.size(); i++)
    	{
    		sendcounts.push_back(facesToSend[i].size());
    		totalRecvSize += recvcounts[i];
    		sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
    		rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );

    		for (int j = 0; j < facesToSend[i]; j++)
    		{
    			sendbuf.push_back(facesToSend[i][j].second);
    			sendbuf.push_back(facesToSend[i][j].first.node0);
    			sendbuf.push_back(facesToSend[i][j].first.node1);
    			sendbuf.push_back(facesToSend[i][j].first.node2);
    		}

    	}

    	// 3) MPI_alltoall key + globalId
    	// ********************************************************************************************8
    	recvbuf.reserve(totalRecvSize);
    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


    	// 4) Mark all missing faces
    	// ********************************************************************************************8
    	// Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
    	int iData = 0;
    	for (int iProc = 0; iProc < _mpihelper.size(); iProc++)
    	{
    		for (int iFace = 0; iFace < recvcounts[iProc]; iFace++)
    		{
    			FemaxxFaceMapKey thisKey;
    			int thisGlobalId = recvbuf[iData++];
    			thisKey.node0 = recvbuf[iData++];
    			thisKey.node1 = recvbuf[iData++];
    			thisKey.node2 = recvbuf[iData++];

    			int thisLocalId = _processInternalMap[thisKey];
    			_faces[thisLocalId].global_id = thisGlobalId;
    		}
    	}
    }


    // Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
    // [FIXME] clear unnecessary arrays
    void globalDistributeMissingEdgeGlobalId(EdgeKey2EdgeIdMap & edgesNonOwned)
    {
    	// 1) MPI_alltoall - tell each process how many edges you want from it
    	// **************************************************************************

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();

    	std::vector<int> edgeNumberRequested;
    	std::vector<int> edgeNumberToSend (_mpihelper.size());
    	std::vector< std::vector< FemaxxEdgeMapKey > >  requestedKeys(_mpihelper.size());

    	for (EdgeKey2EdgeIdMap::iterator iter = edgesNonOwned.begin(); iter != edgesNonOwned.end(); iter++)  { requestedKeys[(*iter).second].push_back((*iter).first); }
    	for (int i = 0; i < _mpihelper.size(); i++) { edgeNumberRequested.push_back(requestedKeys[i].size()); }

    	MPI_Alltoall(edgeNumberRequested.data(), 1, MPI_INT, reinterpret_cast<int*>(edgeNumberToSend.data()), 1, MPI_INT, comm);

    	// 2) MPI_alltoall - send each process the edge keys you want from it
    	// **************************************************************************
    	int totalRecvData = 0;
    	std::vector<int> packageEdgesRequested, sendcounts, sdispls;
    	std::vector<int> packageEdgesToSend,    recvcounts, rdispls;


    	for (int iProc = 0; iProc < _mpihelper.size(); iProc++) {
    		totalRecvData += edgeNumberToSend[iProc];

    		for (int iKey = 0; iKey < requestedKeys[iProc].size(); iKey++)
    		{
    			packageEdgesRequested.push_back(requestedKeys[iProc][iKey].node0);
    			packageEdgesRequested.push_back(requestedKeys[iProc][iKey].node1);

    			sendcounts.push_back(2 * edgeNumberRequested[iProc]);
    			recvcounts.push_back(2 * edgeNumberToSend[iProc]);

        		sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + sendcounts[iProc-1] );
        		rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + recvcounts[iProc-1] );
    		}
    	}

    	packageEdgesToSend.reserve(2 * totalRecvData);
    	MPI_Alltoallv (packageEdgesRequested.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(packageEdgesToSend.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


    	// 3) MPI_alltoall - send each process the global indices of edges it wants from you in the received order
    	// **************************************************************************

    	std::vector<int> packageEdgeGlobalIndicesSend;
    	std::vector<int> packageEdgeGlobalIndicesRecv(totalRecvData);

    	// Package globalIndices of the edges requested by other processes
    	int iData = 0;
       	for (int iProc = 0; iProc < _mpihelper.size(); iProc++) {
    		for (int iFace = 0; iFace < edgeNumberToSend[iProc]; iFace++)
    		{
    			FemaxxEdgeMapKey thisKey;
    			thisKey.node0 = packageEdgesToSend[iData++];
    			thisKey.node1 = packageEdgesToSend[iData++];

    			packageEdgeGlobalIndicesSend.push_back(_edges[_edgemap[thisKey]].global_id);
    		}

    		sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + edgeNumberToSend[iProc-1] );
    		rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + edgeNumberRequested[iProc-1] );
       	}

       	MPI_Alltoallv (packageEdgeGlobalIndicesSend.data(), edgeNumberToSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(packageEdgeGlobalIndicesRecv.data()), edgeNumberRequested.data(), rdispls.data(), MPI_INT, comm );


    	// 4) Mark all missing edges
    	// **************************************************************************
    	iData = 0;
       	for (int iProc = 0; iProc < _mpihelper.size(); iProc++) {
    		for (int iFace = 0; iFace < requestedKeys[iProc].size(); iFace++)
    		{
    			_edges[_edgemap[requestedKeys[iProc][iFace]]].global_id = packageEdgeGlobalIndicesRecv[iData++];
    		}
       	}
    }


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGhostElements()
     * ***************************************************************************/

	// MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
    void ghostDistributeInterpolationOrders(
    		std::vector< std::vector<int> > & thisProcessGhostElementIndices,
    		std::vector< std::vector<int> > & neighborProcessGhostOrder
    )
    {
    	std::vector<int> sendbuf, sendcounts, sdispls;
    	std::vector<int> recvbuf, recvcounts, rdispls;
    	recvbuf.reserve(_processBoundaryFaceIndex.size());

    	for (int i = 0; i < thisProcessGhostElementIndices.size(); i++)
    	{
    		for (int j = 0; j < thisProcessGhostElementIndices[i].size(); j++) { sendbuf.push_back(_elements[thisProcessGhostElementIndices[i][j]].interpOrder); }
    		sendcounts.push_back(thisProcessGhostElementIndices[i].size());
    		sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );

    		// For ghost elements we send/receive same amount to/from each processor
    		recvcounts.push_back(sendcounts[i]);
    		rdispls.push_back(sdispls[i]);
    	}

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


    	int iData = 0;

    	for (int i = 0; i < thisProcessGhostElementIndices.size(); i++)
    	{
    		for (int j = 0; j < thisProcessGhostElementIndices[i].size(); j++)
    		{
    			neighborProcessGhostOrder[i].push_back(recvbuf[iData++]);
    		}
    	}
    }


	// MPI_alltoallv - Package element globalIndex + elementPhysicalTag + assocFaceGlobalIndex + all interpVertex globalIds
    void ghostDistributeGhostElements(
    		std::vector< std::vector<int> > & thisProcessGhostElementIndices,
        	std::vector< std::vector<int> > & thisProcessGhostFaceGlobalIndices,
        	std::vector< std::vector<int> > & neighborProcessGhostOrder,
    		std::vector<int> & recvbuf )
    {
    	std::vector<int> sendbuf, sendcounts, sdispls;
    	std::vector<int> recvcounts, rdispls;

    	// Calculates total amount of integers to receive during DoF communication stage
    	int totalRecvSize = 0;

    	for (int i = 0; i < thisProcessGhostElementIndices.size(); i++)
    	{
    		int thisSendCounts = 0;
    		int thisRecvCounts = 0;

    		for (int j = 0; j < thisProcessGhostElementIndices[i].size(); j++)
    		{
    			int ghostElementIndex = thisProcessGhostElementIndices[i][j];
    			int ghostElementFaceGlobalIndex = thisProcessGhostFaceGlobalIndices[i][j];
    			int thisDofNum = _elements[ghostElementIndex].vertexIds.size();
    			thisSendCounts += 3 + thisDofNum;
    			thisRecvCounts += 3 + Dune::CurvilinearElementInterpolator<ct, 3, 3>::dofPerOrder(neighborProcessGhostOrder[i][j]);
    			totalRecvSize += thisRecvCounts;

    			sendbuf.push_back(_elements[ghostElementIndex].global_id);
    			sendbuf.push_back(_elements[ghostElementIndex].physicalTag);
    			sendbuf.push_back(ghostElementFaceGlobalIndex);
    			for (int iDof = 0; iDof < thisDofNum; iDof++)
    			{
    				sendbuf.push_back(_points[_elements[thisProcessGhostElementIndices[i][j]].vertexIds[iDof]].global_id);
    			}
    		}

    		sendcounts.push_back(thisSendCounts);
    		recvcounts.push_back(thisRecvCounts);
    		sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
    		rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
    	}

    	recvbuf.reserve(totalRecvSize);

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );
    }


    // Add received elements to the mesh
    void ghostInsertGhostElements (
    		std::vector< std::vector<int> > & neighborProcessGhostOrder,
    		std::vector< int > & packageGhostElementData,
    		std::vector<std::set<int> > & missingVertices
    )
    {
    	int iData = 0;
    	Dune::GeometryType tetGeometryType;
    	tetGeometryType.makeTetrahedron();

    	for (int i = 0; i < neighborProcessGhostOrder.size(); i++) {

    		for (int j = 0; j < neighborProcessGhostOrder[i].size(); j++)
    		{
    			FemaxxEntityStorage thisElement;
    			thisElement.geometryType = tetGeometryType;
    			thisElement.global_id = packageGhostElementData[iData++];
    			thisElement.interpOrder = neighborProcessGhostOrder[i][j];
    			thisElement.physicalTag = packageGhostElementData[iData++];

    			int associatedFaceGlobalIndex = packageGhostElementData[iData++];

    			int thisElementDof = Dune::CurvilinearElementInterpolator<ct, 3, 3>::dofPerOrder(thisElement.interpOrder);
    			for (int iDof = 0; iDof < thisElementDof; iDof++)
    			{
    				int thisVertexGlobalIndex = packageGhostElementData[iData++];
    				Index2IndexMap::iterator vertexIter = _vertexGlobal2LocalMap.find(thisVertexGlobalIndex);

    				// If this vertex already exists, just reuse it. Otherwise, create new vertex, and later request its coordinate
    				if (vertexIter != _vertexGlobal2LocalMap.end()) { thisElement.vertexIds.push_back((*vertexIter).second); }
    				else
    				{
    					thisElement.vertexIds.push_back(_points.size());

    					Vertex fakeCoord;
    					insertVertex(thisVertexGlobalIndex, fakeCoord);

    					// Note that this vertex needs communicating
    					missingVertices[i].insert(thisVertexGlobalIndex);
    				}
    			}

    			// Associate a (hopefully processBoundary) face with this ghost element
    			_faces[_faceGlobal2LocalMap[associatedFaceGlobalIndex]].element2_id = _ghostElements.size();
    			_ghostElements.push_back(thisElement);
    		}
    	}
    }


    // Communicates to each process the number of missing vertices out of the ones it had communicated
    // Then communicate the globalId's of all missing vertices
    void ghostDistributeMissingVertexGlobalIndices(
    		std::vector<std::set<int> > & missingVertices,
    		std::vector<int> & recvbuf,
    		std::vector<int> & verticesRequested,
    		std::vector<int> & verticesToSend
    )
    {
    	std::vector<int> sendbuf, sendcounts, sdispls;
    	std::vector<int> recvcounts, rdispls;
    	recvbuf.reserve(_mpihelper.size());

    	// 4.1) MPI_alltoallv - tell each process the number of coordinates you want from it
    	for (int i = 0; i < missingVertices.size(); i++)
    	{
    		sendbuf.push_back(missingVertices[i].size());

    		sendcounts += 1;
    		recvcounts += 1;

    		sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
    		rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
    	}

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


    	// 4.2) MPI_alltoallv - tell each process the list of globalIds of coordinates you want from it
    	// Cleanup
    	sendbuf.clear(); sendcounts.clear(); sdispls.clear();
    	recvcounts.clear(); rdispls.clear();
    	int totalRecvSize = 0;
    	int k = 0;

    	for (int i = 0; i < missingVertices.size(); i++)
    	{
    		int thisSendSize = missingVertices[i].size();
    		int thisRecvSize = recvbuf[k++];

    		recvcounts.push_back(thisRecvSize);
    		sendcounts.push_back(thisSendSize);
    		totalRecvSize += thisRecvSize;

    		for (std::set<int>::iterator mVertIter = missingVertices[i].begin(); mVertIter != missingVertices[i].end(); mVertIter++ )  { sendbuf.push_back(*mVertIter); }

    		sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
    		rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
    	}

    	recvbuf.clear(); recvbuf.reserve(totalRecvSize);
    	MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );

    	// We will require the information about requested and sent vertices when we communicate the coordinates
    	verticesRequested.swap(sendcounts);
    	verticesToSend.swap(recvcounts);
    }


    // Distrubute vertex coordinates and add received coordinates to the mesh
    void ghostDistributeMissingVertexCoordinates (
    		std::vector<std::set<int> > & missingVertices,
    		std::vector<int> & packageMissingVertexGlobalIndices,
            std::vector<int> & verticesRequested,
            std::vector<int> & verticesToSend
    )
    {
    	// 4.3) MPI_alltoallv - package list of globalId+coordinate for each process and send it
    	std::vector<int> sendcounts, sdispls;
    	std::vector<int> recvcounts, rdispls;
    	std::vector<double> recvbuf, sendbuf;

    	int iData = 0;
    	int totalRecvSize = 0;

    	for (int i = 0; i < _mpihelper.size(); i++)
    	{
    		// Go through all vertices requested from this process. Package coordinates
    		for (int j = 0; j < recvcounts[i]; j++)
    		{
    			int thisVertexGlobalIndex = packageMissingVertexGlobalIndices[iData++];
    			int thisVertexLocalIndex = _vertexGlobal2LocalMap[thisVertexGlobalIndex];
    			Vertex p = _points[thisVertexLocalIndex].coord;

    			for (int iDim = 0; iDim < 3; iDim++)  { sendbuf.push_back(p[iDim]); }
    		}

    		// We communicate (coord = 3 doubles) for each sent/received vertex
    		// We now receive the amount we sent before, and send the amount we received before
    		int thisSendSize = 3 * verticesRequested[i];
    		int thisRecvSize = 3 * verticesToSend[i];

    		sendcounts[i] = thisSendSize;
    		recvcounts[i] = thisRecvSize;
    		totalRecvSize += thisRecvSize;

    		sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
    		rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
    	}

    	recvbuf.reserve(totalRecvSize);

    	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    	MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE, reinterpret_cast<double*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_DOUBLE, comm );



    	// Assign coordinates to all missing vertices
    	iData = 0;

    	for (int i = 0; i < missingVertices.size(); i++)
    	{
    		for (std::set<int>::iterator mVertIter = missingVertices[i].begin(); mVertIter != missingVertices[i].end(); mVertIter++ )
    		{
    			Vertex thisCoord;
    			thisCoord[0] = recvbuf[iData++];
    			thisCoord[1] = recvbuf[iData++];
    			thisCoord[2] = recvbuf[iData++];

    			int thisVertexGlobalIndex = *mVertIter;
    			_points[_vertexGlobal2LocalMap[thisVertexGlobalIndex]].coord = thisCoord;
    		}
    	}
    }




    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    // Checks if given point fits into the bounding box of a Tetrahedron
    // FIXME: This method should be superseded by COM-CURVATURE-BOUND from LagrangeGeometry
    bool is_inside_bounding_box_gracious(const Vertex & point) const {
        const double grace_tolerance = 1e-13;
        Vertex node_center, node_extent;

        get_bounding_box(node_center, node_extent);
        node_center -= point;
        node_center.x = fabs(node_center.x);
        node_center.y = fabs(node_center.y);
        node_center.z = fabs(node_center.z);
        return node_center <= node_extent * (1.0 + grace_tolerance);
    }

    // Gets a box in which this Tetrahedron fits
    // FIXME: This method should be superseded by COM-CURVATURE-BOUND from LagrangeGeometry
    void get_bounding_box(Vertex center, Vertex extent) const {
        Vertex min, max, coord;

        min = get_corner(0)->get_coord();
        max = min;
        for (int i = 1; i < 4; i ++) {
            coord = get_corner(i)->get_coord();
            min.make_floor(coord);
            max.make_ceil(coord);
        }
        center = max.mid_point(min);
        extent = 0.5 * (max - min);
    }

    // Checks if a point is inside the element
    // FIXME: This method should be superseded by isInside method from LagrangeGeometry
    bool is_inside_gracious(int id, const Vertex & point) const {
        const double grace_tolerance = 1e-14;
        Vector4 simplex_coord;

        cartesian_to_simplex(point, simplex_coord);
        for (int i = 0; i < 4; i ++)
            if (simplex_coord[i] < -grace_tolerance)
                return false;
        for (int i = 0; i < 4; i ++)
            if (simplex_coord[i] > 1.0 + grace_tolerance)
                return false;
        return true;
    }

    /** Find all tets in mesh which contain the specified point p. The found
        tets are returned in "tets". If p is outside the mesh an empty vector
        is returned. If p is on a shared face/edge/point, "tets" may contain more
        than one entry. */

    // TODO: PBE file allows femaxx to count time. Use alternative in Dune?
    void find_tets_by_point(const Vertex p, std::vector<id_t>& tets) const {
        assert(_octree);

        //pbe_start(132, "Octree traversal");

        // Get list of nodes (Tets) whose bounding box contain point p
        std::vector<id_t> nodes;
        int nof_visited = 0;
        _octree->find_nodes_by_point(p, nodes, nof_visited, &is_inside_bounding_box_gracious);

        // Create list of tets containing point p
        tets.resize(0);
        for (int i = 0; i < nodes.size(); i++) {
        	if (is_inside_gracious(nodes[i], p)) { tets.push_back(nodes[i]); }
        }

        // pbe_stop(132);
        // rDebug("find_tets_by_point: nof_found=%d, nof_visited=%d", static_cast<int>(nodes.size()), nof_visited);
    }


private: // Private members

    bool _verbose;
    bool _withGhostElements;


    // Storage necessary for user access and computation of globalIndex
    int _totalVerticesPerMesh;
    int _totalEdgesPerMesh;
    int _totalFacesPerMesh;
    int _totalElementsPerMesh;

    // Stores vertex coordinates and globalIds
    std::vector<FemaxxPointStorage> _points;

    // Stores all necessary data about edges, faces and elements
    // Minimalism - edges and faces do not store interpolatory vertex ids, but refer to an element subentity
    std::vector<FemaxxEntityStorage> _elements;
    std::vector<FemaxxEntityStorage> _ghostElements;
    std::vector<FemaxxEdgeStorage> _edges;
    std::vector<FemaxxFaceStorage> _faces;

    // Maps from global to local indices
    Index2IndexMap _vertexGlobal2LocalMap;
    Index2IndexMap _edgeGlobal2LocalMap;
    Index2IndexMap _faceGlobal2LocalMap;
    Index2IndexMap _elementGlobal2LocalMap;

    // List of localIds of all faces of different structural types. Speeds up iterators
    std::vector<int> _internalFaceIndex;         // (self -> _faces index)
    std::vector<int> _domainBoundaryFaceIndex;   // (self -> _faces index)
    std::vector<int> _processBoundaryFaceIndex;  // (self -> _faces index)

    // List of all ranks of processors neighboring processorBoundaries. Index the same as _processBoundaryFaceIndex.
    std::vector<int> _processBoundaryNeighborProcess;  // (_processBoundaryFaceIndex -> neighbor rank)

    // Temporary maps necessary to
    EdgeKey2EdgeIdMap _edgemap;				   // (global edgeKey -> _edges index)
    FaceKey2FaceIdMap _internalInternalMap;    // (global faceKey -> _internalFaceIndex)
    FaceKey2FaceIdMap _boundaryInternalMap;    // (global faceKey -> _domainBoundaryFaceIndex)
    FaceKey2FaceIdMap _processInternalMap;     // (global faceKey -> _processBoundaryFaceIndex)


    // Octree used to efficiently locate elements in which the points are located
    LooseOctree<Tet>* _octree;
    /**
     * Pointer to single TetMesh instance (Singleton)
     */
    static TetMesh* instance_ = 0;          // Pointer to TetMesh singleton.

    MPIHelper &_mpihelper;
};

} // namespace mesh

#endif
