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

#include "myrlog.h"
#undef RLOG_SECTION
#define RLOG_SECTION

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include "utility/triple.h"
#include "point.h"
#include "edge.h"
#include "face.h"
#include "tet.h"
#include <fstream>
#include "looseoctree.h"

#include "materials.h"


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
 *  - [FIXME] Need to merge with parallel tetmesh
 *  - [FIXME] Need to wrap for Dune
 *  - [FIXME] Need to match Dune's internal subentity id convention
 *  - [FIXME] Need to rewrite OCTree so it fits new convention
 *  - [FIXME] Need to introduce parallel logging through verbose as in other classes
 *  - [FIXME] Need to add normal and outerNormal
 *
 *  - [FIXME] Split face and edge container. For face, store 2 element pointers, first by definition internal. For edge, it makes sense not to store physical tag
 *
 *
 * CopyPaste log
 *  - [NOTE] ParallelGrid does some weird setting of point ownership
 *
 *  ***************************************************************************/




namespace Dune {





struct FemaxxEntityStorage
{
	Dune::GeometryType geometryType;
	id_t global_id;
	std::vector<int> vertexIds;
	int interpOrder;
	int physicalTag;
};

struct FemaxxSubentityStorage
{
	id_t global_id;
	id_t element_id;
	int element_internalIndex;
	int physicalTag;
};

struct FemaxxEdgeMapKey
{
	int node0;
	int node1;
};

struct FemaxxFaceMapKey
{
	int node0;
	int node1;
	int node2;
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

    typedef std::map<FemaxxEdgeMapKey, id_t> EdgeKey2EdgeIdMap;
    typedef std::map<FemaxxFaceMapKey, id_t> FaceKey2FaceIdMap;

    /** STL random iterator for entities in TetMesh. */
    typedef std::vector<FemaxxEntityStorage>::iterator EntityIterator;

    typedef Dune::LagrangeGeometry<ct, 1, 3> EdgeGeometry;
    typedef Dune::LagrangeGeometry<ct, 2, 3> FaceGeometry;
    typedef Dune::LagrangeGeometry<ct, 3, 3> ElementGeometry;





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
    virtual ~TetMesh() {
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
    virtual void init_points(int nofPoint)
    {
        _points.reserve(nofPoint);
        _pointGlobalIds.reserve(nofPoint);
    }

    virtual void init_elements(int nofElement)
    {
        _elements.reserve(nofElement);
    }


    /** Add a new point to the mesh */
    virtual void insertVertex(int id, Vertex p)
    {
    	// ??
    	//assert(static_cast<size_t>(id) == _points.size());
    	_points.push_back(p);
    	_pointGlobalIds.push_back(id);
    }

    /** Insert an element into the mesh */
    virtual void insertElement(Dune::GeometryType gt, id_t global_id, std::vector<int> vertexIds, int order, int physicalTag)
    {
    	// TODO: Assert geometryType = tetrahedron
    	// TODO: Assert vertex length = expected for this order and geometry type    assert(id == _tets.size());

    	FemaxxEntityStorage thisElement;

    	thisElement.geometryType = gt;
    	thisElement.global_id = global_id;
    	thisElement.interpOrder = order;
    	thisElement.physicalTag = physicalTag;
    	thisElement.vertexIds = vertexIds;

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
	// TODO: Assert geometryType = triangle
	// TODO: Assert vertex length = expected for this order and geometry type    assert(id == _tets.size());

    // FIXME: Must store faces in the map (faceId -> elemId + faceElemIndex)
    virtual void insertBoundarySegment(Dune::GeometryType gt, id_t assoc_element_id, std::vector<int> vertexIds, int order, int physicalTag)
    {
    	insertBoundaryFace(gt, assoc_element_id, vertexIds, order, true, physicalTag);
    }

    virtual void insertProcessBoundary(Dune::GeometryType gt, id_t assoc_element_id, std::vector<int> vertexIds, int order, int physicalTag)
    {
    	insertBoundaryFace(gt, assoc_element_id, vertexIds, order, false, physicalTag);
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
    			if (thisKey.node0 > thisKey.node1)  { std::swap(thisKey.node0, thisKey.node1); }

    			// If this edge has not been added already, add it
    			if (_edgemap.find(thisKey) == _edgemap.end())
    			{
    				// Store map (key -> edgeId)
    				_edgemap[thisKey] = _edges.size();

    				// Store vector (edgeId -> elemId + edgeElemIndex)
    				FemaxxSubentityStorage thisEdge;
    				thisEdge.element_id = i;
    				thisEdge.element_internalIndex = j;
    				thisEdge.physicalTag = -1;	// At the moment physicalTags have no meaning for edges but that may change
    				_edges.push_back(thisEdge);
    			}


    		}

    	}
    }

    // Generates all edges
    // FIXME: Original code seems to give edges some orientation
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
    void generateFaces()
    {
    	// Loop over all elements and their faces
    	for (int i = 0; i < _elements.size(); i++)
    	{
    		std::vector<id_t> elementCorners = getEntityCorners(_elements[i]);
    		std::vector<std::vector<int> > elementInternalIndicesFaces = internalFacesOfTetrahedron();

    		for (int j = 0; j < elementInternalIndicesFaces.size(); j++)
    		{
    			// Define (key = sorted localIds of corners)
    			FemaxxFaceMapKey thisKey;
    			thisKey.node0 = elementCorners[elementInternalIndicesFaces[j][0]];
    			thisKey.node1 = elementCorners[elementInternalIndicesFaces[j][1]];
    			thisKey.node2 = elementCorners[elementInternalIndicesFaces[j][2]];

    			// Sort in ascending order
    			if (thisKey.node0 > thisKey.node1)  { std::swap(thisKey.node0, thisKey.node1); }
    			if (thisKey.node1 > thisKey.node2)  { std::swap(thisKey.node1, thisKey.node2); }
    			if (thisKey.node0 > thisKey.node1)  { std::swap(thisKey.node0, thisKey.node1); }

    			// If this face has not been added already, add it
    			bool newFace =
    					   (_boundaryInternalMap.find(thisKey) == _boundaryInternalMap.end())
    					&& (_processInternalMap.find(thisKey) == _processInternalMap.end())
    					&& (_internalInternalMap.find(thisKey) == _internalInternalMap.end());

    			if (newFace)
    			{
    				// Store map (key -> internalFaceId)
    				_internalInternalMap[thisKey] = _faces.size();

    				// Store vector (internalFaceId -> elemId + faceElemIndex)
    				FemaxxSubentityStorage thisFace;
    				thisFace.global_id = 0;       // GlobalId is defined at a later stage
    				thisFace.element_id = i;
    				thisFace.element_internalIndex = j;
    				thisFace.physicalTag = -1;    // At the moment physicalTag of an internal face is not defined as it could be inbetween two different elements

    				_faces.push_back(thisFace);
    			}
    		}
    	}
    }

    // TODO: Call after generating Ghost Elements
    void generateGlobalIndices()
    {
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

    // TODO: Implement me!!! :)
    void generateGhostElements()
    {
    	// Stupid strategy
    	// 1) Each process must somehow learn which other process it shares each process boundary with
    	// 1.1) Every process tells everyone how many interprocessor boundaries he has
    	// 1.2) Calculate max = max(procbnd per process), loops 1 to max
    	// 1.3) Per each loop allgather - send 1 face global key from each to all. If all faces already sent, send empty face
    	// 1.4) From received keys, select ones present on this process, assign neighbor process to process boundaries

    	// 2) Communicate Ghost elements
    	// 2.1) MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
    	// 2.2) MPI_alltoallv - Package all interpVertex globalIds + elementPhysicalTag

    	// 3) Each process will be missing some coordinates, they need to be communicated.
    	// 3.1) For data received from each process, see which coordinates are missing
    	// 3.2) MPI_alltoallv - tell each process the number of coordinates you want from it
    	// 3.3) MPI_alltoallv - tell each process the list of globalIds of coordinates you want from it
    	// 3.4) MPI_alltoallv - package list of globalId+coordinate for each process and send it

    	// Problems:
    	// 1) Excessive communication. Package sent to all can only be used by 1 other process
    	// 2) Imbalance bottleneck. All processes have to wait until the ones with biggest number of procBoundaries communicate
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
    // TODO: Implement static subentityIdSet in CurvilinearGeometry class
    FemaxxEntityStorage getEdgeData(id_t id)
    {
    	FemaxxEntityStorage assocElement = _elements[_edges[id].element_id];
    	FemaxxEntityStorage thisEdge;
    	thisEdge.geometryType.makeLine();
    	thisEdge.global_id = _edges[id].global_id;
    	thisEdge.interpOrder = assocElement.interpOrder;
    	thisEdge.physicalTag = _edges[id].physicalTag;		// Note that at the moment the physical tag of the edges is meaningless

    	// Get the internal element vertex indices associated with this face as a subentity
    	std::vector<int> subentityVertexIndices = Dune::CurvilinearElementInterpolator<ct, 1, 3>::subentityIdSet(assocElement.geometryType, thisEdge.geometryType, _edges[id].element_internalIndex );

    	// Calculate the localId's of vertices of this face by extracting them from the element vertex Ids
    	for(int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIds.push_back(assocElement.vertexIds[subentityVertexIndices[i]]); }

    	return thisEdge;

    }


    // Returns the storage data related to this edge, except of explicit vertex coordinates
    // This data is stored only partially - vertexIds are extracted from associated element
    // TODO: Implement static subentityIdSet in CurvilinearGeometry class
    FemaxxEntityStorage getFaceData(id_t id)
    {
    	FemaxxEntityStorage assocElement = _elements[_faces[id].element_id];
    	FemaxxEntityStorage thisFace;
    	thisFace.geometryType.makeTriangle();
    	thisFace.global_id = _faces[id].global_id;
    	thisFace.interpOrder = assocElement.interpOrder;
    	thisFace.physicalTag = _faces[id].physicalTag;

    	// Get the internal element vertex indices associated with this face as a subentity
    	std::vector<int> subentityVertexIndices = Dune::CurvilinearElementInterpolator<ct, 2, 3>::subentityIdSet(assocElement.geometryType, thisFace.geometryType, _faces[id].element_internalIndex );

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

        for (id_t i = 1; i < get_nof_points(); i ++) {
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
    virtual id_t lookupEdge(id_t node0, id_t node1) {
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
    virtual id_t lookupFace(id_t node0, id_t node1, id_t node2){
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

    /* Deletes all temporary memory */
    virtual void finalize_mesh() {
    	typedef FaceKey2FaceIdMap::iterator IterType;

        _edgemap.clear();

        _internalFaceIds.reserve(_internalInternalMap.size());
        _domainBoundaryFaceIds.reserve(_boundaryInternalMap.size());
        _processBoundaryFaceIds.reserve(_processInternalMap.size());

        for (IterType it = _internalInternalMap.begin(); it < _internalInternalMap.end();  it++) { _internalFaceIds.push_back((*it).second); }
        for (IterType it = _boundaryInternalMap.begin(); it < _boundaryInternalMap.end();  it++) { _domainBoundaryFaceIds.push_back((*it).second); }
        for (IterType it = _processInternalMap.begin();  it < _processInternalMap.end();   it++) { _processBoundaryFaceIds.push_back((*it).second); }

        _internalInternalMap.clear();
        _boundaryInternalMap.clear();
        _processInternalMap.clear();
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

    /** Return true if p is in inside the mesh, i.e. inside at least one tet of the mesh. */
    bool is_inside(const Vertex & p) const {
        std::vector<id_t> containerIds;
        find_tets_by_point(p, containerIds);
        return containerIds.size() > 0;
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
























    // Iterators for different storage cells of the mesh
    EntityIterator elementBegin()             { return _elements.end(); }
    EntityIterator domainBoundaryFaceBegin()  { return _domainboundaryfaces.end(); }
    EntityIterator processBoundaryFaceBegin() { return _processboundaryfaces.end(); }
    EntityIterator edgeBegin()                { return _edges.end(); }
    EntityIterator vertexBegin()              { return _points.end(); }

    EntityIterator elementEnd()             { return _elements.begin(); }
    EntityIterator domainBoundaryFaceEnd()  { return _domainboundaryfaces.end(); }
    EntityIterator processBoundaryFaceEnd() { return _processboundaryfaces.end(); }
    EntityIterator edgeEnd()                { return _edges.end(); }
    EntityIterator vertexEnd()              { return _points.end(); }

    /** Return pointer to single TetMesh instance. */
    static TetMesh* get_instance() { return instance_; }


protected:

    // Writes debug info to the command line
    // TODO: Use IFDEF to manipulate between no output, all output, or only master process output
    void print_debug(std::string s)
    {
        if (verbose) { std::cout << "Process_" << _mpihelper.rank() << ": " << s << std::endl; }
    }

    // Inserts a boundary face into the mesh
    // TODO: Must throw error if boundary not found matching corners in associated element
    virtual void insertBoundaryFace(Dune::GeometryType gt, id_t assoc_element_id, std::vector<int> vertexIds, int order, bool isDomainBoundary, int physicalTag)
    {
    	// Get corners of this face
    	// **********************************************************************************
    	std::vector<id_t> faceCorners = getEntityCorners(gt, vertexIds, order);
    	FemaxxFaceMapKey thisFaceKey;
    	thisFaceKey.node0 = faceCorners[0];
    	thisFaceKey.node1 = faceCorners[1];
    	thisFaceKey.node2 = faceCorners[2];

		// Sort in ascending order
		if (thisFaceKey.node0 > thisFaceKey.node1)  { std::swap(thisFaceKey.node0, thisFaceKey.node1); }
		if (thisFaceKey.node1 > thisFaceKey.node2)  { std::swap(thisFaceKey.node1, thisFaceKey.node2); }
		if (thisFaceKey.node0 > thisFaceKey.node1)  { std::swap(thisFaceKey.node0, thisFaceKey.node1); }


		// Take associated element, get all its corners, get all keys, compare to face key
		// **********************************************************************************
		std::vector<id_t> elementCorners = getEntityCorners(
				_elements[assoc_element_id].geometryType,
				_elements[assoc_element_id].vertexIds,
				_elements[assoc_element_id].interpOrder
		);
		std::vector<std::vector<int> > elementInternalIndicesFaces = internalFacesOfTetrahedron();

		for (int j = 0; j < elementInternalIndicesFaces.size(); j++)
		{
			// Define (key = sorted localIds of corners)
			FemaxxFaceMapKey thisKey;
			thisKey.node0 = elementCorners[elementInternalIndicesFaces[j][0]];
			thisKey.node1 = elementCorners[elementInternalIndicesFaces[j][1]];
			thisKey.node2 = elementCorners[elementInternalIndicesFaces[j][2]];

			// Sort in ascending order
			if (thisKey.node0 > thisKey.node1)  { std::swap(thisKey.node0, thisKey.node1); }
			if (thisKey.node1 > thisKey.node2)  { std::swap(thisKey.node1, thisKey.node2); }
			if (thisKey.node0 > thisKey.node1)  { std::swap(thisKey.node0, thisKey.node1); }

			// By comparison find internalIndex of this face
			if (thisKey == thisFaceKey)
			{

		    	// Store domain and process boundaries separately for faster iterators
				// Store Map (key -> faceId)
		        if (isDomainBoundary)  { _boundaryInternalMap[thisKey] = _faces.size(); }
		        else                   { _processInternalMap[thisKey] = _faces.size(); }

				// Store Vector (faceId -> associated element)
				FemaxxSubentityStorage thisFaceAsSubentity;
				thisFaceAsSubentity.global_id  = 0;				  // At this stage the globalId is not known yet
				thisFaceAsSubentity.element_id = assoc_element_id;
				thisFaceAsSubentity.element_internalIndex = j;
				thisFaceAsSubentity.physicalTag = physicalTag;    // Here physical tag is very important as it need not match the tag of the element

				_faces.push_back(thisFaceAsSubentity);
			}
		}
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
    	int cornerNo = 4;
    	//cornerNumber = ReferenceElements::general(geomType).size( geomType.dim() );

    	// Get corners
    	for (int j = 0; j < cornerNo; j++) {
    		int internalId = Dune::CurvilinearElementInterpolator::cornerID(gt, order, j );
    		corners.push_back(vertexIds[internalId]);
    	}

    	return corners;
	}

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

private: // Private members

    bool _verbose;
    bool _withGhostElements;


    // Storage necessary for user access and computation of globalIndex
    int _totalVerticesPerMesh;
    int _totalEdgesPerMesh;
    int _totalFacesPerMesh;
    int _totalElementsPerMesh;

    // Stores vertex coordinates and globalIds
    std::vector<int> _pointGlobalIds;
    std::vector<Vertex> _points;

    // Stores all necessary data about edges, faces and elements
    // Minimalism - edges and faces do not store interpolatory vertex ids, but refer to an element subentity
    std::vector<FemaxxEntityStorage> _elements;
    std::vector<FemaxxEntityStorage> _ghostElements;
    std::vector<FemaxxSubentityStorage> _edges;
    std::vector<FemaxxSubentityStorage> _faces;

    // List of localIds of all faces of different structural types
    // Note: Unnecessary, but speeds up iterators
    std::vector<int> _internalFaceIds;
    std::vector<int> _domainBoundaryFaceIds;
    std::vector<int> _processBoundaryFaceIds;

    // process boundary local (in terms of _processBoundaryFaceIds) id -> Ghost element localId
    // If mesh has no hanging nodes, this array would be = {0, 1, 2, 3, 4, 5, ... etc}, otherwise different elements may point to the same ghost
    std::vector<int> _processBoundaryGhostIds;

    // Temporary maps necessary to
    EdgeKey2EdgeIdMap _edgemap;
    FaceKey2FaceIdMap _internalInternalMap;
    FaceKey2FaceIdMap _boundaryInternalMap;
    FaceKey2FaceIdMap _processInternalMap;


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
