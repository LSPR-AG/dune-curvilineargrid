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

//#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>


/* ***************************************************************************
 * Specifications: Parallel Curvilinear Mesh Manager
 *
 * Existent Functionality:
 *     - All vertices, edges, faces, elements are accessible through their globalIndex and as subentities
 *     - All elements and faces possess physicalTag - integer associated to their material property
 *     - InterProcessBoundaries automatically generate
 *     - GhostElements automatically communicated and generated (optional)
 *     - GlobalIndex automatically generated for edges and faces
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
 *
 *
 *
 *
 * Development log
 *  - [FIXME] Need global index for elements since globalIndex Element+Boundary is not sufficient to separate
 *
 *  - [FIXME] Need to rewrite OCTree so it fits new convention
 *  - [FIXME] Need to introduce parallel logging through verbose as in other classes
 *  - [FIXME] Need to add normal and outerNormal
 *  - [FIXME] Need to wrap for Dune
 *  - [FIXME] Need to match Dune's internal subentity id convention
 *  - [FIXME] When returning Ghost elements, must check if they are defined, and throw error if not
 *
 *
 *
 *  ***************************************************************************/




namespace Dune {









// Wraps the storage types used in the mesh
template <class ct>
struct CurvilinearGridStorage
{

  struct Vertex
  {
    int globalIndex;
    Dune::FieldVector<ct, 3> coord;
  };

  struct Edge
  {
      int globalIndex;
      int elementIndex;
      int subentityIndex;
  };

  // Face stores indices to 2 intersecting elements, and subentity index for the first one
  // Note: element1Index is always an internal element index
  // Note: element2Index is an internal element index only for internal faces, it is a ghost element index for process boundaries, and it is -1 for domain boundaries.
  struct Face
  {
      int globalIndex;
      int element1Index;
      int element2Index;
      int element1SubentityIndex;
      int physicalTag;
  };

  struct Entity
  {
      Dune::GeometryType geometryType;
      int globalIndex;
      std::vector<int> vertexIndexSet;
      int interpOrder;
      int physicalTag;
  };

};

// Wraps the map keys necessary to construct Id for entities with minimal memory usage
struct CurvilinearEntityMapKey
{
    // This is a minimal info necessary to recognize an edge among all processes
    // The sorted globalId's of vertices of an edge
    struct EdgeKey
    {
        int node0;
        int node1;

        void sort() {
            if (node0 > node1)  { std::swap(node0, node1); }
        }

        // Allows comparing keys. Not automatically defined by C++
        bool operator==(const EdgeKey& A) const
        {
            return (A.node0 == node0) && (A.node1 == node1);
        }

        // Necessary to construct a set/map of keys
        // Compare with priority on the lower registers
        bool operator<(const EdgeKey& A) const
        {
            if (A.node0 == node0)  { return A.node1 < node1; }
            else                   { return A.node0 < node0; }
        }
    };

    // This is a minimal info necessary to recognize a face among all processes
    // The sorted globalId's of vertices of an edge
    struct FaceKey
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

        // Allows comparing keys. Not automatically defined by C++
        bool operator==(const FaceKey& A) const
        {
            return (A.node0 == node0) && (A.node1 == node1) && (A.node2 == node2);
        }

        // Necessary to construct a set/map of keys
        // Compare with priority on the lower registers
        bool operator<(const FaceKey& A) const
        {
            if (A.node0 == node0)
            {
            	if (A.node1 == node1) { return A.node2 < node2; }
            	else { return A.node1 < node1; }
            }
            else { return A.node0 < node0; }
        }
    };
};


// Enumerates the structural type of faces of the mesh
struct CurvilinearGridFaceType
{
  enum {
    Internal = 0,
    DomainBoundary = 1,
    ProcessBoundary = 2
  };
};




template <class ct>
class CurvilinearGridBase {
public:

    /* public types */
    typedef Dune::FieldVector<ct, 3>               Vertex;
    typedef std::vector<Vertex>                    VertexVector;

    typedef typename Dune::CurvilinearGridStorage<ct>::Vertex     VertexStorage;
    typedef typename Dune::CurvilinearGridStorage<ct>::Edge       EdgeStorage;
    typedef typename Dune::CurvilinearGridStorage<ct>::Face       FaceStorage;
    typedef typename Dune::CurvilinearGridStorage<ct>::Entity     EntityStorage;

    typedef CurvilinearEntityMapKey::EdgeKey       EdgeKey;
    typedef CurvilinearEntityMapKey::FaceKey       FaceKey;

    typedef std::map<EdgeKey, int>                 EdgeKey2EdgeIdMap;
    typedef std::map<FaceKey, int>                 FaceKey2FaceIdMap;

    typedef Dune::CurvilinearGeometry<ct, 1, 3>    EdgeGeometry;
    typedef Dune::CurvilinearGeometry<ct, 2, 3>    FaceGeometry;
    typedef Dune::CurvilinearGeometry<ct, 3, 3>    ElementGeometry;

    typedef std::map<int, int> Index2IndexMap;



public: /* public methods */

    /** Default constructor - DO NOT USE*/
    CurvilinearGridBase() { }

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridBase(bool withGhostElements, bool verbose,  MPIHelper &mpihelper ) :
    	nVertexTotal_(0),
        nEdgeTotal_(0),
        nFaceTotal_(0),
        nElementTotal_(0),
        withGhostElements_(withGhostElements),
        verbose_(verbose),
        octree_(0),
        mpihelper_(mpihelper)
    {
        //assert(instance_ == 0);
        //instance_ = this;

        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();
    }

    /** Destructor */
    ~CurvilinearGridBase() {
        // ** delete octree **
        // if (octree_)  { delete octree_; }
        //instance_ = 0;
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
        VertexStorage point;
        point.coord = p;
        point.globalIndex = globalIndex;

        vertexGlobal2LocalMap_[globalIndex] = point_.size();
        point_.push_back(point);

    }

    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] globalId         the global index of this element
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * Note: Even though we pass the globalId as a parameter from GMSH, it is a globalIndex for the set Elements+BoundaryFaces,
     * therefore obtaining globalIndex for elements from it is not possible, it will have to be communicated
     *
     * */
    void insertElement(Dune::GeometryType gt, int globalId, const std::vector<int> & vertexIndexSet, int order, int physicalTag)
    {
        if (!gt.isTetrahedron() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
            DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertElement() unexpected element type or number of interpolatory points");
        }

        EntityStorage thisElement;

        thisElement.geometryType = gt;
        thisElement.globalIndex = 0;        // At this stage globalIndex is not known yet
        thisElement.interpOrder = order;
        thisElement.physicalTag = physicalTag;
        thisElement.vertexIndexSet = vertexIndexSet;

        elementGlobal2LocalMap_[globalId] = element_.size();
        element_.push_back(thisElement);
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
        if (!gt.isTriangle() || (vertexIndexSet.size() != Dune::CurvilinearGeometryHelper::dofPerOrder(gt, order)))  {
            DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertBoundarySegment() unexpected number of interpolatory points");
        }


        // Get corners of this face
        // **********************************************************************************
        std::vector<int> faceCorners = getEntityCorners(gt, vertexIndexSet, order);
        FaceKey thisFaceKey;
        thisFaceKey.node0 = point_[faceCorners[0]].globalIndex;
        thisFaceKey.node1 = point_[faceCorners[1]].globalIndex;
        thisFaceKey.node2 = point_[faceCorners[2]].globalIndex;

        // Sort in ascending order
        thisFaceKey.sort();


        // Take associated element, get all its corners, get all keys, compare to face key
        // **********************************************************************************
        std::vector<int> elementCorners = getEntityCorners(
                element_[associatedElementIndex].geometryType,
                element_[associatedElementIndex].vertexIndexSet,
                element_[associatedElementIndex].interpOrder
        );



        // Search for the face among subentities of the element
        // **********************************************************************************
        int j = 0;
        int nFacePerTetrahedron = 4;
        bool found_face = false;

        while (!found_face)
        {
            if (j == nFacePerTetrahedron)  {
                DUNE_THROW(Dune::IOError, "CurvilinearGrid: insertBoundarySegment() did not find the face in the associated element");
            }

            std::vector<int> internalLinearSubentityIndices = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(element_[associatedElementIndex].geometryType, 1, j);

            // Define (key = sorted localIndices of corners)
            FaceKey thisKey;
            thisKey.node0 = point_[elementCorners[internalLinearSubentityIndices[0]]].globalIndex;
            thisKey.node1 = point_[elementCorners[internalLinearSubentityIndices[1]]].globalIndex;
            thisKey.node2 = point_[elementCorners[internalLinearSubentityIndices[2]]].globalIndex;


            // Sort in ascending order
            thisKey.sort();

            // By comparison find internalIndex of this face
            if (thisKey == thisFaceKey)
            {
                found_face = true;

                // Store domain and process boundaries separately for faster iterators
                // Store Map (key -> domainBoundaryIndex), Vector (domainBoundaryIndex -> faceIndex)

                boundaryInternalMap_[thisKey] = domainBoundaryFaceIndex_.size();
                domainBoundaryFaceIndex_.push_back(face_.size());

                // Store Vector (faceId -> associated element)
                FaceStorage thisFaceAsSubentity;
                thisFaceAsSubentity.globalIndex  = 0;                  // At this stage the globalId is not known yet
                thisFaceAsSubentity.element1Index = associatedElementIndex;
                thisFaceAsSubentity.element2Index = -1;              // Boundary Segments do not have a 2nd neighbor
                thisFaceAsSubentity.element1SubentityIndex = j;
                thisFaceAsSubentity.physicalTag = physicalTag;    // Here physical tag is very important as it need not match the tag of the element

                face_.push_back(thisFaceAsSubentity);
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
    void generateMesh(int nVertexTotal, int nElementTotal) {
        typedef FaceKey2FaceIdMap::iterator IterType;

        // Store the total number of entities obtained from the reader
    	nVertexTotal_ = nVertexTotal;
    	nElementTotal_ = nElementTotal;

        // Construct missing parts of the mesh
        generateEdges();
        generateFaces();
#if HAVE_MPI
        generateProcessBoundaryNeighbors();
        if (withGhostElements_) { generateGhostElements(); }
        generateGlobalIndices();
#else
        // In sequential case:
        // * Boundary Neighbors not necessary, since all boundaries are domain boundaries
        // * No ghost elements, even if requested by user
        // * Fake globalIndex by making it equal to localIndex

        withGhostElements_ = false;
        for (int i = 0; i < edge_.size(); i++)     { edge_[i].globalIndex = i;     edgeGlobal2LocalMap_[i] = i; }
        for (int i = 0; i < face_.size(); i++)     { face_[i].globalIndex = i;     faceGlobal2LocalMap_[i] = i; }
        for (int i = 0; i < element_.size(); i++)  { element_[i].globalIndex = i;  elementGlobal2LocalMap_[i] = i; }
#endif
        //constructoctree_();


        // Deletes all temporary memory
        edgemap_.clear();

        internalFaceIndex_.reserve(internalInternalMap_.size());
        domainBoundaryFaceIndex_.reserve(boundaryInternalMap_.size());
        processBoundaryFaceIndex_.reserve(processInternalMap_.size());

        for (IterType it = internalInternalMap_.begin(); it != internalInternalMap_.end();  it++) { internalFaceIndex_.push_back((*it).second); }
        for (IterType it = boundaryInternalMap_.begin(); it != boundaryInternalMap_.end();  it++) { domainBoundaryFaceIndex_.push_back((*it).second); }
        for (IterType it = processInternalMap_.begin();  it != processInternalMap_.end();   it++) { processBoundaryFaceIndex_.push_back((*it).second); }

        internalInternalMap_.clear();
        boundaryInternalMap_.clear();
        processInternalMap_.clear();
    }

    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    /** Get total number of entities in a mesh  */
    int nVertexTotal()   { return nVertexTotal_; }
    int nEdgeTotal()     { return nEdgeTotal_; }
    int nFaceTotal()     { return nFaceTotal_; }
    int nElementTotal()  { return nElementTotal_; }

    /** Get total number of entities on this process  */
    int nVertex()   { return point_.size(); }
    int nEdge()     { return edge_.size(); }
    int nFace()     { return face_.size(); }
    int nElement()  { return element_.size(); }


    /** Vertex coordinate
     *  \param[in] localIndex            local vertex index (insertion index)
     * */
    Vertex vertex(int localIndex) const { return point_[localIndex]; }


    /** Storage data related to this edge, except of explicit vertex coordinates
     *  \param[in] localIndex            local vertex index (insertion index)
     *
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    EntityStorage edgeData(int localIndex)
    {
        EntityStorage & assocElement = element_[edge_[localIndex].elementIndex];
        EntityStorage thisEdge;
        thisEdge.geometryType.makeLine();
        thisEdge.globalIndex = edge_[localIndex].globalIndex;
        thisEdge.interpOrder = assocElement.interpOrder;
        thisEdge.physicalTag = -1;        // Note: Edges do not have a physical tag

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<int> subentityVertexIndices = Dune::CurvilinearElementInterpolator<ct, 1, 3>::subentityInternalCoordinates(assocElement.geometryType, thisEdge.interpOrder, 2, edge_[localIndex].subentityIndex );

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisEdge;

    }


    /** Storage data related to this face, except of explicit vertex coordinates
     *  \param[in] localIndex            local vertex index (insertion index)
     *
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    EntityStorage faceData(int localIndex)
    {
        EntityStorage assocElement = element_[face_[localIndex].element1Index];
        EntityStorage thisFace;
        thisFace.geometryType.makeTriangle();
        thisFace.globalIndex = face_[localIndex].globalIndex;
        thisFace.interpOrder = assocElement.interpOrder;
        thisFace.physicalTag = face_[localIndex].physicalTag;

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<int> subentityVertexIndices = Dune::CurvilinearElementInterpolator<ct, 2, 3>::subentityInternalCoordinates(assocElement.geometryType, thisFace.interpOrder, 1, face_[localIndex].element1SubentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisFace;
    }


    /** Storage data related to this element, except of explicit vertex coordinates
     *  \param[in] localIndex            local vertex index (insertion index)
     *
     * */
    EntityStorage elementData(int localIndex) { return element_[localIndex]; }


    // Construct and retrieve geometry class associated to this edge
    EdgeGeometry edgeGeometry(int localIndex)
    {
        EntityStorage edgeData = edgeData(localIndex);
        return entityGeometry<1>(edgeData.geometryType, edgeData.vertexIndexSet, edgeData.interpOrder);
    }

    // Construct and retrieve geometry class associated to this face
    FaceGeometry faceGeometry(int localIndex)
    {
        EntityStorage faceData = faceData(localIndex);
        return entityGeometry<2>(faceData.geometryType, faceData.vertexIndexSet, faceData.interpOrder);
    }

    // Construct and retrieve geometry class associated to this element
    ElementGeometry elementGeometry(int localIndex)
    {
        EntityStorage elementData = elementData(localIndex);
        return entityGeometry<3>(elementData.geometryType, elementData.vertexIndexSet, elementData.interpOrder);
    }

    /** Return pointer to Octree or 0 if it is not constructed. */
    //LooseOctree<Tet>* getoctree_() const { return octree_; }


    /** Compute center and extent (halved) of the bounding box of the whole mesh.
     *  \param[in] center            output: coordinate of the center of the box
     *  \param[in] extent            output: coordinate-displacement from the center
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    void meshBoundingBox(Vertex & center, Vertex & extent) const
    {
        Vertex min = point_[0];
        Vertex max = min;

        for (int i = 1; i < point_.size(); i ++) {
            min[0] = std::min(min[0], point_[i][0]);
            min[1] = std::min(min[1], point_[i][1]);
            min[2] = std::min(min[2], point_[i][2]);

            max[0] = std::max(max[0], point_[i][0]);
            max[1] = std::max(max[1], point_[i][1]);
            max[2] = std::max(max[2], point_[i][2]);
        }
        center = min + max;  center *= 0.5;
        extent = max - min;  extent *= 0.5;
    }


    /** Return edge connecting node0 and node1. A runtime_error is
     *  thrown if the edge does not exist.
     */

    // FIXME: The maps will be deleted to save space. If the lookup functionality is not necessary it should just be deleted. Check which routines still use it
    int lookupEdge(int node0, int node1) {
        EdgeKey thisEdgeKey;
        thisEdgeKey.node0 = node0;
        thisEdgeKey.node1 = node1;

        if (thisEdgeKey.node0 > thisEdgeKey.node1)  { std::swap(thisEdgeKey.node0, thisEdgeKey.node1); }

        EdgeKey2EdgeIdMap::iterator iter = edgemap_.find(thisEdgeKey);
        if (iter == edgemap_.end())  { throw std::runtime_error("Accessed non-existing edge."); }

        return (*iter).second;
    }

    /** Return face connecting by node0, node1 and node2. A runtime_error is
     *  thrown if the face does not exist.
     */
    int lookupFace(int node0, int node1, int node2){
        FaceKey thisFaceKey;
        thisFaceKey.node0 = node0;
        thisFaceKey.node1 = node1;
        thisFaceKey.node2 = node2;

        // sort node ids in ascending order
        if (thisFaceKey.node0 > thisFaceKey.node1)  { std::swap(thisFaceKey.node0, thisFaceKey.node1); }
        if (thisFaceKey.node1 > thisFaceKey.node2)  { std::swap(thisFaceKey.node1, thisFaceKey.node2); }
        if (thisFaceKey.node0 > thisFaceKey.node1)  { std::swap(thisFaceKey.node0, thisFaceKey.node1); }

        // Check if the face is present in any of the 3 face maps, otherwise throw error
        FaceKey2FaceIdMap::iterator iterInt = internalInternalMap_.find(thisFaceKey);

        if (iterInt != internalInternalMap_.end()) { return (*iterInt).second; }
        else
        {
            FaceKey2FaceIdMap::iterator iterProc = processInternalMap_.find(thisFaceKey);
            if (iterProc != processInternalMap_.end()) { return (*iterProc).second; }
            else
            {
                FaceKey2FaceIdMap::iterator iterBnd = boundaryInternalMap_.find(thisFaceKey);
                if (iterBnd != processInternalMap_.end()) { return (*iterBnd).second; }
                else { throw std::runtime_error("Accessed non-existing face."); }
            }
        }
    }


    /** Return true if p is in inside the mesh, i.e. inside at least one tet of the mesh. */
    bool is_inside(const Vertex & p) const {
        std::vector<int> containerIds;
        find_tets_by_point(p, containerIds);
        return containerIds.size() > 0;
    }

    // Iterators over local indices of the mesh
    // NOTE: There are no iterators over entities because there is no entity object in the core mesh
    // There will be generic entity in the wrapper because the wrapper will define an entity object
    Index2IndexMap::iterator elementIndexBegin()  { return elementGlobal2LocalMap_.begin(); }
    Index2IndexMap::iterator faceIndexBegin()     { return faceGlobal2LocalMap_.begin(); }
    Index2IndexMap::iterator edgeIndexBegin()     { return edgeGlobal2LocalMap_.begin(); }
    Index2IndexMap::iterator vertexIndexBegin()   { return vertexGlobal2LocalMap_.begin(); }

    Index2IndexMap::iterator elementIndexEnd()  { return elementGlobal2LocalMap_.end(); }
    Index2IndexMap::iterator faceIndexEnd()     { return faceGlobal2LocalMap_.end(); }
    Index2IndexMap::iterator edgeIndexEnd()     { return edgeGlobal2LocalMap_.end(); }
    Index2IndexMap::iterator vertexIndexEnd()   { return vertexGlobal2LocalMap_.end(); }

    // This construction allows fast iteration over faces of different type
    std::vector<int>::iterator faceInternalBegin()         { return internalFaceIndex_.begin(); }
    std::vector<int>::iterator faceDomainBoundaryBegin()   { return domainBoundaryFaceIndex_.begin(); }
    std::vector<int>::iterator faceProcessBoundaryBegin()  { return processBoundaryFaceIndex_.begin(); }

    std::vector<int>::iterator faceInternalEnd()         { return internalFaceIndex_.end(); }
    std::vector<int>::iterator faceDomainBoundaryEnd()   { return domainBoundaryFaceIndex_.end(); }
    std::vector<int>::iterator faceProcessBoundaryEnd()  { return processBoundaryFaceIndex_.end(); }



    /** Return pointer to single CurvilinearGridBase instance. */
    //static CurvilinearGridBase* get_instance() { return instance_; }


protected:


    /* ***************************************************************************
     * Section: Auxiliary Methods
     * ***************************************************************************/

    // Writes debug info to the command line
    // TODO: Use IFDEF to manipulate between no output, all output, or only master process output
    void print_debug(std::string s)
    {
        if (verbose_) { std::cout << "Process_" << rank_ << ": " << s << std::endl; }
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
    Dune::CurvilinearGeometry<ct, mydim, 3> entityGeometry(
            Dune::GeometryType gt,
            std::vector<int> & vertexIndexSet,
            int order) const
    {
        VertexVector entityVertices;
        for (int i = 0; i < vertexIndexSet.size(); i++) { entityVertices.push_back(point_[vertexIndexSet[i]]); }

        return Dune::CurvilinearGeometry<ct, mydim, 3> (gt, entityVertices, order);
    }

    // Returns corner id's of this entity
    // TODO: Move functionality to the curvilinear geometry
    std::vector<int> getEntityCorners(
            Dune::GeometryType gt,
            const std::vector<int> & vertexIndexSet,
            int order) const
    {
        std::vector<int> corners;

        // Get corner number
        int cornerNo = gt.dim() + 1;
        //cornerNumber = ReferenceElements::general(geomType).size( geomType.dim() );

        // Get corners
        for (int j = 0; j < cornerNo; j++) {
            int internalId = Dune::CurvilinearGeometryHelper::cornerID(gt, order, j );
            corners.push_back(vertexIndexSet[internalId]);
        }

        return corners;
    }



    /* ***************************************************************************
     * Section: Constructing the mesh
     * ***************************************************************************/

    // Generates all edges
    // FIXME: Original code seems to give edges some orientation
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: Use more generic functions when extending to any mesh other than tetrahedral
    void generateEdges()
    {
        // Loop over all elements and their edges
        for (int i = 0; i < element_.size(); i++)
        {
        	int nEdgePerTetrahedron = 6;
        	std::vector<int> elementCorners = getEntityCorners(element_[i].geometryType, element_[i].vertexIndexSet, element_[i].interpOrder);

            for (int j = 0; j < nEdgePerTetrahedron; j++)
            {
                std::vector<int> internalLinearSubentityIndices = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(element_[i].geometryType, 2, j);

            	// Define (key = sorted localIndices of corners)
                EdgeKey thisKey;
                thisKey.node0 = elementCorners[internalLinearSubentityIndices[0]];
                thisKey.node1 = elementCorners[internalLinearSubentityIndices[1]];

                // Sort in ascending order
                thisKey.sort();

                // If this edge has not been added already, add it
                if (edgemap_.find(thisKey) == edgemap_.end())
                {
                    // Store map (key -> edgeId)
                    edgemap_[thisKey] = edge_.size();

                    // Store vector (edgeId -> elemId + edgeElemIndex)
                    // Note: Edges do not have physical tag at all so we do not even store it
                    EdgeStorage thisEdge;
                    thisEdge.globalIndex = 0;        // GlobalId for edge determined later using global communication
                    thisEdge.elementIndex = i;
                    thisEdge.subentityIndex = j;
                    edge_.push_back(thisEdge);
                }


            }

        }
    }

    // Generates Internal and ProcessBoundary Faces. (!!!) Assumes that all Domain Boundary Faces have been added.
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
    void generateFaces()
    {
        typedef std::map<FaceKey, std::vector<int>> tmpFace2InfoMap;
        tmpFace2InfoMap tmpFaceMap;

        // Loop over all elements and their faces
        for (int i = 0; i < element_.size(); i++)
        {
        	int nFacePerTetrahedron = 6;
            std::vector<int> elementCorners = getEntityCorners(element_[i].geometryType, element_[i].vertexIndexSet, element_[i].interpOrder);

            // Store info for all faces except of domain boundaries
            // Store it in a map, not to store internal faces twice
            for (int j = 0; j < nFacePerTetrahedron; j++)
            {
            	std::vector<int> internalLinearSubentityIndices = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(element_[i].geometryType, 1, j);

                // Define (key = sorted localIndices of corners)
                FaceKey thisKey;
                thisKey.node0 = point_[elementCorners[internalLinearSubentityIndices[0]]].globalIndex;
                thisKey.node1 = point_[elementCorners[internalLinearSubentityIndices[1]]].globalIndex;
                thisKey.node2 = point_[elementCorners[internalLinearSubentityIndices[2]]].globalIndex;

                // Sort in ascending order
                thisKey.sort();

                // Only do sth if this is not the Domain Boundary
                if (boundaryInternalMap_.find(thisKey) == boundaryInternalMap_.end())
                {
                    std::vector<int> connectedFaceInfo;
                    tmpFace2InfoMap::iterator iter = tmpFaceMap.find(thisKey);

                    if (iter == tmpFaceMap.end()) { connectedFaceInfo = (*iter).second; }

                    connectedFaceInfo.push_back(i);
                    connectedFaceInfo.push_back(j);

                    tmpFaceMap[thisKey] = connectedFaceInfo;
                }
            }

            // Add internal and process boundary faces to the mesh
            for (tmpFace2InfoMap::iterator iter = tmpFaceMap.begin();  iter != tmpFaceMap.end();  iter++)
            {
                FaceStorage thisFace;
                std::vector<int> connectedFaceInfo = (*iter).second;

                thisFace.globalIndex = 0;       // GlobalId is defined at a later stage
                thisFace.element1Index = connectedFaceInfo[0];
                thisFace.element1SubentityIndex = connectedFaceInfo[1];
                thisFace.physicalTag = -1;    // At the moment physicalTag of an internal face is not defined as it could be inbetween two different elements

                // Store internal, domain and process boundaries separately for faster iterators
                if (connectedFaceInfo.size() == 2)
                {
                    // Store Map (key -> processBoundaryFaceIndex), Vector (processBoundaryFaceIndex -> faceIndex)
                    processInternalMap_[(*iter).first] = processBoundaryFaceIndex_.size();
                    processBoundaryFaceIndex_.push_back(face_.size());

                    thisFace.element2Index = 0;    // Eventually this will be the Ghost Element Index
                }
                else
                {
                    // Store Map (key -> internalFaceIndex), Vector (internalFaceIndex -> faceIndex)
                    internalInternalMap_[(*iter).first] = internalFaceIndex_.size();
                    internalFaceIndex_.push_back(face_.size());

                    thisFace.element2Index = connectedFaceInfo[2];  // This is the 2nd neighbor of this internal face
                }

                face_.push_back(thisFace);
            }
        }
    }


    /** Construct octree for locating tetrahedrons in mesh */
    /*
    void constructoctree_() {
        rAssert(octree_ == 0);

        // bounding box of whole mesh
        Vector3 center, extent;
        meshBoundingBox(center, extent);

        // octree length is the largest component of extent
        double length = extent.x;
        if (extent.y > length)
            length = extent.y;
        if (extent.z > length)
            length = extent.z;

        // construct LooseOctree with large max depth
        octree_ = new LooseOctree<Tet>(center, length, 100);

        // loop over all tets and insert them in the octree
        for (int t = 0; t < get_nof_tets(); t ++) {
            octree_->add_node(get_tet(t));
        }
        int max_depth, nof_octants, nof_nodes;
        double avg_node_depth;
        octree_->get_statistics(max_depth, avg_node_depth, nof_octants, nof_nodes);
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
        octree_->find_nodes_by_point(p, nodes, nof_visited,
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
        octree_->find_nodes_by_point(p, nodes, nof_visited,
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
        octree_->find_nodes_by_point(p, nodes, nof_visited,
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
        OctreeNode* node = octree_->find_one_node_by_point(p, nof_visited, (LooseOctree::filter)(&Tet::is_inside));
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
            get_tet(6488)->elementBoundingBox(center, extent);
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
            octree_->add_node(get_tet(6488))->dump(cout);
            cout << "eps_mach=" << std::numeric_limits<double>::epsilon() << "\n";
        }
    #endif
    }
    */


#if HAVE_MPI

    // [TODO] Possibly inefficient implementation:
    // 1) Excessive communication. Package sent to all can only be used by 1 other process
    // 2) Imbalance bottleneck. All processes have to wait until the ones with biggest number of procBoundaries communicate
    void generateProcessBoundaryNeighbors()
    {
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // 0) Make space in the neighbor process rank vector
        processBoundaryNeighborProcess_.reserve(processInternalMap_.size());

        // 1) Each process must somehow learn which other process it shares each process boundary with
        // 1.1) Every process tells everyone how many interprocessor boundaries he has. Calculate max = max(procbnd per process), loops 1 to max
        int thisProcBoundarySize = processInternalMap_.size();
        int maxProcBoundarySize = collective_comm.max(thisProcBoundarySize);


        int thisProcFace = 0;
        FaceKey2FaceIdMap::iterator iter = processInternalMap_.begin();

        // 1.2) Per each loop allgather - send 1 face global key from each to all. If all faces already sent, send empty face
        for (int thisProcFace = 0; thisProcFace < maxProcBoundarySize; thisProcFace++)
        {
            FaceKey thisFaceKey;
            if (thisProcFace < thisProcBoundarySize)  { thisFaceKey = (*iter).first; iter++; }
            else {
                // All process boundaries have been sent by this process. Send fake faces now.
                // Do NOT increase iterator
                thisFaceKey.node0 = -1;
                thisFaceKey.node1 = -1;
                thisFaceKey.node2 = -1;
            }

            std::vector<FaceKey> keyset (size_);

            collective_comm.allgather(&thisFaceKey, 1, reinterpret_cast<FaceKey*> (keyset.data()) );

            // 1.2.2) Loop over all received faces. If a face is present, note its sender rank
            for (int iProc = 0; iProc < keyset.size(); iProc++)
            {
                // Skip the face if it is fake, or if it is the one sent by this same process
                if ((iProc != rank_) && (keyset[iProc].node0 >= 0))
                {
                    FaceKey2FaceIdMap::iterator iter2 = processInternalMap_.find(keyset[iProc]);

                    // If this face is present, note its sender process
                    if (iter2 != processInternalMap_.end()) { processBoundaryNeighborProcess_[(*iter2).second] = iProc; }
                }
            }
        }
    }


    // NOTE: Call after generating Ghost Elements
    void generateGlobalIndices()
    {
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // 1) Get edges and faces on this process that are not owned by this process
        // *************************************************************************
        EdgeKey2EdgeIdMap edgesNonOwned;  globalComputeNonOwnedEdges(edgesNonOwned);
        FaceKey2FaceIdMap facesNonOwned;  globalComputeNonOwnedFaces(facesNonOwned);

        int edgesOwned = edge_.size() - edgesNonOwned.size();
        int facesOwned = face_.size() - facesNonOwned.size();
        int elementsOwned = element_.size();


        // 2) Communicate number of edges and faces owned by each process to all
        // *************************************************************************
        std::vector<int> edgesOnProcess(size_);      // owned edges [rank]
        std::vector<int> facesOnProcess(size_);      // owned faces [rank]
        std::vector<int> elementsOnProcess(size_);  // owned elements [rank]

        collective_comm.allgather (&edgesOwned, 1, reinterpret_cast<int*> (edgesOnProcess.data()));
        collective_comm.allgather (&facesOwned, 1, reinterpret_cast<int*> (facesOnProcess.data()));
        collective_comm.allgather (&elementsOwned, 1, reinterpret_cast<int*> (elementsOnProcess.data()));

        int edgesBeforeMe = 0;      // Sum(edgesOwned : rank < thisRank)
        int facesBeforeMe = 0;      // Sum(facesOwned : rank < thisRank)
        int elementsBeforeMe = 0;   // Sum(elementsOwned : rank < thisRank)

        for (int i = 0; i < size_; i++)
        {
            nEdgeTotal_ += edgesOnProcess[i];
            nFaceTotal_ += facesOnProcess[i];

            if (i < rank_)
            {
                edgesBeforeMe += edgesOnProcess[i];
                facesBeforeMe += facesOnProcess[i];
                elementsBeforeMe += elementsOnProcess[i];
            }
        }

        // 3) Enumerate all edges, faces and elements that you own
        // *************************************************************************

        int iEdgeGlobalId = edgesBeforeMe;
        int iFaceGlobalId = facesBeforeMe;

        // Enumerating elements is simply shifting the local index, since all elements on this process are owned by it
        for (int i = 0; i < element_.size(); i++)  { element_[i].globalIndex = elementsBeforeMe + i; }

        // Faces that are not shared with other processes are automatically owned by this process
        for (int i = 0; i < internalFaceIndex_.size(); i++)        { face_[internalFaceIndex_[i]].globalIndex = iFaceGlobalId++; }
        for (int i = 0; i < domainBoundaryFaceIndex_.size(); i++)  { face_[domainBoundaryFaceIndex_[i]].globalIndex = iFaceGlobalId++; }

        for (FaceKey2FaceIdMap::iterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            // This process owns this face if it is in the processInternalMap and if it is not in the non-owned face map
            if (facesNonOwned.find((*iter).first) == facesNonOwned.end())  { face_[(*iter).second].globalIndex = iFaceGlobalId++; }
        }

        for (EdgeKey2EdgeIdMap::iterator iter = edgemap_.begin(); iter != edgemap_.end(); iter++)
        {
            // This process owns this edge if it is in the edgemap and if it is not in the non-owned edge map
            if (edgesNonOwned.find((*iter).first) == edgesNonOwned.end())  { edge_[(*iter).second].globalIndex = iEdgeGlobalId++; }
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
        std::vector< std::vector<int> > thisProcessGhostElementIndices (size_, std::vector<int>() );
        std::vector< std::vector<int> > thisProcessGhostFaceGlobalIndices (size_, std::vector<int>() );

        for (int iPBFace = 0; iPBFace < processBoundaryFaceIndex_.size(); iPBFace++ )  {
            int thisFaceIndex = processBoundaryFaceIndex_[iPBFace];

            int thisGhostIndex = face_[thisFaceIndex].element1Index;
            int thisFaceGlobalIndex = face_[thisFaceIndex].globalIndex;
            int thisNeighborRank = processBoundaryNeighborProcess_[iPBFace];
            thisProcessGhostElementIndices[thisNeighborRank].push_back(thisGhostIndex);
            thisProcessGhostFaceGlobalIndices[thisNeighborRank].push_back(thisFaceGlobalIndex);
        }


        // 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
        // *************************************************************************************
        std::vector< std::vector<int> > neighborProcessGhostOrder (size_, std::vector<int>() );
        ghostDistributeInterpolationOrders(thisProcessGhostElementIndices, neighborProcessGhostOrder);


        // 3) MPI_alltoallv - Package element globalIndex + elementPhysicalTag + all interpVertex globalIds
        // *************************************************************************************
        std::vector<int> packageGhostElementData;
        ghostDistributeGhostElements(thisProcessGhostElementIndices, thisProcessGhostFaceGlobalIndices, neighborProcessGhostOrder, packageGhostElementData );

        thisProcessGhostElementIndices.clear();
        thisProcessGhostFaceGlobalIndices.clear();


        // 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
        // *************************************************************************************
        std::vector<std::set<int> > missingVertices (size_, std::set<int>());
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

#endif

    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGlobalIndices()
     * ***************************************************************************/

#if HAVE_MPI
    // Compute all edges that you do not own
    // Return pairs (edge local index -> owning process)
    void globalComputeNonOwnedEdges(EdgeKey2EdgeIdMap & edgesNonOwned)
    {
        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        EdgeKey2EdgeIdMap             processBoundaryEdges;                // Edges that are shared by exactly 2 processes
        EdgeKey2EdgeIdMap             processBoundaryComplicatedEdges;     // Edges that are shared by more than 2 processes

        // Mark each process boundary edge as either simple or complicated
        // For simple edges immediately mark them as owned or non-owned
        // *************************************************************************
        for(FaceKey2FaceIdMap::iterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            FaceKey thisFace = (*iter).first;
            int neighborRank = processBoundaryNeighborProcess_[(*iter).second];

            EdgeKey thisEdge[3];
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
                    else if (neighborRank < rank_)  { edgesNonOwned[thisEdge[i]] = neighborRank; }
                }
            }
        }


        // Communicate number of own complicated edges to all processes
        // *************************************************************************
        int nThisComplicatedEdges = processBoundaryComplicatedEdges.size();
        int totalRecvKeys = 0;
        std::vector<int> processComplicatedEdgeNumbers (size_);
        std::vector<int> processComplicatedEdgeDispls (size_);
        collective_comm.allgather(&nThisComplicatedEdges, 1, reinterpret_cast<int*> (processComplicatedEdgeNumbers.data()));

        // Communicate complicate edge keys to all
        // *************************************************************************
        for (int i = 0; i < size_; i++)
        {
            totalRecvKeys += processComplicatedEdgeNumbers[i];
            processComplicatedEdgeDispls[i] = (i == 0) ? 0 : processComplicatedEdgeDispls[i - 1] + processComplicatedEdgeNumbers[i - 1];
        }

        std::vector<EdgeKey> packageProcessComplicatedEdgeKeysSend;
        std::vector<EdgeKey> packageProcessComplicatedEdgeKeysRecv(totalRecvKeys);

        for (EdgeKey2EdgeIdMap::iterator faceIter = processBoundaryComplicatedEdges.begin(); faceIter != processBoundaryComplicatedEdges.end(); faceIter++)
        {
        	packageProcessComplicatedEdgeKeysSend.push_back((*faceIter).first);
        }

        collective_comm.allgatherv(
                packageProcessComplicatedEdgeKeysSend.data(),
                processBoundaryComplicatedEdges.size(),
                reinterpret_cast<EdgeKey*>(packageProcessComplicatedEdgeKeysRecv.data()),
                processComplicatedEdgeNumbers.data(), processComplicatedEdgeDispls.data()
                );


        // Mark complicated edges as non-owned if they are on this process and also on a process with lower rank
        // We only care about processes with lower rank when checking who owns the edge
        // *************************************************************************
        int iData = 0;
        for (int i = 0; i < rank_; i++)
        {
            for (int j = 0; j < processComplicatedEdgeNumbers[i]; j++)
            {
                EdgeKey tmpEdgeKey = packageProcessComplicatedEdgeKeysRecv[iData++];
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
        for (FaceKey2FaceIdMap::iterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            int thisNeighborRank = processBoundaryNeighborProcess_[(*iter).second];
            if (thisNeighborRank < rank_)  { facesNonOwned[(*iter).first] = thisNeighborRank; }
        }
    }

    // Communicates all process boundary face global Id's to the neighbors if owned
    void globalDistributeMissingFaceGlobalId()
    {
        typedef std::pair<FaceKey, int>  FaceInfo;
        std::vector< std::vector< FaceInfo > > facesToSend (size_);

        int totalRecvSize = 0;
        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvbuf, recvcounts, rdispls;

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************8
        for (FaceKey2FaceIdMap::iterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            int localIndex = (*iter).second;
            int neighborRank = processBoundaryNeighborProcess_[localIndex];

            // If the neighbor of this face has lower rank, then add it to recv, else to send
            if (neighborRank < rank_)  { recvcounts[neighborRank]++; }
            else
            {
                int thisGlobalIndex = face_[localIndex].globalIndex;
                facesToSend[neighborRank].push_back(FaceInfo ((*iter).first, thisGlobalIndex ));
            }
        }


        // 2) Fill in communication arrays
        // ********************************************************************************************8
        for (int i = 0; i < size_; i++)
        {
            sendcounts.push_back(facesToSend[i].size());
            totalRecvSize += recvcounts[i];
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );

            for (int j = 0; j < facesToSend[i].size(); j++)
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
        for (int iProc = 0; iProc < size_; iProc++)
        {
            for (int iFace = 0; iFace < recvcounts[iProc]; iFace++)
            {
                FaceKey thisKey;
                int thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];
                thisKey.node2 = recvbuf[iData++];

                int thislocalIndex = processInternalMap_[thisKey];
                face_[thislocalIndex].globalIndex = thisGlobalId;
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
        std::vector<int> edgeNumberToSend (size_);
        std::vector< std::vector< EdgeKey > >  requestedKeys(size_);

        for (EdgeKey2EdgeIdMap::iterator iter = edgesNonOwned.begin(); iter != edgesNonOwned.end(); iter++)  { requestedKeys[(*iter).second].push_back((*iter).first); }
        for (int i = 0; i < size_; i++) { edgeNumberRequested.push_back(requestedKeys[i].size()); }

        MPI_Alltoall(edgeNumberRequested.data(), 1, MPI_INT, reinterpret_cast<int*>(edgeNumberToSend.data()), 1, MPI_INT, comm);

        // 2) MPI_alltoall - send each process the edge keys you want from it
        // **************************************************************************
        int totalRecvData = 0;
        std::vector<int> packageEdgesRequested, sendcounts, sdispls;
        std::vector<int> packageEdgesToSend,    recvcounts, rdispls;


        for (int iProc = 0; iProc < size_; iProc++) {
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
           for (int iProc = 0; iProc < size_; iProc++) {
            for (int iFace = 0; iFace < edgeNumberToSend[iProc]; iFace++)
            {
                EdgeKey thisKey;
                thisKey.node0 = packageEdgesToSend[iData++];
                thisKey.node1 = packageEdgesToSend[iData++];

                packageEdgeGlobalIndicesSend.push_back(edge_[edgemap_[thisKey]].globalIndex);
            }

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + edgeNumberToSend[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + edgeNumberRequested[iProc-1] );
           }

           MPI_Alltoallv (packageEdgeGlobalIndicesSend.data(), edgeNumberToSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(packageEdgeGlobalIndicesRecv.data()), edgeNumberRequested.data(), rdispls.data(), MPI_INT, comm );


        // 4) Mark all missing edges
        // **************************************************************************
        iData = 0;
           for (int iProc = 0; iProc < size_; iProc++) {
            for (int iFace = 0; iFace < requestedKeys[iProc].size(); iFace++)
            {
                edge_[edgemap_[requestedKeys[iProc][iFace]]].globalIndex = packageEdgeGlobalIndicesRecv[iData++];
            }
           }
    }
#endif


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGhostElements()
     * ***************************************************************************/

#if HAVE_MPI
    // MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
    void ghostDistributeInterpolationOrders(
            std::vector< std::vector<int> > & thisProcessGhostElementIndices,
            std::vector< std::vector<int> > & neighborProcessGhostOrder
    )
    {
        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvbuf, recvcounts, rdispls;
        recvbuf.reserve(processBoundaryFaceIndex_.size());

        for (int i = 0; i < thisProcessGhostElementIndices.size(); i++)
        {
            for (int j = 0; j < thisProcessGhostElementIndices[i].size(); j++) { sendbuf.push_back(element_[thisProcessGhostElementIndices[i][j]].interpOrder); }
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
    // TODO: If planning to use with non-tetrahedral meshes, need to pass element type as well
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

        Dune::GeometryType meshGeometryType;
        meshGeometryType.makeTetrahedron();

        for (int i = 0; i < thisProcessGhostElementIndices.size(); i++)
        {
            int thisSendCounts = 0;
            int thisRecvCounts = 0;

            for (int j = 0; j < thisProcessGhostElementIndices[i].size(); j++)
            {
                int ghostElementIndex = thisProcessGhostElementIndices[i][j];
                int ghostElementFaceGlobalIndex = thisProcessGhostFaceGlobalIndices[i][j];
                int thisDofNum = element_[ghostElementIndex].vertexIndexSet.size();
                thisSendCounts += 3 + thisDofNum;
                thisRecvCounts += 3 + Dune::CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, neighborProcessGhostOrder[i][j]);
                totalRecvSize += thisRecvCounts;

                sendbuf.push_back(element_[ghostElementIndex].globalIndex);
                sendbuf.push_back(element_[ghostElementIndex].physicalTag);
                sendbuf.push_back(ghostElementFaceGlobalIndex);
                for (int iDof = 0; iDof < thisDofNum; iDof++)
                {
                    sendbuf.push_back(point_[element_[thisProcessGhostElementIndices[i][j]].vertexIndexSet[iDof]].globalIndex);
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
        Dune::GeometryType meshGeometryType;
        meshGeometryType.makeTetrahedron();

        for (int i = 0; i < neighborProcessGhostOrder.size(); i++) {

            for (int j = 0; j < neighborProcessGhostOrder[i].size(); j++)
            {
                EntityStorage thisElement;
                thisElement.geometryType = meshGeometryType;
                thisElement.globalIndex = packageGhostElementData[iData++];
                thisElement.interpOrder = neighborProcessGhostOrder[i][j];
                thisElement.physicalTag = packageGhostElementData[iData++];

                int associatedFaceGlobalIndex = packageGhostElementData[iData++];

                int thisElementDof = Dune::CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, thisElement.interpOrder);
                for (int iDof = 0; iDof < thisElementDof; iDof++)
                {
                    int thisVertexGlobalIndex = packageGhostElementData[iData++];
                    Index2IndexMap::iterator vertexIter = vertexGlobal2LocalMap_.find(thisVertexGlobalIndex);

                    // If this vertex already exists, just reuse it. Otherwise, create new vertex, and later request its coordinate
                    if (vertexIter != vertexGlobal2LocalMap_.end()) { thisElement.vertexIndexSet.push_back((*vertexIter).second); }
                    else
                    {
                        thisElement.vertexIndexSet.push_back(point_.size());

                        Vertex fakeCoord;
                        insertVertex(fakeCoord, thisVertexGlobalIndex);

                        // Note that this vertex needs communicating
                        missingVertices[i].insert(thisVertexGlobalIndex);
                    }
                }

                // Associate a (hopefully processBoundary) face with this ghost element
                face_[faceGlobal2LocalMap_[associatedFaceGlobalIndex]].element2Index = ghostElement_.size();
                ghostElement_.push_back(thisElement);
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
        std::vector<int> sendbuf;
        recvbuf.reserve(size_);

        // 4.1) MPI_alltoallv - tell each process the number of coordinates you want from it
        for (int i = 0; i < missingVertices.size(); i++)  { sendbuf.push_back(missingVertices[i].size()); }

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoall(sendbuf.data(), 1, MPI_INT, reinterpret_cast<int*>(recvbuf.data()), 1, MPI_INT, comm);


        // 4.2) MPI_alltoallv - tell each process the list of globalIds of coordinates you want from it
        // Cleanup
        std::vector<int> sendcounts, sdispls, recvcounts, rdispls;
        sendbuf.clear();
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

        for (int i = 0; i < size_; i++)
        {
            // Go through all vertices requested from this process. Package coordinates
            for (int j = 0; j < recvcounts[i]; j++)
            {
                int thisVertexGlobalIndex = packageMissingVertexGlobalIndices[iData++];
                int thisVertexLocalIndex = vertexGlobal2LocalMap_[thisVertexGlobalIndex];
                Vertex p = point_[thisVertexLocalIndex].coord;

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
                point_[vertexGlobal2LocalMap_[thisVertexGlobalIndex]].coord = thisCoord;
            }
        }
    }
#endif



    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    // Checks if given point fits into the bounding box of a Tetrahedron
    // FIXME: This method should be superseded by COM-CURVATURE-BOUND from CurvilinearGeometry
    /*
    bool is_inside_bounding_box_gracious(const Vertex & point) const {
        const double grace_tolerance = 1e-13;
        Vertex node_center, node_extent;

        meshBoundingBox(node_center, node_extent);
        node_center -= point;
        node_center.x = fabs(node_center.x);
        node_center.y = fabs(node_center.y);
        node_center.z = fabs(node_center.z);
        return node_center <= node_extent * (1.0 + grace_tolerance);
    }
    */

    // Gets a box in which this Tetrahedron fits
    // FIXME: This method should be superseded by COM-CURVATURE-BOUND from CurvilinearGeometry
    /*
    void elementBoundingBox(Vertex center, Vertex extent) const {
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
    */

    // Checks if a point is inside the element
    // FIXME: This method should be superseded by isInside method from CurvilinearGeometry
    /*
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
    */

    /** Find all tets in mesh which contain the specified point p. The found
        tets are returned in "tets". If p is outside the mesh an empty vector
        is returned. If p is on a shared face/edge/point, "tets" may contain more
        than one entry. */

    // TODO: PBE file allows femaxx to count time. Use alternative in Dune?
    /*
    void find_tets_by_point(const Vertex p, std::vector<int>& tets) const {
        assert(octree_);

        //pbe_start(132, "Octree traversal");

        // Get list of nodes (Tets) whose bounding box contain point p
        std::vector<int> nodes;
        int nof_visited = 0;
        octree_->find_nodes_by_point(p, nodes, nof_visited, &is_inside_bounding_box_gracious);

        // Create list of tets containing point p
        tets.resize(0);
        for (int i = 0; i < nodes.size(); i++) {
            if (is_inside_gracious(nodes[i], p)) { tets.push_back(nodes[i]); }
        }

        // pbe_stop(132);
        // rDebug("find_tets_by_point: nof_found=%d, nof_visited=%d", static_cast<int>(nodes.size()), nof_visited);
    }
    */


private: // Private members

    bool verbose_;
    bool withGhostElements_;


    // Storage necessary for user access and computation of globalIndex
    int nVertexTotal_;
    int nEdgeTotal_;
    int nFaceTotal_;
    int nElementTotal_;

    // Stores vertex coordinates and globalIds
    std::vector<VertexStorage> point_;

    // Stores all necessary data about edges, faces and elements
    // Minimalism - edges and faces do not store interpolatory vertex ids, but refer to an element subentity
    std::vector<EntityStorage> element_;
    std::vector<EntityStorage> ghostElement_;
    std::vector<EdgeStorage> edge_;
    std::vector<FaceStorage> face_;

    // Maps from global to local indices
    Index2IndexMap vertexGlobal2LocalMap_;
    Index2IndexMap edgeGlobal2LocalMap_;
    Index2IndexMap faceGlobal2LocalMap_;
    Index2IndexMap elementGlobal2LocalMap_;

    // List of localIndices of all faces of different structural types. Speeds up iterators
    std::vector<int> internalFaceIndex_;         // (self -> face_ index)
    std::vector<int> domainBoundaryFaceIndex_;   // (self -> face_ index)
    std::vector<int> processBoundaryFaceIndex_;  // (self -> face_ index)

    // List of all ranks of processors neighboring processorBoundaries. Index the same as processBoundaryFaceIndex_.
    std::vector<int> processBoundaryNeighborProcess_;  // (processBoundaryFaceIndex_ -> neighbor rank)

    // Temporary maps necessary to
    EdgeKey2EdgeIdMap edgemap_;                // (global edgeKey -> edge_ index)
    FaceKey2FaceIdMap internalInternalMap_;    // (global faceKey -> internalFaceIndex_)
    FaceKey2FaceIdMap boundaryInternalMap_;    // (global faceKey -> domainBoundaryFaceIndex_)
    FaceKey2FaceIdMap processInternalMap_;     // (global faceKey -> processBoundaryFaceIndex_)


    // Octree used to efficiently locate elements in which the points are located
    //LooseOctree<Tet>* octree_;
    int octree_;
    
    /** Pointer to single CurvilinearGridBase instance (Singleton) */
    //static CurvilinearGridBase* instance_ = 0;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDBASE_HH
