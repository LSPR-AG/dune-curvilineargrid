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
      int structuralType;
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

    typedef Dune::CurvilinearGridBase<ct>          GridBaseType;

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


    typedef Dune::CurvilinearOctreeNode<ct>                       NodeType;
    typedef Dune::CurvilinearLooseOctree<ct, 3, NodeType>         CurvilinearLooseOctree;

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridBase(bool withGhostElements, bool verbose, bool processVerbose, MPIHelper &mpihelper ) :
        nVertexTotal_(0),
        nEdgeTotal_(0),
        nFaceTotal_(0),
        nElementTotal_(0),
        withGhostElements_(withGhostElements),
        verbose_(verbose),
        processVerbose_(processVerbose),
        mpihelper_(mpihelper),
        octree_(0)
    {
        //assert(instance_ == 0);
        //instance_ = this;

        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        std::string log_string = "Initialized CurvilinearGridBase withGhostElements=" + std::to_string(withGhostElements);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
    }

    /** Destructor */
    ~CurvilinearGridBase() {
        //instance_ = 0;
    	delete octree_;
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

        std::string log_string = "CurvilinearGridBase: Inserted vertex LocalIndex=" + std::to_string(point_.size()-1) + " GlobalIndex=" + std::to_string(globalIndex);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
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

        std::string log_string = "CurvilinearGridBase: Inserted Element Type=" + Dune::CurvilinearGeometryHelper::geometryName(gt);
        log_string += " LocalIndex=" + std::to_string(element_.size()-1);
        log_string += " Order=" + std::to_string(order);
        log_string += " PhysicalTag=" + std::to_string(physicalTag);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
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
                thisFaceAsSubentity.structuralType = CurvilinearGridFaceType::DomainBoundary;
                thisFaceAsSubentity.element1Index = associatedElementIndex;
                thisFaceAsSubentity.element2Index = -1;              // Boundary Segments do not have a 2nd neighbor
                thisFaceAsSubentity.element1SubentityIndex = j;
                thisFaceAsSubentity.physicalTag = physicalTag;    // Here physical tag is very important as it need not match the tag of the element

                face_.push_back(thisFaceAsSubentity);


                std::string log_string = "CurvilinearGridBase: Inserted BoundarySegment Type=" + Dune::CurvilinearGeometryHelper::geometryName(gt);
                log_string += " LocalIndex=" + std::to_string(face_.size()-1);
                log_string += " Order=" + std::to_string(order);
                log_string += " PhysicalTag=" + std::to_string(physicalTag);
                log_string += " AssociatedElementIndex=" + std::to_string(associatedElementIndex);
                log_string += " InternalSubentityIndex=" + std::to_string(j);
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
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

    void generateMesh(int nVertexTotalMesh, int nElementTotalMesh) {

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Initializing mesh");


        // Temporary maps for Global Index construction
        // ************************************************************
        typedef FaceKey2FaceIdMap::iterator IterType;

        Index2IndexMap processBoundaryCornerMap;                            // (vertex global index -> processBoundaryCornerNeighborRank_ index)
        EdgeKey2EdgeIdMap processBoundaryEdgeMap;                           // (EdgeKey             -> processBoundaryEdgeNeighborRank_ index)


        // Construct missing parts of the mesh
        // ************************************************************

        nVertexTotal_ = nVertexTotalMesh;
        nElementTotal_ = nElementTotalMesh;

        computeProcessBoundingBox();
        generateEdges();
        generateFaces();

        if (size_ > 1)
        {
#if HAVE_MPI
        	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Assembling Parallel Grid");

        	// Parallel case
            generateProcessBoundaryCorners(processBoundaryCornerMap);
            generateProcessBoundaryEdges(processBoundaryEdgeMap);

            generateGlobalIndices(processBoundaryCornerMap, processBoundaryEdgeMap);
            if (withGhostElements_) { generateGhostElements(); }
#endif
        }
        else
        {
        	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Assembling Serial Grid");
            // Serial case:
            // * Boundary Neighbors not necessary, since all boundaries are domain boundaries
            // * No ghost elements, even if requested by user
            // * Fake globalIndex by making it equal to localIndex

            withGhostElements_ = false;
            for (int i = 0; i < edge_.size(); i++)     { edge_[i].globalIndex = i;     edgeGlobal2LocalMap_[i] = i; }
            for (int i = 0; i < face_.size(); i++)     { face_[i].globalIndex = i;     faceGlobal2LocalMap_[i] = i; }
            for (int i = 0; i < element_.size(); i++)  { element_[i].globalIndex = i;  elementGlobal2LocalMap_[i] = i; }

            nEdgeTotal_ = nEdge();
            nFaceTotal_ = nFace();
        }

        constructOctree();


        // Deletes all temporary memory
        // ************************************************************
        edgemap_.clear();
        internalInternalMap_.clear();
        boundaryInternalMap_.clear();
        processInternalMap_.clear();


        int log_tmp;
        std::string log_string = "CurvilinearGridBase: Constructed Mesh ";
        log_tmp = nVertexTotal();            log_string += " nVertexPerMesh="       + std::to_string(log_tmp);
        log_tmp = nEdgeTotal();              log_string += " nEdgePerMesh="         + std::to_string(log_tmp);
        log_tmp = nFaceTotal();              log_string += " nFacePerMesh="         + std::to_string(log_tmp);
        log_tmp = nElementTotal();           log_string += " nElementPerMesh="      + std::to_string(log_tmp);
        log_tmp = nVertex();                 log_string += " nVertex="              + std::to_string(log_tmp);
        log_tmp = nEdge();                   log_string += " nEdge="                + std::to_string(log_tmp);
        log_tmp = nFace();                   log_string += " nFace="                + std::to_string(log_tmp);
        log_tmp = nFaceInternal();           log_string += " nFaceInternal="        + std::to_string(log_tmp);
        log_tmp = nFaceProcessBoundary();    log_string += " nFaceProcessBoundary=" + std::to_string(log_tmp);
        log_tmp = nFaceDomainBoundary();     log_string += " nFaceDomainBoundary="  + std::to_string(log_tmp);
        log_tmp = nElement();                log_string += " nElement="             + std::to_string(log_tmp);
        log_tmp = nGhost();                  log_string += " nGhostElement="        + std::to_string(log_tmp);

        // Diagnostics output
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
    }

    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    /** Get total number of entities in a mesh  */
    int nVertexTotal() const   { return nVertexTotal_; }
    int nEdgeTotal() const     { return nEdgeTotal_; }
    int nFaceTotal() const     { return nFaceTotal_; }
    int nElementTotal() const  { return nElementTotal_; }

    /** Get total number of entities on this process  */
    int nVertex() const   { return point_.size(); }
    int nEdge() const     { return edge_.size(); }
    int nFace() const     { return face_.size(); }
    int nElement() const  { return element_.size(); }
    int nGhost() const    { return ghostElement_.size(); }

    int nFaceInternal() const        { return internalFaceIndex_.size(); }
    int nFaceProcessBoundary() const { return processBoundaryFaceIndex_.size(); }
    int nFaceDomainBoundary() const  { return domainBoundaryFaceIndex_.size(); }


    /** Vertex coordinate
     *  \param[in] localIndex            local vertex index (insertion index)
     * */
    Vertex vertex(int localIndex) const { return point_[localIndex].coord; }


    int facePhysicalTag(int localIndex) const { return face_[localIndex].physicalTag;}

    int elementPhysicalTag(int localIndex) const { return element_[localIndex].physicalTag;}

    int ghostElementPhysicalTag(int localIndex) const { return ghostElement_[localIndex].physicalTag;}

    int faceStructuralType(int localIndex) const { return face_[localIndex].structuralType; }

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
        case 0 : rez = face_[localIndex].element1Index;  break;
        case 1 : rez = face_[localIndex].element2Index;  break;
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
        EntityStorage & assocElement = element_[edge_[localIndex].elementIndex];
        EntityStorage thisEdge;
        thisEdge.geometryType.makeLine();
        thisEdge.globalIndex = edge_[localIndex].globalIndex;
        thisEdge.interpOrder = assocElement.interpOrder;
        thisEdge.physicalTag = -1;        // Note: Edges do not have a physical tag

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<int> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, 3>(assocElement.geometryType, thisEdge.interpOrder, 2, edge_[localIndex].subentityIndex );

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
        EntityStorage assocElement = element_[face_[localIndex].element1Index];
        EntityStorage thisFace;
        thisFace.geometryType.makeTriangle();
        thisFace.globalIndex = face_[localIndex].globalIndex;
        thisFace.interpOrder = assocElement.interpOrder;
        thisFace.physicalTag = face_[localIndex].physicalTag;

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<int> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, 3>(assocElement.geometryType, thisFace.interpOrder, 1, face_[localIndex].element1SubentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisFace;
    }


    /** Storage data related to this element, except of explicit vertex coordinates
     *  \param[in] localIndex            local element index
     *
     * */
    EntityStorage elementData(int localIndex) const { return element_[localIndex]; }


    /** Storage data related to this element, except of explicit vertex coordinates
     *  \param[in] localIndex            local ghost element index
     *
     * */
    EntityStorage ghostElementData(int localIndex) const { return ghostElement_[localIndex];}

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
    const CurvilinearLooseOctree & octree() const { return &octree_; }

    void processBoundingBox(Vertex & center, Vertex & extent) const
    {
        center = boundingBoxCenter_;
        extent = boundingBoxExtent_;
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
            octree_->findNode(globalC, elementIndices, nNodeVisited, &isInsideProcessBoundingBoxGracious);



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

    // Get curved geometry of an entity
    // TODO: assert mydim == element geometry type dim
    template<int mydim>
    Dune::CurvilinearGeometry<ct, mydim, 3> entityGeometry(
            Dune::GeometryType gt,
            std::vector<int> & vertexIndexSet,
            int order) const
    {
        VertexVector entityVertices;
        for (int i = 0; i < vertexIndexSet.size(); i++) { entityVertices.push_back(point_[vertexIndexSet[i]].coord); }

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

    // Takes two sorted arrays with non-repeating entries
    // Returns an array which only has entries found in both input arrays
    std::vector<int> sortedSetIntersection(std::vector<int> A, std::vector<int> B)
    {
        std::vector<int>  rez;

        int indA = 0;
        int indB = 0;

        while ((indA < A.size()) && (indB < B.size()))
        {
                  if (A[indA] < B[indB]) { indA++; }
             else if (A[indA] > B[indB]) { indB++; }
             else {
                 rez.push_back(A[indA]);
                 indA++;
                 indB++;
             }
        }

        return rez;
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
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started constructing edges");

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


                std::cout << "process_" << rank_ << "GenerateEdgeKey=(" << thisKey.node0 << "," << thisKey.node1 << ")" << std::endl;

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


                    std::string log_string = "CurvilinearGridBase: Added Edge";
                    log_string += " LocalIndex=" + std::to_string(edge_.size());
                    log_string += " AssociatedElementIndex=" + std::to_string(i);
                    log_string += " InternalSubentityIndex=" + std::to_string(j);
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                    edge_.push_back(thisEdge);
                }
            }
        }
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished constructing edges");
    }

    // Generates Internal and ProcessBoundary Faces. (!!!) Assumes that all Domain Boundary Faces have been added.
    // FIXME: Currently subentity orientation does not match the one of Dune
    // TODO: If assume mesh with non-uniform p-refinement, it may make sense to point the face to the element which has the higher refinement
    void generateFaces()
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started generating faces");

        typedef std::map<FaceKey, std::vector<int>> tmpFace2InfoMap;
        tmpFace2InfoMap tmpFaceMap;

        // Loop over all elements and their faces
        for (int i = 0; i < element_.size(); i++)
        {
            int nFacePerTetrahedron = 4;
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

                    if (iter != tmpFaceMap.end()) { connectedFaceInfo = std::vector<int> ( (*iter).second ); }

                    connectedFaceInfo.push_back(i);
                    connectedFaceInfo.push_back(j);

                    tmpFaceMap[thisKey] = connectedFaceInfo;


                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridBase: Adding FaceKey=(" << thisKey.node0 << ", " << thisKey.node1 << ", " << thisKey.node2 << ") attached to total of " << connectedFaceInfo.size() / 2 << " elements";
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
                }


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


            std::stringstream log_stream;
            log_stream << "CurvilinearGridBase: Added Face";
            log_stream << " LocalIndex=" << face_.size();
            log_stream << " AssociatedElementIndex=" << thisFace.element1Index;
            log_stream << " InternalSubentityIndex=" << thisFace.element1SubentityIndex;



            // Store internal, domain and process boundaries separately for faster iterators
            if (connectedFaceInfo.size() == 2)
            {
                thisFace.structuralType = CurvilinearGridFaceType::ProcessBoundary;

                // Store Map (key -> processBoundaryFaceIndex), Vector (processBoundaryFaceIndex -> faceIndex)
                processInternalMap_[(*iter).first] = processBoundaryFaceIndex_.size();
                processBoundaryFaceIndex_.push_back(face_.size());

                thisFace.element2Index = 0;    // Eventually this will be the Ghost Element Index

                log_stream << " StructuralType=processBoundary";
            }
            else
            {
                thisFace.structuralType = CurvilinearGridFaceType::Internal;

                // Store Map (key -> internalFaceIndex), Vector (internalFaceIndex -> faceIndex)
                internalInternalMap_[(*iter).first] = internalFaceIndex_.size();
                internalFaceIndex_.push_back(face_.size());

                thisFace.element2Index = connectedFaceInfo[2];  // This is the 2nd neighbor of this internal face

                log_stream << " StructuralType=internal";
            }


            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
            face_.push_back(thisFace);
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished generating faces");
    }


    // TODO: To make work with arbitrary geometries, must replace numbers, as well as make FaceKeys more abstract
    void generateProcessBoundaryCorners(Index2IndexMap & processBoundaryCornerMap)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started generating BoundaryCorneers");

        // Construct the set of EdgeKeys corresponding to edges of processBoundaries
        // ********************************************************
        for (FaceKey2FaceIdMap::iterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++)
        {
            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;

            int thisVertexKey[3] = {thisFaceKey.node0, thisFaceKey.node1, thisFaceKey.node2};

            for (int i = 0; i < 3; i++)
            {
                // For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (processBoundaryCornerMap.find(thisVertexKey[i]) == processBoundaryCornerMap.end())
                {
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    int processBoundaryCornerIndex = processBoundaryCornerMap.size();
                	processBoundaryCornerMap[thisVertexKey[i]] = processBoundaryCornerIndex;

                    std::string log_string = "CurvilinearGridBase: -- Adding boundary corner GlobalIndex=" + std::to_string(thisVertexKey[i]);
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
                }
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished generating BoundaryCorneers");
    }


    void generateProcessBoundaryEdges(EdgeKey2EdgeIdMap & processBoundaryEdgeMap)
    {
        // Construct the set of process boundary corners - corners necessary to make process boundary faces on this process
        // ********************************************************
        for (FaceKey2FaceIdMap::iterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++)
        {
            // Get global indices of the associated vertices from the map
            FaceKey thisFaceKey = (*faceIter).first;
            EdgeKey thisEdgeKey[3];

            thisEdgeKey[0].node0 = thisFaceKey.node0;  thisEdgeKey[0].node1 = thisFaceKey.node1;
            thisEdgeKey[1].node0 = thisFaceKey.node0;  thisEdgeKey[1].node1 = thisFaceKey.node2;
            thisEdgeKey[2].node0 = thisFaceKey.node1;  thisEdgeKey[2].node1 = thisFaceKey.node2;

            for (int i = 0; i < 3; i++)
            {
                // For each vertex, if it has not yet been added to the map, add it, mapping to a new entry
                // in the processBoundaryCornerNeighborRank_
                if (processBoundaryEdgeMap.find(thisEdgeKey[i]) == processBoundaryEdgeMap.end())
                {
                	// DO NOT combine these two lines into one - then the map returns the size after insertion which is +1 to expected
                    int processBoundaryEdgeIndex = processBoundaryEdgeMap.size();
                    processBoundaryEdgeMap[thisEdgeKey[i]]  = processBoundaryEdgeIndex;

                    std::stringstream log_stream;
                    log_stream << "CurvilinearGridBase: -- Adding boundary EdgeKey= (" << thisEdgeKey[i].node0 << ", " << thisEdgeKey[i].node1 << ")";
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
                }
            }
        }
    }



    // TODO: Inefficient Algorithm
    // The convention of owning an entity based on rank priority implies that processes with lower rank
    // have to do most of the work, and then perform a lot of communication. To balance out the workload
    // one would derive a more balanced owning paradigm
    //
    // Propose balanced paradigm:
    // Sum(Key[i]) mod nNeighbors
    void generateGlobalIndices(
        Index2IndexMap & processBoundaryCornerMap,
        EdgeKey2EdgeIdMap & processBoundaryEdgeMap
    )
    {
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Constructing Global Indices");


        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // 1) Communicate process ranks associated with each process boundary corner
        // Then compute that for edges and faces using set intersection
        // *************************************************************************
        std::vector<std::vector<int> > processBoundaryCornerNeighborRank( processBoundaryCornerMap.size(), std::vector<int>() );   // List of ranks of all other processes sharing this vertex
        std::vector<std::vector<int> > processBoundaryEdgeNeighborRank( processBoundaryEdgeMap.size(), std::vector<int>() );       // List of ranks of all other processes sharing this edge
        processBoundaryNeighborProcess_.resize(processInternalMap_.size());

        globalCommunicateVertexNeighborRanks(processBoundaryCornerMap, processBoundaryCornerNeighborRank);
        globalComputeEdgeNeighborRanks(processBoundaryCornerMap, processBoundaryEdgeMap, processBoundaryCornerNeighborRank, processBoundaryEdgeNeighborRank);
        globalComputeFaceNeighborRank(processBoundaryCornerMap, processBoundaryCornerNeighborRank);

        processBoundaryCornerMap.clear();
        processBoundaryCornerNeighborRank.clear();

        std::cout << "1" << std::endl;

        // 2) Get edges and faces on this process that are not owned by this process
        // *************************************************************************
        EdgeKey2EdgeIdMap edgeNonOwned;
        FaceKey2FaceIdMap faceNonOwned;

        std::cout << "2" << std::endl;

        // Edges
        for (EdgeKey2EdgeIdMap::iterator edgeIter = processBoundaryEdgeMap.begin(); edgeIter != processBoundaryEdgeMap.end(); edgeIter++ )
        {
            int edgeOwnerCandidateRank = processBoundaryEdgeNeighborRank[(*edgeIter).second][0];
            if (edgeOwnerCandidateRank < rank_) { edgeNonOwned[(*edgeIter).first] = edgeOwnerCandidateRank; }
        }

        std::cout << "3" << std::endl;

        // Faces
        for (FaceKey2FaceIdMap::iterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++ )
        {
            int faceOwnerCandidateRank = processBoundaryNeighborProcess_[(*faceIter).second];
            if (faceOwnerCandidateRank < rank_) { faceNonOwned[(*faceIter).first] = faceOwnerCandidateRank; }
        }

        int nEdgeOwned = edge_.size() - edgeNonOwned.size();
        int nFaceOwned = face_.size() - faceNonOwned.size();
        int elementsOwned = element_.size();


        std::cout << "4" << std::endl;

        // 3) Communicate number of edges and faces owned by each process to all
        // *************************************************************************
        std::vector<int> edgesOnProcess(size_);      // owned edges [rank]
        std::vector<int> facesOnProcess(size_);      // owned faces [rank]
        std::vector<int> elementsOnProcess(size_);  // owned elements [rank]

        collective_comm.allgather (&nEdgeOwned, 1, reinterpret_cast<int*> (edgesOnProcess.data()));
        collective_comm.allgather (&nFaceOwned, 1, reinterpret_cast<int*> (facesOnProcess.data()));
        collective_comm.allgather (&elementsOwned, 1, reinterpret_cast<int*> (elementsOnProcess.data()));

        std::cout << "5" << std::endl;

        int edgesBeforeMe = 0;      // Sum(edgesOwned : rank < thisRank)
        int facesBeforeMe = 0;      // Sum(facesOwned : rank < thisRank)
        int elementsBeforeMe = 0;   // Sum(elementsOwned : rank < thisRank)

        for (int iProc = 0; iProc < size_; iProc++)
        {
            nEdgeTotal_ += edgesOnProcess[iProc];
            nFaceTotal_ += facesOnProcess[iProc];

            if (iProc < rank_)
            {
                edgesBeforeMe += edgesOnProcess[iProc];
                facesBeforeMe += facesOnProcess[iProc];
                elementsBeforeMe += elementsOnProcess[iProc];
            }
        }

        std::cout << "6" << std::endl;

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
            if (faceNonOwned.find((*iter).first) == faceNonOwned.end())  { face_[(*iter).second].globalIndex = iFaceGlobalId++; }
        }

        for (EdgeKey2EdgeIdMap::iterator iter = edgemap_.begin(); iter != edgemap_.end(); iter++)
        {
            // This process owns this edge if it is in the edgemap and if it is not in the non-owned edge map
            if (edgeNonOwned.find((*iter).first) == edgeNonOwned.end())  { edge_[(*iter).second].globalIndex = iEdgeGlobalId++; }
        }

        std::cout << "7" << std::endl;


        // 4) Communicate missing edge and face global indices to their corresponding neighbors. Receive them and assign.
        // *************************************************************************

        globalDistributeMissingEdgeGlobalIndex(processBoundaryEdgeMap, processBoundaryEdgeNeighborRank);

        std::cout << "8.1" << std::endl;
        globalDistributeMissingFaceGlobalIndex();

        std::cout << "8.2" << std::endl;
    }


    /** Communicates the Ghost Elements
     *
     * Prerequisites:
     * * Requires all processBoundaries to know neighboring process rank
     * * Requires existence of face global index
     *
     * Aspects:
     * * Ghost elements can have different interpolation order
     *
     * Algorithm:
     * 1) Compute number of process boundaries shared with each process, as well as and the exact interpolation orders of Ghost Elements
     * 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the Ghost Elements it is going to receive
     * 3) MPI_alltoallv - Package and send element globalIndex + elementPhysicalTag + all interpolatory Vertex global indices
     * 4) Add all ghost elements to the mesh. Calculate which vertices are missing from which processes
     * 5) MPI_alltoall - Communicates to each process the number of missing vertices out of the ones it had communicated.
     * 5.1) Then communicate the globalIndices's of all missing vertices
     * 5.2) Then communicate the vertex coordinates corresponding to received global indices
     * 6) Distrubute vertex coordinates and add received coordinates to the mesh
     *
     *
     * TODO: Only supports tetrahedral ghost elements at the moment
     * [FIXME] Check if everything makes sense in terms of self-to-self communication
     * [FIXME] Check if there is balance between send and receive
     *
     * */
    void generateGhostElements()
    {

        // 1) Compute number of process boundaries shared with each process, as well as and the exact interpolation orders of Ghost Elements
        // *************************************************************************************

        // For each other process stores the set of element indices local to this process. These elements will be Ghost Elements on the other processes
        std::vector< std::vector<int> > thisProcessGhostElementIndexSet (size_, std::vector<int>() );
        // For each process store the set of face global indices, to which the corresponding Ghost Elements are associated
        std::vector< std::vector<int> > thisProcessGhostFaceGlobalIndexSet (size_, std::vector<int>() );

        for (int iPBFace = 0; iPBFace < processBoundaryFaceIndex_.size(); iPBFace++ )  {
            int thisFaceIndex = processBoundaryFaceIndex_[iPBFace];
            int thisNeighborRank = processBoundaryNeighborProcess_[iPBFace];

            int thisGhostIndex = face_[thisFaceIndex].element1Index;
            int thisFaceGlobalIndex = face_[thisFaceIndex].globalIndex;
            thisProcessGhostElementIndexSet[thisNeighborRank].push_back(thisGhostIndex);
            thisProcessGhostFaceGlobalIndexSet[thisNeighborRank].push_back(thisFaceGlobalIndex);
        }


        // 2) MPI_alltoallv - communicate to each process a list of interpolation orders of the ghost elements it is going to receive
        // *************************************************************************************

        // For each other process stores the set of interpolation orders of Ghost Elements that process wishes to communicate to this process
        std::vector< std::vector<int> > neighborProcessGhostOrder (size_, std::vector<int>() );
        ghostDistributeInterpolationOrders(thisProcessGhostElementIndexSet, neighborProcessGhostOrder);


        // 3) MPI_alltoallv - Package element globalIndex + elementPhysicalTag + all interpVertex globalIds
        // *************************************************************************************
        std::vector<int> packageGhostElementData;
        ghostDistributeGhostElements(thisProcessGhostElementIndexSet, thisProcessGhostFaceGlobalIndexSet, neighborProcessGhostOrder, packageGhostElementData );

        thisProcessGhostElementIndexSet.clear();
        thisProcessGhostFaceGlobalIndexSet.clear();


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


    /** Construct OCTree for locating tetrahedrons in mesh */

    // TODO: Use standard logging message
    // TODO: Original octree has diagnostics output under #if 0, can append at later stage

    // FIXME: Replace OCTree pointer by just an instance
    void constructOctree() {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started OCTree construction");

        // bounding box of whole mesh
        Vertex center, extent;
        processBoundingBox(center, extent);

        // octree length is the largest component of extent
        double length = extent[0];
        if (extent[1] > length)  { length = extent[1]; }
        if (extent[2] > length)  { length = extent[2]; }

        // construct LooseOctree with large max depth
        octree_ = new CurvilinearLooseOctree(center, length, 100, verbose_, processVerbose_, mpihelper_);

        // loop over all tets and insert them in the octree
        for (int iElem = 0; iElem < element_.size(); iElem ++)
        {
            const GridBaseType & thisGrid = *this;

            NodeType* thisNode = new NodeType(thisGrid, iElem);
            octree_->addNode(thisNode);
        }


        int maxDepth, nOctant, nNode;
        double avgNodeDepth;
        octree_->statistics(maxDepth, avgNodeDepth, nOctant, nNode);

        std::stringstream outputString;
        outputString << "CurvilinearGridBase: Constructed OCTree MaxDepth=" << maxDepth;
        outputString << ", #octants=" << nOctant;
        outputString << ", #nodes=" << nNode;
        outputString << ", avg. node depth=" << avgNodeDepth;
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, outputString.str());
    }


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGlobalIndices()
     * ***************************************************************************/

    /** Communicate the process ranks of neighbor processes for all process boundary vertices
     *
     * Algorithm:
     * 1) collective_comm.max() - find the maximal number of process boundary corners per process
     * 2) Loop over maximal number of process boundary corners per process
     * 2.1) collective_comm.allgather() - communicate a global index of your process boundary corner to all other processes
     * 2.2) If all process boundary corners of this process have already been communicated, create and communicate fake indices (negative)
     *      This is necessary to keep the protocol going until the last process communicates all its corners
     * 2.3) From received vertices, select ones that are on this process, and mark the sender as the neighbor
     *
     * [TODO] Algorithm possibly inefficient. No ideas how to improve at the moment
     * * Every process boundary corners is communicated to all processes, but can be used only by few
     * * All processes have to wait until the process with the largest number of corners finishes communicating, since they could receive sth from it
     *
     * */

    void globalCommunicateVertexNeighborRanks (
            Index2IndexMap & processBoundaryCornerMap,
            std::vector<std::vector<int> > & processBoundaryCornerNeighborRank
    )
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started communicating process boundary neighbors");

        // 1) collective_comm.max() - find the maximal number of process boundary corners per process
        // ********************************************************

        Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

        // Reserve memory for saving ranks associated to process boundary corners
        int thisProcessBoundarySize = processBoundaryCornerMap.size();

        int maxProcessBoundarySize = collective_comm.max(thisProcessBoundarySize);

        // 2) collective_comm.allgather() - communicate global index of your process boundary corner to all other processes
        // Repeat this process until every process has communicated all its corners
        // If you run out of corners, communicate fake corners
        // ********************************************************

        Index2IndexMap::iterator procCornerIter = processBoundaryCornerMap.begin();

        for (int iCorner = 0; iCorner < maxProcessBoundarySize; iCorner++)
        {
            // If all process boundary corners have been sent start sending fake corners
            int thisCornerInd = (iCorner < thisProcessBoundarySize) ? (*(procCornerIter++)).first : -1;

            // Communicate
            std::vector<int> procCornerIndexSet (size_);
            collective_comm.allgather(&thisCornerInd, 1, reinterpret_cast<int*> (procCornerIndexSet.data()) );

            // Loop over corners sent by other processes. If this corner present, note its sender rank
            for (int iProc = 0; iProc < size_; iProc++)
            {
                // Only consider non-fake corners sent by other processes
                if ((iProc != rank_) && (procCornerIndexSet[iProc] >= 0))
                {
                    // Attempt to find this corner global id among process boundary corners of this process
                    Index2IndexMap::iterator tmpIter = processBoundaryCornerMap.find(procCornerIndexSet[iProc]);

                    // If this corner is present, note its sender process
                    if (tmpIter != processBoundaryCornerMap.end()) {
                    	processBoundaryCornerNeighborRank[(*tmpIter).second].push_back(iProc);
                    }
                }
            }
        }

        // 3) Sort all neighbor rank sets, to accelerate set intersection algorithm in future
        // ********************************************************
        for (int i = 0; i < processBoundaryCornerNeighborRank.size(); i++)
        {
            std::sort(processBoundaryCornerNeighborRank[i].begin(), processBoundaryCornerNeighborRank[i].end());
        }


        // Testing output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridBase: -- Process boundary corner";
        for (Index2IndexMap::iterator cornerIter = processBoundaryCornerMap.begin(); cornerIter != processBoundaryCornerMap.end(); cornerIter++)
        {
        	log_stream << " GlobalIndex=" << (*cornerIter).first;
        	log_stream << " has Neighbors=(" << vector2string(processBoundaryCornerNeighborRank[(*cornerIter).second]) << ")";
        }
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished process boundary neighbors");
    }


    /** Compute the process ranks of neighbor processes for all process boundary edges
     *
     * Algorithm:
     * 1) Loop over all edges in the edge map
     * 1.1) For each corner in the EdgeKey get associated neighbor ranks from provided vertex neighbor ranks
     * 1.2) Perform intersection on the two sets
     * 1.3) Following edge map write that intersection to the output array
     *
     * */
    void globalComputeEdgeNeighborRanks(
            Index2IndexMap & processBoundaryCornerMap,
            EdgeKey2EdgeIdMap & processBoundaryEdgeMap,
            std::vector<std::vector<int> > & processBoundaryCornerNeighborRank,
            std::vector<std::vector<int> > & processBoundaryEdgeNeighborRank)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started computing edge process boundary neighbors");

        for (EdgeKey2EdgeIdMap::iterator edgeIter = processBoundaryEdgeMap.begin(); edgeIter != processBoundaryEdgeMap.end(); edgeIter++ )
        {
            // Get corners of the edge
            EdgeKey thisEdgeKey = (*edgeIter).first;

            // Get neighbor processes associated with each corner
            std::vector<int> corner0neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisEdgeKey.node0]];
            std::vector<int> corner1neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisEdgeKey.node1]];

            // Find neighbors common to both edge corners
            std::vector<int> edgeneighborset = sortedSetIntersection(corner0neighborset, corner1neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGrid: --";
            log_stream << "Neighbors[0]=(" << vector2string(corner0neighborset) << ")";
            log_stream << " Neighbors[1]=(" << vector2string(corner1neighborset) << ")";
            log_stream << " Intersection=" << vector2string(edgeneighborset);
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());

            if (edgeneighborset.size() < 1) { DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found no neighbor processes to an edge "); }

            // Store the edge neighbor rank set
            edgeneighborset.swap(processBoundaryEdgeNeighborRank[(*edgeIter).second]);
        }
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished computing edge process boundary neighbors");
    }


    /** Compute the process ranks of neighbor processes for all process boundary faces
     *
     *  Algorithm:
     *  1) Loop over all process boundary faces in the face map
     *  1.1) For each face corner in the FaceKey get associated neighbor ranks from provided vertex neighbor ranks
     *  1.2) Perform intersection on the three sets
     *  1.3) Ideally the intersection should result in one single rank, which is this face's neighbor. Otherwise throw error
     *
     * */
    void globalComputeFaceNeighborRank(
            Index2IndexMap & processBoundaryCornerMap,
            std::vector<std::vector<int> > & processBoundaryCornerNeighborRank)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started computing face process boundary neighbors");

        for (FaceKey2FaceIdMap::iterator faceIter = processInternalMap_.begin(); faceIter != processInternalMap_.end(); faceIter++ )
        {
            // Get corners of the edge
            FaceKey thisFaceKey = (*faceIter).first;

            // Get neighbor processes associated with each corner
            std::vector<int> corner0neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisFaceKey.node0]];
            std::vector<int> corner1neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisFaceKey.node1]];
            std::vector<int> corner2neighborset = processBoundaryCornerNeighborRank[processBoundaryCornerMap[thisFaceKey.node2]];

            // Find neighbors common to all 3 face corners. Need to intersect sets twice
            std::vector<int> faceneighborset;
            faceneighborset = sortedSetIntersection(corner0neighborset, corner1neighborset);
            faceneighborset = sortedSetIntersection(faceneighborset,    corner2neighborset);

            std::stringstream log_stream;
            log_stream << "CurvilinearGrid: --";
            log_stream << "Neighbors[0]=(" << vector2string(corner0neighborset) << ")";
            log_stream << " Neighbors[1]=(" << vector2string(corner1neighborset) << ")";
            log_stream << " Neighbors[2]=(" << vector2string(corner2neighborset) << ")";
            log_stream << " Intersection=" << vector2string(faceneighborset);
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());

            if (faceneighborset.size() != 1) { DUNE_THROW(Dune::IOError, "CurvilinearGrid: Found wrong number of neighbor processes to a face"); }

            // Store the face neighbor rank. Face is only allowed to have exactly one neighbor
            processBoundaryNeighborProcess_[(*faceIter).second] = faceneighborset[0];
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished computing face process boundary neighbors");
    }


    /** Communicates all process boundary face global Id's to the neighbors if owned
     *
     * Algorithm:
     *
     * 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
     * 1.1) Assemble a global index array to send to each process
     * 1.2) Note how many faces will be received from each process
     * 2) Assemble one big send array from small arrays (FaceKey + globalIndex)
     * 3) MPI_Alltoallv - communicate this array
     * 4) Save global indices for non-owned faces. Find the exact face by using the communicated FaceKey
     *
     * Optimization Proposal:
     * In principle communication of the FaceKey is not necessary. Instead, the natural FaceKey "<" operator
     * can be used to sort all communicated faces, thus allowing the receiving process to "figure out" what are
     * the faces sent to it by sorting its own faces.
     * 1) Sort all process boundary faces wrt FaceKey
     * 2) Fill arrays to send according to this sorted order
     * 3) Make map for each process from rank & received face to local face index accoridng to the sorted FaceKey order
     * 3) Communicate only globalIndices
     * 4) Use constructed map to fill in received global Indices
     *
     * Will decrease the global communication at the expense of increasing local computation time
     *
     * */
    void globalDistributeMissingFaceGlobalIndex()
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started distributing missing face GlobalIndices");

        typedef std::pair<FaceKey, int>  FaceInfo;
        std::vector< std::vector< FaceInfo > > facesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_FACEINFO = 4;
        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvbuf, recvcounts, rdispls;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Checking which faces are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        for (FaceKey2FaceIdMap::iterator iter = processInternalMap_.begin(); iter != processInternalMap_.end(); iter++)
        {
            int localProcessBoundaryIndex = (*iter).second;
            int localFaceIndex = processBoundaryFaceIndex_[localProcessBoundaryIndex];
            int neighborRank = processBoundaryNeighborProcess_[localProcessBoundaryIndex];

            // If the neighbor of this face has lower rank, then add it to recv, else to send
            if (neighborRank < rank_)  { recvcounts[neighborRank] += N_INTEGER_FACEINFO; }
            else
            {
                int thisGlobalIndex = face_[localFaceIndex].globalIndex;
                facesToSend[neighborRank].push_back(FaceInfo((*iter).first, thisGlobalIndex ));
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
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

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Sending");

        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Extracting missing indices");

        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
            int nThisFaceInfo = recvcounts[iProc] / N_INTEGER_FACEINFO;

            for (int iFace = 0; iFace < nThisFaceInfo; iFace++)
            {
                FaceKey thisKey;
                int thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];
                thisKey.node2 = recvbuf[iData++];

                FaceKey2FaceIdMap::iterator faceIter = processInternalMap_.find(thisKey);

                if (faceIter == processInternalMap_.end()) { DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated FaceKey does not correspond to any face on this process "); }
                else
                {
                    int localFaceIndex = processBoundaryFaceIndex_[(*faceIter).second];
                    face_[localFaceIndex].globalIndex = thisGlobalId;
                }
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished distributing missing face GlobalIndices");
    }


    /** Communicates all process boundary face global Id's to the neighbors if owned
     *
     * Algorithm:
     *
     * 1) Loop over all process boundary edges, split edges into ones to be sent and to be received
     * 1.1) If this edge rank lower than all other neighbor ranks, note to send it to all neighbors,
     *      Otherwise note which neighbor to receive it from
     * 1.1) Assemble a global index array to send to each process
     * 1.2) Note how many edges will be received from each process
     * 2) Assemble one big send array from small arrays (EdgeKey + globalIndex)
     * 3) MPI_Alltoallv - communicate this array
     * 4) Save global indices for non-owned edges. Find the exact edge by using the communicated EdgeKey
     *
     * Optimization Proposal:
     * Same as for FaceGlobalIndex algorithm
     *
     * */
    void globalDistributeMissingEdgeGlobalIndex(
            EdgeKey2EdgeIdMap & processBoundaryEdgeMap,
            std::vector<std::vector<int> > & processBoundaryEdgeNeighborRank)
    {
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Started distributing missing edge GlobalIndices");

        typedef std::pair<EdgeKey, int>  EdgeInfo;
        std::vector< std::vector< EdgeInfo > > edgesToSend (size_);

        int totalRecvSize = 0;
        int N_INTEGER_EDGEINFO = 3;
        std::vector<int> sendbuf, sendcounts(size_), sdispls;
        std::vector<int> recvbuf, recvcounts(size_), rdispls;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Checking which edges are missing");

        // 1) Loop over all process boundary faces, split faces on the ones to be sent and to be received
        // ********************************************************************************************
        for (EdgeKey2EdgeIdMap::iterator iter = processBoundaryEdgeMap.begin(); iter != processBoundaryEdgeMap.end(); iter++)
        {
            EdgeKey thisEdgeKey = (*iter).first;
            int localProcessBoundaryIndex = (*iter).second;

            EdgeKey2EdgeIdMap::iterator thisEdgeMapIter = edgemap_.find(thisEdgeKey);
            if (thisEdgeMapIter == edgemap_.end()) {
            	std::cout << ":((" << std::endl;
            	std::cout << "process_" << rank_ << "ThisEdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ")" << std::endl;
            	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Process boundary edge not found among all edges");
            }

            int localEdgeIndex = (*thisEdgeMapIter).second;
            int candidateOwnerRank = processBoundaryEdgeNeighborRank[localProcessBoundaryIndex][0];

            std::cout << "process_" << rank_ <<  " EdgeKey=(" << thisEdgeKey.node0 << "," << thisEdgeKey.node1 << ") localIndex=" << localEdgeIndex <<  std::endl;

            // If the one of the neighbors of this edge has lower rank, then note one more received edge from that process
            // else note to send it to all other neighbors
            if (candidateOwnerRank < rank_)  { recvcounts[candidateOwnerRank] += N_INTEGER_EDGEINFO; }
            else
            {
                int thisGlobalIndex = edge_[localEdgeIndex].globalIndex;

                EdgeInfo thisEdgeInfo(thisEdgeKey, thisGlobalIndex);

                for (int iNeighbor = 0; iNeighbor < processBoundaryEdgeNeighborRank[localProcessBoundaryIndex].size(); iNeighbor++)
                {
                    int thisNeighborRank = processBoundaryEdgeNeighborRank[localProcessBoundaryIndex][iNeighbor];
                    edgesToSend[thisNeighborRank].push_back(thisEdgeInfo);
                };
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Assembling arrays to send");


        // 2) Fill in communication arrays
        // ********************************************************************************************
        for (int i = 0; i < size_; i++)
        {
            sendcounts[edgesToSend[i].size()];
            totalRecvSize += recvcounts[i];
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );

            for (int j = 0; j < edgesToSend[i].size(); j++)
            {
                sendbuf.push_back(edgesToSend[i][j].second);
                sendbuf.push_back(edgesToSend[i][j].first.node0);
                sendbuf.push_back(edgesToSend[i][j].first.node1);
            }

        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Sending");


        std::cout << "process_" << rank_ <<  " sendbuf: " << vector2string(sendbuf) << " sendcounts: " << vector2string(sendcounts) << " recvcounts: " << vector2string(recvcounts) << std::endl;


        // 3) MPI_alltoall key + globalId
        // ********************************************************************************************
        recvbuf.resize(totalRecvSize, 0);   // There are 3 integers per sent FaceInfo
        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: -- Extracting missing indices");

        // 4) Mark all missing faces
        // ********************************************************************************************8
        // Note that recvcounts should be 0 for this rank, as it was assembled by only considering the neighbor ranks
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
            int nThisEdgeInfo = recvcounts[iProc] / N_INTEGER_EDGEINFO;

            for (int iEdge = 0; iEdge < nThisEdgeInfo; iEdge++)
            {
                EdgeKey thisKey;
                int thisGlobalId = recvbuf[iData++];
                thisKey.node0 = recvbuf[iData++];
                thisKey.node1 = recvbuf[iData++];

                EdgeKey2EdgeIdMap::iterator edgeIter = edgemap_.find(thisKey);

                if (edgeIter == edgemap_.end()) {
                	DUNE_THROW(Dune::IOError, "CurvilinearGrid: Communicated EdgeKey does not correspond to any edge on this process "); }
                else
                {
                    int localEdgeIndex = (*edgeIter).second;
                    edge_[localEdgeIndex].globalIndex = thisGlobalId;
                }
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Finished distributing missing edge GlobalIndices");
    }


    /* ***************************************************************************
     * SubSection: Auxiliary methods of generateGhostElements()
     * ***************************************************************************/


    /** Communicate to each process a list of interpolation orders of the ghost elements it is going to receive
     *
     * Algorithm:
     *
     * 1) Put interpolation orders of elements of this process that will become Ghost Elements on other processes
     *    into send array in the correct order
     * 2) Communicate
     * 3) Store interpolation orders of Ghost Elements that will be sent to this process from each other process
     *
     * Optimization Proposal:
     * Same as for FaceGlobalIndex algorithm
     *
     * */
    void ghostDistributeInterpolationOrders(
            const std::vector< std::vector<int> > & thisProcessGhostElementIndexSet,
            std::vector< std::vector<int> > & neighborProcessGhostOrder
    )
    {
        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvbuf, recvcounts, rdispls;

        // We should receive the number interpolation orders equal to the number of process boundaries
        recvbuf.resize(processBoundaryFaceIndex_.size(), 0);

        for (int iProc = 0; iProc < size_; iProc++)
        {
            for (int iElem = 0; iElem < thisProcessGhostElementIndexSet[iProc].size(); iElem++) {
                sendbuf.push_back(element_[thisProcessGhostElementIndexSet[iProc][iElem]].interpOrder);
            }
            sendcounts.push_back(thisProcessGhostElementIndexSet[iProc].size());
            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + sendcounts[iProc-1] );

            // For ghost elements we send/receive same amount to/from each processor
            recvcounts.push_back(sendcounts[iProc]);
            rdispls.push_back(sdispls[iProc]);
        }

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );


        // Parse the received data
        // Note that with Ghost Elements, we receive from each neighbor process as many entities as we send to it
        int iData = 0;
        for (int iProc = 0; iProc < size_; iProc++)
        {
            for (int iElem = 0; iElem < thisProcessGhostElementIndexSet[iProc].size(); iElem++)
            {
                neighborProcessGhostOrder[iProc].push_back(recvbuf[iData++]);
            }
        }
    }


    //


    /** Communicate to all ghost element information except of the explicit vertex coordinates
     * MPI_alltoallv - Package element globalIndex + elementPhysicalTag + assocFaceGlobalIndex + all interpVertex globalIds
     *
     * Optimization Proposal:
     * It is possible not to send the associated face global index, and instead figure out the order of elements
     * according to the order of associated FaceKeys. The gain however is very small
     *
     * TODO: If planning to use with non-tetrahedral meshes, need to pass element type as well
     *
     * */
    void ghostDistributeGhostElements(
            std::vector< std::vector<int> > & thisProcessGhostElementIndexSet,
            std::vector< std::vector<int> > & thisProcessGhostFaceGlobalIndexSet,
            std::vector< std::vector<int> > & neighborProcessGhostOrder,
            std::vector<int> & recvbuf )
    {
        std::vector<int> sendbuf, sendcounts, sdispls;
        std::vector<int> recvcounts, rdispls;

        // Calculates total amount of integers to receive during DoF communication stage
        int totalRecvSize = 0;

        Dune::GeometryType meshGeometryType;
        meshGeometryType.makeTetrahedron();

        for (int i = 0; i < thisProcessGhostElementIndexSet.size(); i++)
        {
            int thisSendCounts = 0;
            int thisRecvCounts = 0;

            for (int j = 0; j < thisProcessGhostElementIndexSet[i].size(); j++)
            {
                int ghostElementIndex = thisProcessGhostElementIndexSet[i][j];
                int ghostElementFaceGlobalIndex = thisProcessGhostFaceGlobalIndexSet[i][j];
                int thisDofNum = element_[ghostElementIndex].vertexIndexSet.size();
                thisSendCounts += 3 + thisDofNum;
                thisRecvCounts += 3 + Dune::CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, neighborProcessGhostOrder[i][j]);
                totalRecvSize += thisRecvCounts;

                sendbuf.push_back(element_[ghostElementIndex].globalIndex);
                sendbuf.push_back(element_[ghostElementIndex].physicalTag);
                sendbuf.push_back(ghostElementFaceGlobalIndex);
                for (int iDof = 0; iDof < thisDofNum; iDof++)
                {
                    sendbuf.push_back(point_[element_[thisProcessGhostElementIndexSet[i][j]].vertexIndexSet[iDof]].globalIndex);
                }
            }

            sendcounts.push_back(thisSendCounts);
            recvcounts.push_back(thisRecvCounts);
            sdispls.push_back((i == 0) ? 0 : sdispls[i-1] + sendcounts[i-1] );
            rdispls.push_back((i == 0) ? 0 : rdispls[i-1] + recvcounts[i-1] );
        }

        recvbuf.resize(totalRecvSize, 0);

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );
    }


    /** Add received elements to the mesh. For each vertex global index, find if coordinate is already present on this process
     *  If not, mark this vertex as a missing vertex for further communication.
     *
     *
     * */
    void ghostInsertGhostElements (
            std::vector< std::vector<int> > & neighborProcessGhostOrder,
            std::vector< int > & packageGhostElementData,
            std::vector<std::set<int> > & missingVertices
    )
    {
        int iData = 0;
        Dune::GeometryType meshGeometryType;
        meshGeometryType.makeTetrahedron();

        for (int iProc = 0; iProc < neighborProcessGhostOrder.size(); iProc++) {

            for (int iGhost = 0; iGhost < neighborProcessGhostOrder[iProc].size(); iGhost++)
            {
                EntityStorage thisElement;
                thisElement.geometryType = meshGeometryType;
                thisElement.globalIndex = packageGhostElementData[iData++];
                thisElement.interpOrder = neighborProcessGhostOrder[iProc][iGhost];
                thisElement.physicalTag = packageGhostElementData[iData++];

                int associatedFaceGlobalIndex = packageGhostElementData[iData++];

                int thisElementDof = Dune::CurvilinearGeometryHelper::dofPerOrder(meshGeometryType, thisElement.interpOrder);
                for (int iDof = 0; iDof < thisElementDof; iDof++)
                {
                    int thisVertexGlobalIndex = packageGhostElementData[iData++];
                    Index2IndexMap::iterator vertexIter = vertexGlobal2LocalMap_.find(thisVertexGlobalIndex);

                    // If this vertex already exists, just reuse it. Otherwise, create new vertex, and later request its coordinate
                    if (vertexIter != vertexGlobal2LocalMap_.end()) {
                        int thisVerexLocalIndex = (*vertexIter).second;
                        thisElement.vertexIndexSet.push_back(thisVerexLocalIndex);
                    }
                    else
                    {
                        // Create a new vertex with local index pointing to the end of current vertex array
                        thisElement.vertexIndexSet.push_back(point_.size());

                        // Insert the fake vertex into the mesh
                        Vertex fakeCoord;
                        insertVertex(fakeCoord, thisVertexGlobalIndex);

                        // Note that this vertex needs communicating
                        missingVertices[iProc].insert(thisVertexGlobalIndex);
                    }
                }

                // Associate a (hopefully processBoundary) face with this ghost element
                face_[faceGlobal2LocalMap_[associatedFaceGlobalIndex]].element2Index = ghostElement_.size();
                ghostElement_.push_back(thisElement);
            }
        }
    }


    //
    //

    /** Communicate to each process the number of missing vertices out of the ones it had provided with Ghost Elements
     *  Then communicate the globalIndices of all missing vertices
     *
     * */
    void ghostDistributeMissingVertexGlobalIndices(
            std::vector<std::set<int> > & missingVertices,
            std::vector<int> & recvbuf,
            std::vector<int> & verticesRequested,
            std::vector<int> & verticesToSend
    )
    {
        std::vector<int> sendbuf;
        recvbuf.resize(size_, 0);

        // 4.1) MPI_alltoallv - tell each process the number of coordinates you want from it
        for (int iProc = 0; iProc < size_; iProc++)  { sendbuf.push_back(missingVertices[iProc].size()); }

        MPI_Comm comm = Dune::MPIHelper::getCommunicator();
        MPI_Alltoall(sendbuf.data(), 1, MPI_INT, reinterpret_cast<int*>(recvbuf.data()), 1, MPI_INT, comm);


        // 4.2) MPI_alltoallv - tell each process the list of global Indices of coordinates you want from it
        // Cleanup
        std::vector<int> sendcounts, sdispls, recvcounts, rdispls;
        sendbuf.clear();
        int totalRecvSize = 0;
        int iData = 0;

        for (int iProc = 0; iProc < size_; iProc++)
        {
            int thisSendSize = missingVertices[iProc].size();
            int thisRecvSize = recvbuf[iData++];

            recvcounts.push_back(thisRecvSize);
            sendcounts.push_back(thisSendSize);
            totalRecvSize += thisRecvSize;

            for (std::set<int>::iterator mVertIter = missingVertices[iProc].begin(); mVertIter != missingVertices[iProc].end(); mVertIter++ )  { sendbuf.push_back(*mVertIter); }

            sdispls.push_back((iProc == 0) ? 0 : sdispls[iProc-1] + sendcounts[iProc-1] );
            rdispls.push_back((iProc == 0) ? 0 : rdispls[iProc-1] + recvcounts[iProc-1] );
        }

        recvbuf.clear(); recvbuf.resize(totalRecvSize, 0);
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

        recvbuf.resize(totalRecvSize, 0);

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



    /* ***************************************************************************
     * Section: Methods of the mesh
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
        Vertex min = point_[0].coord;
        Vertex max = min;

        for (int i = 1; i < point_.size(); i ++) {
            min[0] = std::min(min[0], point_[i].coord[0]);
            min[1] = std::min(min[1], point_[i].coord[1]);
            min[2] = std::min(min[2], point_[i].coord[2]);

            max[0] = std::max(max[0], point_[i].coord[0]);
            max[1] = std::max(max[1], point_[i].coord[1]);
            max[2] = std::max(max[2], point_[i].coord[2]);
        }
        boundingBoxCenter_ = min + max;  boundingBoxCenter_ *= 0.5;
        boundingBoxExtent_ = max - min;  boundingBoxExtent_ *= 0.5;
    }

    // Checks if given point fits into the bounding box of the mesh
    bool isInsideProcessBoundingBoxGracious(const Vertex & point) const {
        const double grace_tolerance = 1e-13;
        Vertex center, extent;

        bool isInside = true;

        isInside &= fabs(boundingBoxCenter_[0] - point[0]) <= boundingBoxExtent_[0] * (1.0 + grace_tolerance);
        isInside &= fabs(boundingBoxCenter_[1] - point[1]) <= boundingBoxExtent_[1] * (1.0 + grace_tolerance);
        isInside &= fabs(boundingBoxCenter_[2] - point[2]) <= boundingBoxExtent_[2] * (1.0 + grace_tolerance);

        return isInside;
    }

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::string tmp_str;
        for (int i = 0; i < V.size(); i++) { tmp_str += std::to_string(V[i]) + " "; }
        return tmp_str;
    }

private: // Private members

    bool verbose_;
    bool processVerbose_;
    bool withGhostElements_;


    // Storage of process bounding box, since its computation is expensive
    Vertex boundingBoxCenter_;
    Vertex boundingBoxExtent_;

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

    // Temporary maps necessary to locate and communicate entities during grid base construction
    EdgeKey2EdgeIdMap edgemap_;                // (global edgeKey -> edge_ index)
    FaceKey2FaceIdMap internalInternalMap_;    // (global faceKey -> internalFaceIndex_)
    FaceKey2FaceIdMap boundaryInternalMap_;    // (global faceKey -> domainBoundaryFaceIndex_)
    FaceKey2FaceIdMap processInternalMap_;     // (global faceKey -> processBoundaryFaceIndex_)

    // Octree used to efficiently locate elements in which the points are located
    CurvilinearLooseOctree * octree_;
    
    /** Pointer to single CurvilinearGridBase instance (Singleton) */
    //static CurvilinearGridBase* instance_ = 0;

    MPIHelper &mpihelper_;
    int rank_;
    int size_;
};

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDBASE_HH
