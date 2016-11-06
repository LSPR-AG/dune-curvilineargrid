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

#include <dune/grid/common/gridenums.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>

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
 *  - [TODO] Does NOT support front/overlap elements at the moment
 *  - [TODO] DOes NOT support periodic boundaries at the moment
 *  - [TODO] Does NOT support non-tetrahedral meshes. Generalization to arbitrary and mixed geometry meshes is possible but will be cumbersome
 *
 * Development log
 *  - [TODO] To make some member functions const, need to introduce const iterators
 *  - [TODO] Consider making GridConstructor a pointer and deleting it at the end of construction phase.
 *  - [TODO] Convert all code to unsigned ints wherever the number is definitely non-negative
 *  - Currently nEntity counts ghost elements. Is this expected [Yes]
 *  - Currently nEntity counts only corners. Is this expected [Yes]
 *
 * Usage:
 *  - [TODO] Disable all the vertex2string output for multiprocessor case - too much output
 *  - [TODO] Implement timing (and perhaps memory usage) for each of the construction operations
 *  - [TODO] Implement memory in LoggingMessage as is done in Hades
 *  - [TODO] Implement TimeSync logging message so processes do not comment on top of each other. Check Boost if already exists
 *
 *
 * Testing:
 *  - [TODO] Constructor run check if all non-owned entities have been successfully enumerated at the end
 *
 *
 * Additional Functionality - Extended Physical Tags
 *  - Solution 1: Implement Vector Physical Tags, read from GMSH. Involves changing tag reading in GMSH and tag communication in grid construction
 *  - Solution 2: Keep single integer tag. Involves implementing a tag->info mapper in the derived code.
 *
 * Additional Functionality - Internal Surfaces
 *  - [TODO] When reading GMSH surfaces, stop assuming they are boundary surfaces. Insert them as simple BoundarySegments
 *  - [TODO] When constructing faces, do not check for existence of Domain Boundaries.
 *  - [TODO] When finding face neighbor, allow 0 neighbor-intersect and mark as Domain Boundary
 *  - [TODO] When cross-checking fake faces and edges, may not assume that having 1 neighbor -> neighbor is correct. Have to cross check all faces except of those who have no neighbors at all
 *
 * Optimization of CurvilinearVTKWriter:
 *    - Problem: Output too many vertices, many vertices are used several times. However, want to keep functionality to insert
 *               elements completely unrelated to the mesh if necessary
 *  - [TODO] Insert element together with both coordinates and their global indices. If global index already mapped, compare coordinates and throw errors if difference is significant
 *  - [TODO] In writer all vertices have their own local index, however, there are 3 types of them
 *      - [TODO] Elements inserted with vertices without specifying their global indices, then these vertices only possess local index
 *      - [TODO] Elements inserted with specifying global indices, then local index  and  global2local map
 *      - [TODO] Discretization vertices. A further optimization stage is to store vertices based on their ParentKey+DiscretizationKey
 *              - Example for the key: Have 3 maps: 2point, 3point and 4point
 *              - The key is 2,3 or 4 vertices + 2,3 or 4 fractions of each point used for interpolation as seen from the sorted key orientation
 *              - So for each element, calculate fractions, if some of them are 0, then reduce to lower key. This way no vertex will be stored twice
 *
 *
 *  ***************************************************************************/




namespace Dune {

// *******************************************
// Forwards-declaration
// *******************************************
template<class ct, int cdim, bool isCached>
class CurvilinearGridStorage;



template <class ct, int cdim, bool isCached>
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

	typedef ct        ctype;
	static const int dimension = cdim;
	static const int dimensionworld = cdim;
	static const int is_cached = isCached;

    typedef Dune::CurvilinearGridBase<ct, cdim, isCached>               GridBaseType;
    typedef typename Dune::CurvilinearGridStorage<ct, cdim, isCached>   GridStorageType;
    typedef typename Dune::CurvilinearGridConstructor<GridBaseType>     GridConstructorType;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    typedef typename GridStorageType::Vertex                 Vertex;
    typedef typename GridStorageType::VertexStorage          VertexStorage;
    typedef typename GridStorageType::EdgeStorage            EdgeStorage;
    typedef typename GridStorageType::FaceStorage            FaceStorage;
    typedef typename GridStorageType::EntityStorage          EntityStorage;

    typedef typename GridStorageType::EdgeKey                EdgeKey;
    typedef typename GridStorageType::FaceKey                FaceKey;
    typedef typename GridStorageType::IdType                 IdType;

    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::Global2LocalConstIterator      Global2LocalConstIterator;
    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;

    typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;

    typedef typename GridStorageType::EntityNeighborRankVector  EntityNeighborRankVector;

    typedef typename Dune::LoggingMessage                       LoggingMessage;
    typedef typename Dune::LoggingTimer<LoggingMessage>         LoggingTimer;

    // Defines the Elementary geometry
    template <int codim>
    struct Codim
    {
    	typedef typename std::conditional<
    	  isCached,
    	  Dune::CachedCurvilinearGeometry<ct, cdim - codim, cdim>,
    	  Dune::CurvilinearGeometry<ct, cdim - codim, cdim>
    	>::type  EntityGeometry;
    };

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;
    static const unsigned int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;




public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridBase(bool withGhostElements, bool withElementGlobalIndex, MPIHelper &mpihelper) :
        gridstage_(0),
        mpihelper_(mpihelper),
        gridstorage_(withGhostElements, withElementGlobalIndex)
    {
    	gridconstructor_ = new GridConstructorType(gridstorage_, *this, mpihelper);

        rank_ = mpihelper_.rank();
        size_ = mpihelper_.size();

        std::string log_string = "Initialized CurvilinearGridBase withGhostElements=" + std::to_string(withGhostElements);
        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
    }

private:
    /** Copy constructor: private, undefined: disallow copy */
    CurvilinearGridBase(const CurvilinearGridBase&);

    /** Assignment operator: private, undefined: disallow assignment */
    //CurvilinearGridBase& operator=(const CurvilinearGridBase&);

public:

    /* ***************************************************************************
     * Section: Loading the mesh
     * ***************************************************************************/


    /** \brief Add a new vertex to the mesh
     * \param[in] globalIndex      global index of this vertex
     * \param[in] p                coordinate of this vertex
     * */
    void insertVertex(Vertex p, GlobalIndexType globalIndex)
    {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridconstructor_->insertVertex(p, globalIndex);
    }

    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * */
    void insertElement(
        	Dune::GeometryType gt,
        	const std::vector<LocalIndexType> & vertexIndexSet,
        	GlobalIndexType globalIndex,
        	InterpolatoryOrderType order,
        	PhysicalTagType physicalTag)
    {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridconstructor_->insertElement(gt, vertexIndexSet, globalIndex, order, physicalTag);
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
     *  \param[in] associatedElementIndex   local index of the element this face is associated to
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
    	assertStage(Stage::GRID_CONSTRUCTION);
    	// Note: associatedElementIndex no longer necessary
    	// gridconstructor_.insertBoundarySegment(gt, associatedElementIndex, vertexIndexSet, order, physicalTag, isDomainBoundary);
    	gridconstructor_->insertBoundarySegment(gt, vertexIndexSet, order, physicalTag, isDomainBoundary);
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

    void generateMesh() {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridstage_ = Stage::GRID_OPERATION;

        LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "[[CurvilinearGridBase: Generating curvilinear mesh...");

        gridconstructor_->generateMesh();

        // Free up the memory taken by the construction procedure
        delete(gridconstructor_);

        // Diagnostics output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridBase: Constructed Mesh ";
        log_stream << " nVertexPerMesh="             << nEntityTotal(VERTEX_CODIM);
        log_stream << " nEdgePerMesh="               << nEntityTotal(EDGE_CODIM);
        log_stream << " nFacePerMesh="               << nEntityTotal(FACE_CODIM);
        log_stream << " nElementPerMesh="            << nEntityTotal(ELEMENT_CODIM);

        log_stream << std::endl; "    *** ";
        log_stream << " nCorner="                    << nEntity(VERTEX_CODIM);
        log_stream << " nCornerInterior="            << nEntity(VERTEX_CODIM, PartitionType::InteriorEntity);
        log_stream << " nCornerBorder="              << nEntity(VERTEX_CODIM, PartitionType::BorderEntity);
        log_stream << " nCornerGhost="               << nEntity(VERTEX_CODIM, PartitionType::GhostEntity);

        log_stream << std::endl; "    *** ";
        log_stream << " nEdge="                      << nEntity(EDGE_CODIM);
        log_stream << " nEdgeInterior="              << nEntity(EDGE_CODIM, PartitionType::InteriorEntity);
        log_stream << " nEdgeBorder="                << nEntity(EDGE_CODIM, PartitionType::BorderEntity);
        log_stream << " nEdgeGhost="                 << nEntity(EDGE_CODIM, PartitionType::GhostEntity);

        log_stream << std::endl; "    *** ";
        log_stream << " nFace="                      << nEntity(FACE_CODIM);
        log_stream << " nFaceInterior="              << nEntity(FACE_CODIM, PartitionType::InteriorEntity);
        log_stream << " nFaceBoundarySegment="       << nEntity(FACE_CODIM, PartitionType::InteriorEntity, DOMAIN_BOUNDARY_TYPE);
        log_stream << " nFaceBorder="                << nEntity(FACE_CODIM, PartitionType::BorderEntity);
        log_stream << " nFaceGhost="                 << nEntity(FACE_CODIM, PartitionType::GhostEntity);

        log_stream << std::endl; "    *** ";
        log_stream << " nElement="                   << nEntity(ELEMENT_CODIM);
        log_stream << " nElementInterior="           << nEntity(ELEMENT_CODIM, PartitionType::InteriorEntity);
        log_stream << " nElementGhost="              << nEntity(ELEMENT_CODIM, PartitionType::GhostEntity);

        LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, log_stream.str());
        LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "...CurvilinearGridBase: Finished generating curvilinear mesh]]");
    }



    /* ***************************************************************************
     * Section: Setting user constants
     * ***************************************************************************/

    /** \brief Sets the geometry tolerance. Geometry tolerance determines the precision of curvilinear volumes
     *  Returned by the grid
     *  */
    void   geometryRelativeTolerance(double tolerance)  { gridstorage_.GEOMETRY_TOLERANCE = tolerance; }

    /** \brief Retrieves the geometry tolerance */
    double geometryRelativeTolerance() const            { return gridstorage_.GEOMETRY_TOLERANCE; }

    MPIHelper & mpiHelper() const { return mpihelper_; }

    bool withGhostElements() const { return gridstorage_.withGhostElements_; }




    /****************************************************************************
     * Section: (Fake) Refinement Methods of the mesh
     *
     * \note that there is no refinement functionality implemented at the moment
     ****************************************************************************/

    int entityLevel(int codim, LocalIndexType localIndex)  const { return 0; }

    // Checks if the grid is located on a single process
    bool isSerial () { return size_ == 1; }


    /* ***************************************************************************
     * Section: Methods of the mesh
     * ***************************************************************************/

    /** Get total number of entities in a mesh  */
    int nEntityTotal(int codim) const  { return gridstorage_.nEntityTotal_[codim]; }


    /** Get total number of entities on this process
     *
     * \note Currently interpolation points are not counted towards this number.
     * One should not use this number to loop over entities
     * */
    int nEntity(int codim) const  { return gridstorage_.entityAllIndexSet_[codim].size(); }


    /** Get total number of entities of specific type on this process  */
    int nEntity(int codim, PartitionType ptype, StructuralType btype = NO_BOUNDARY_TYPE) const
    {
    	return entityIndexSetSelect(codim, ptype, btype).size();
    }


    /** Get the GeometryType of entities on this process  */
    Dune::GeometryType entityGeometryType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : return Dune::GeometryType(Dune::Impl::SimplexTopology< 0 >::type::id);   break;
    	case EDGE_CODIM    : return Dune::GeometryType(Dune::Impl::SimplexTopology< 1 >::type::id);   break;
    	case FACE_CODIM    : return gridstorage_.face_[localIndex].geometryType;                                  break;
    	case ELEMENT_CODIM : return gridstorage_.element_[localIndex].geometryType;                               break;
    	}
    }


    /** Get physical tag based on codimension  */
    PhysicalTagType physicalTag(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : return 0;                                              break;
    	case EDGE_CODIM    : return 0;                                              break;
    	case FACE_CODIM    : return gridstorage_.face_[localIndex].physicalTag;     break;
    	case ELEMENT_CODIM : return gridstorage_.element_[localIndex].physicalTag;  break;
    	}
    }


    /** Finds global index using local index and codimension of entity. Returns false if requested local index does not correspond to an entity on this process  */
    bool findEntityGlobalIndex(int codim, LocalIndexType localIndex, GlobalIndexType & globalIndex) const
    {
    	if (localIndex < 0)  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Received negative index");  }

    	switch(codim)
    	{
    	case VERTEX_CODIM  : if (localIndex >= gridstorage_.point_.size())    { return false; }  else { globalIndex = gridstorage_.point_[localIndex].globalIndex; }    break;
    	case EDGE_CODIM    : if (localIndex >= gridstorage_.edge_.size())     { return false; }  else { globalIndex = gridstorage_.edge_[localIndex].globalIndex; }     break;
    	case FACE_CODIM    : if (localIndex >= gridstorage_.face_.size())     { return false; }  else { globalIndex = gridstorage_.face_[localIndex].globalIndex; }     break;
    	case ELEMENT_CODIM : if (localIndex >= gridstorage_.element_.size())  { return false; }  else { globalIndex = gridstorage_.element_[localIndex].globalIndex; }  break;
    	}
    	return true;
    }


    /** Finds local index using global index and codimension of entity. Returns false if requested global index does not correspond to an entity on this process  */
    bool findEntityLocalIndex(int codim, GlobalIndexType globalIndex, LocalIndexType & localIndex) const
    {
    	if (globalIndex < 0)  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Received negative index");  }

    	Global2LocalConstIterator tmpIter = gridstorage_.entityIndexMap_[codim].find(globalIndex);
    	if (tmpIter != gridstorage_.entityIndexMap_[codim].end())  { localIndex = (*tmpIter).second;  return true; }
    	else  { return false; }
    }


    /** Checks if the entity with specified codim and local index exists, and if it has the specified structtype  */
    /*
    bool verifyEntity(int codim, LocalIndexType localIndex, StructuralType structtype) const
    {
    	LocalIndexSet & thisSet = entityIndexSetSelect(codim, structtype);

    	if (thisSet.find(localIndex) == thisSet.end())  { return false; }
    	else  { return true; }
    }
    */


    /** Returns the globalId of the entity  */
    IdType globalId(int codim, LocalIndexType localIndex) const
    {
    	GlobalIndexType globalIndex;
    	if (! findEntityGlobalIndex(codim, localIndex, globalIndex))  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: requested local index does not point to an entity");  }
    	return IdType(globalIndex);
    }


    /** Get boundary segment index of this entity if it is a Domain Boundary Face */
    LocalIndexType boundarySegmentIndex(LocalIndexType localIndex) const
    {
    	assert(gridstorage_.face_[localIndex].boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary);
    	return gridstorage_.boundarySegmentIndexMap_.at(localIndex);
    }


    /** Get DomainBoundaryFace index given its boundary segment index */
    // [TODO] The map in this case can be replaced by a simple array/vector, since index is continuous
    LocalIndexType boundarySegment2LocalIndex(LocalIndexType boundarySegmentIndex) const
    {
    	assert(boundarySegmentIndex < gridstorage_.boundarySegment2LocalIndexMap_.size());  // Boundary segment index is continuous
    	return gridstorage_.boundarySegment2LocalIndexMap_.at(boundarySegmentIndex);
    }


    /** Get Structural Type of an entity */
    PartitionType entityPartitionType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : assert(localIndex < gridstorage_.point_.size());    return gridstorage_.point_[localIndex].ptype;    break;
    	case EDGE_CODIM    : assert(localIndex < gridstorage_.edge_.size());     return gridstorage_.edge_[localIndex].ptype;     break;
    	case FACE_CODIM    : assert(localIndex < gridstorage_.face_.size());     return gridstorage_.face_[localIndex].ptype;     break;
    	case ELEMENT_CODIM : assert(localIndex < gridstorage_.element_.size());  return gridstorage_.element_[localIndex].ptype;  break;
    	default :
    	{
    		LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridBase: unexpected codim " + std::to_string(codim));
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected codimension");
    		break;
    	}
    	}
    }

    /** Get Boundary type of an face */
    StructuralType faceBoundaryType(LocalIndexType localIndex) const
    {
    	return gridstorage_.face_[localIndex].boundaryType;
    }


    /** Get Interpolatory order of the entity */
    InterpolatoryOrderType entityInterpolationOrder(int codim, LocalIndexType localIndex) const
    {
    	if (codim >= VERTEX_CODIM)  { return 0; }

    	LocalIndexType localAssocElementIndex;

    	switch (codim)
    	{
    	case ELEMENT_CODIM : localAssocElementIndex = localIndex;  break;
    	case FACE_CODIM    : localAssocElementIndex = gridstorage_.face_[localIndex].element1Index;  break;
    	case EDGE_CODIM    : localAssocElementIndex = gridstorage_.edge_[localIndex].elementIndex;  break;
    	}

    	return gridstorage_.element_[localAssocElementIndex].interpOrder;
    }


    /** Get unique local index that enumerates only corners, disregarding other interpolatory vertices
     * The only purpose of this index is to satisfy dune-standard for consecutive corner index, it is not used internally in the grid
     * */
    LocalIndexType cornerUniqueLocalIndex(LocalIndexType localVertexIndex) const {

    	typedef typename Local2LocalMap::const_iterator  Local2LocalConstIter;
    	Local2LocalConstIter tmp = gridstorage_.cornerIndexMap_.find(localVertexIndex);

    	if (tmp == gridstorage_.cornerIndexMap_.end())
    	{
    		std::cout << "Error: CurvilinearGridBase: CornerUniqueIndex: Vertex with localIndex=" << localVertexIndex << " is not marked as a corner" << std::endl;
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected vertex local index for a corner");
    	}

    	return (*tmp).second;
    }


    // The reverse of cornerUniqueLocalIndex. Given the unique corner in
    LocalIndexType cornerUnique2LocalIndex(LocalIndexType localCornerIndex) const {

    	typedef typename Local2LocalMap::const_iterator  Local2LocalConstIter;
    	Local2LocalConstIter tmp = gridstorage_.cornerIndexMapRev_.find(localCornerIndex);

    	if (tmp == gridstorage_.cornerIndexMapRev_.end())
    	{
    		std::cout << "Error: CurvilinearGridBase: CornerUnique2LocalIndex: Vertex with cornerIndex=" << localCornerIndex << " is not marked as a corner" << std::endl;
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected vertex local index for a corner");
    	}

    	return (*tmp).second;
    }



    /** Generates a vector of local coordinates of corners of requested entity
     * The order of corners corresponds to the one provided by the ref.subEntity()
     *   */
    std::vector<LocalIndexType> entityCornerLocalIndex(int codim, LocalIndexType entityLocalIndex) const
	{
    	std::vector<LocalIndexType> rez;

    	if (codim == VERTEX_CODIM)  { rez.push_back(entityLocalIndex);  return rez; }

    	LocalIndexType    localElementIndex;
    	InternalIndexType subentityIndex;

    	switch(codim)
    	{
    	case ELEMENT_CODIM : localElementIndex = entityLocalIndex;  break;
    	case FACE_CODIM : {
    		localElementIndex = gridstorage_.face_[entityLocalIndex].element1Index;
    		subentityIndex  = gridstorage_.face_[entityLocalIndex].element1SubentityIndex;

    	} break;
    	case EDGE_CODIM : {
    		localElementIndex = gridstorage_.edge_[entityLocalIndex].elementIndex;
    		subentityIndex  = gridstorage_.edge_[entityLocalIndex].subentityIndex;
    	} break;
    	}


    	const EntityStorage & thisElem = gridstorage_.element_[localElementIndex];
        std::vector<LocalIndexType> elementCornerLocalIndexSet = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ct, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

        // If we are interested in corners of the element, we are interested in all corners
        if (codim == ELEMENT_CODIM)  { return elementCornerLocalIndexSet; }

        // Otherwise, we need to calculate the subset of corners wrt selected subentity
        Dune::GeometryType elemGT;  elemGT.makeSimplex(cdim);
        int thisEntityCornerNumber = Dune::ReferenceElements<ct,cdim>::general(elemGT).size(0, codim, cdim);

        for (int i = 0; i < thisEntityCornerNumber; i++)  {
        	//std::cout << "-- Attempting ref.subentity() using gt.dim=" << elemGT.dim() << " cdim=" << cdim << " codim="<<codim << " subcodim="<< VERTEX_CODIM << " subIndex=" << subentityIndex << " subsubindex=" << i << std::endl;

        	InternalIndexType thisCornerSubIndex = Dune::ReferenceElements<ct,cdim>::general(elemGT).subEntity(subentityIndex, codim, i, VERTEX_CODIM);
        	rez.push_back(elementCornerLocalIndexSet[thisCornerSubIndex]);
        }

        return rez;
	}


    /** Check if edge is a complex edge */
    bool isComplex(LocalIndexType localIndex) const
    {
    	Local2LocalIterator tmp = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].find(localIndex);
    	if (tmp == gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end())  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected local edge index"); }

    	LocalIndexType edgePBIndex = (*tmp).second;
    	return gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][edgePBIndex].size() > 1;
    }


    /** Get the neighbour ranks of this communication entity  */
    std::vector<int> & commEntityNeighborRankSet(
    	int codim, LocalIndexType localIndex, PartitionType structSend, PartitionType structRecv
    )
	{
    	typedef typename std::map<LocalIndexType,LocalIndexType>::const_iterator Local2LocalConstIter;

    	Local2LocalMap & tmpMap = selectCommMap(codim, structSend);
    	Local2LocalConstIter tmp = tmpMap.find(localIndex);

    	if (tmp == tmpMap.end())  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected local index"); }

    	LocalIndexType entityLocalSubIndex = (*tmp).second;
    	return selectCommRankVector(codim, structSend, structRecv)[entityLocalSubIndex];
	}


    /** Returns the local index of a subentity of a given entity
     *  \param[in] entityIndex              local index of the entity
     *  \param[in] codim                    codimension of the entity
     *  \param[in] subcodim                 codimension of the subentity
     *  \param[in] subentityInternalIndex   subentity internal index wrt entity
     *
     *  \note subcodim > codim required
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
    LocalIndexType subentityLocalIndex (LocalIndexType entityIndex, int codim, int subcodim, InternalIndexType subentityInternalIndex) const
    {
    	// Stage 1) Find this subentity as an element subentity
    	// **************************************************************************
    	Dune::GeometryType tetrahedronGeometry;
    	tetrahedronGeometry.makeTetrahedron();
    	const Dune::ReferenceElement<ct,cdim> & thisRefElement = Dune::ReferenceElements<ct,cdim>::general(tetrahedronGeometry);

    	if (subcodim == codim)  { return entityIndex; }  // In this case return itself as own subentity
    	if (subcodim < codim) {                          // Wrong by definition
    		LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair = (" + std::to_string(codim) + "," + std::to_string(subcodim) + ")");
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair");
    	}

    	LocalIndexType elementLocalIndex;
    	InternalIndexType elementSubentityInternalIndex1;

    	switch (codim)
    	{
    		case ELEMENT_CODIM :
    		{
    			elementLocalIndex = entityIndex;
    			elementSubentityInternalIndex1 = subentityInternalIndex;
    		} break;
    		case FACE_CODIM :
    		{
    			InternalIndexType elementSubentityInternalIndex2 = gridstorage_.face_[entityIndex].element1SubentityIndex;

    			elementLocalIndex = gridstorage_.face_[entityIndex].element1Index;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    		case EDGE_CODIM :
    		{
    			InternalIndexType elementSubentityInternalIndex2 = gridstorage_.edge_[entityIndex].subentityIndex;

    			elementLocalIndex = gridstorage_.edge_[entityIndex].elementIndex;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    	}


    	// Stage 2) Find index of the element subentity
    	// **************************************************************************

    	int rez;

    	switch (subcodim)
    	{
    		// Face
    		case FACE_CODIM   :  rez = gridstorage_.elementSubentityCodim1_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Edge
    		case EDGE_CODIM   :  rez = gridstorage_.elementSubentityCodim2_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Corner
    		case VERTEX_CODIM :
    		{
    			InterpolatoryOrderType interpolationOrder = gridstorage_.element_[elementLocalIndex].interpOrder;
    			InternalIndexType cornerInternalIndex = Dune::CurvilinearGeometryHelper::cornerIndex(tetrahedronGeometry, interpolationOrder, elementSubentityInternalIndex1);
    			rez = gridstorage_.element_[elementLocalIndex].vertexIndexSet[cornerInternalIndex];
    		} break;
    	}

    	//std::cout << "Requested subentity index of entity=" << entityIndex << " codim=" << codim << " subcodim=" << subcodim << " internalIndex=" << subentityInternalIndex << " rez=" << rez << std::endl;
    	return rez;
    }


    /** \brief local index of the element that is neighbour to this edge.
     * [TODO] Bad name, better say EdgeOwner, because this is simply one of the elements that owns this edge
     *  */
    LocalIndexType edgeNeighbor(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.edge_.size());
    	return gridstorage_.edge_[localIndex].elementIndex;
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
     * Conventions of internalNeighborIndex for face types
     * * For Domain Boundary there is only one neighbor
     * * For Process Boundary 2nd neighbor is always the Ghost Element
     * * For Internal Face there is no convention on order of the neighbors
     *
     *  */
    LocalIndexType faceNeighbor(LocalIndexType localIndex, InternalIndexType internalNeighborIndex) const
    {
    	LocalIndexType rez;

    	assert(localIndex < gridstorage_.face_.size());

        switch(internalNeighborIndex)
        {
        case 0 : rez = gridstorage_.face_[localIndex].element1Index;  break;
        case 1 : rez = gridstorage_.face_[localIndex].element2Index;  break;
        default:
        {
        	LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearPostConstructor: Unexpected neighbor subentity index =" + std::to_string(internalNeighborIndex));
        	DUNE_THROW(Dune::IOError, "CurvilinearGrid: faceNeighbor() unexpected neighbor index");  break;
        }
        }

        return rez;
    }


    /** \brief Coordinate of a requested interpolatory vertex
     *  \param[in] localIndex            local vertex index (insertion index)
     * */
    Vertex vertex(LocalIndexType localIndex) const { return gridstorage_.point_[localIndex].coord; }

    /** Storage data related to this entity, except of explicit vertex coordinates
     *  \param[in] codim                 codimension of the entity
     *  \param[in] localIndex            local edge index
     *
     *  Note: This data is stored only partially - vertex indices are extracted from associated element
     * */
    EntityStorage entityData(int codim, LocalIndexType localIndex) const
    {
    	switch (codim)
    	{
    	case VERTEX_CODIM  : return pointData(localIndex);     break;
    	case EDGE_CODIM    : return edgeData(localIndex);      break;
    	case FACE_CODIM    : return faceData(localIndex);      break;
    	case ELEMENT_CODIM : return elementData(localIndex);   break;
    	default:
    	{
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Requested unexpected codim for entity data ");
    	} break;

    	}
    }

    /** \brief Retrieves the geometry class of the entity. This procedure is expensive, especially for cached geometries */
    template<int codim>
    typename Codim<codim>::EntityGeometry
    entityGeometry(LocalIndexType localIndex) const
    {
    	//std::cout << "attempting to create a codim " << codim << " entity with localIndex=" << localIndex << std::endl;
    	return entityGeometryConstructor<codim>(entityData(codim, localIndex));
    }


    /** Return pointer to Octree or 0 if it is not constructed. */
    const CurvilinearLooseOctree & octree() const { return *gridstorage_.octree_; }


    /** \brief Minimal bounding box for set of elements on this process */
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
    	DUNE_THROW(Dune::IOError, "Using non-initialised OCTree");

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
            for (unsigned int i = 0; i < elementIndices.size(); i++) {
            	typename Codim<ELEMENT_CODIM>::EntityGeometry thisGeometry = entityGeometry<ELEMENT_CODIM>(elementIndices[i]);

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

    // Iterator based on local index of the entity
    // Iterator all entities of a given codimension
    IndexSetIterator entityIndexIterator(int codim, LocalIndexType localIndex) const
    {
    	return gridstorage_.entityAllIndexSet_[codim].find(localIndex);
    }

    // Iterator for entities of a given codimension and structural type only
    IndexSetIterator entityIndexIterator(int codim, PartitionType ptype, LocalIndexType localIndex, StructuralType btype = NO_BOUNDARY_TYPE) const
    {
    	return entityIndexSetSelect(codim, ptype, btype).find(localIndex);
    }

    // Iterator for entities of a given codimension and Dune partition type only
    IndexSetIterator entityIndexDuneIterator(int codim, Dune::PartitionIteratorType pitype, LocalIndexType localIndex) const
    {
    	return entityIndexSetDuneSelect(codim, pitype).find(localIndex);
    }


    // Iterators over all entities of a given codimension
    IndexSetIterator entityIndexBegin(int codim)  { return gridstorage_.entityAllIndexSet_[codim].begin(); }

    IndexSetIterator entityIndexEnd(int codim)    { return gridstorage_.entityAllIndexSet_[codim].end(); }

    // Iterators over specific entities of a given codimension
    // This construction allows fast iteration over entities of specific structural type
    IndexSetIterator entityIndexBegin(int codim, PartitionType ptype, StructuralType btype = NO_BOUNDARY_TYPE) const
    {
    	return entityIndexSetSelect(codim, ptype, btype).begin();
    }

    IndexSetIterator entityIndexEnd(int codim, PartitionType ptype, StructuralType btype = NO_BOUNDARY_TYPE) const
    {
    	return entityIndexSetSelect(codim, ptype, btype).end();
    }


    IndexSetIterator entityDuneIndexBegin(int codim, Dune::PartitionIteratorType pitype) const {
    	return entityIndexSetDuneSelect(codim, pitype).begin();
    }

    IndexSetIterator entityDuneIndexEnd(int codim, Dune::PartitionIteratorType pitype)   const {
    	return entityIndexSetDuneSelect(codim, pitype).end();
    }



    const GridStorageType & gridstorage() const { return gridstorage_; }

    /* ***************************************************************************
     * Section: Public Auxiliary Methods
     * ***************************************************************************/

    // Checks if entities of a given codim are allowed to be of a given structural type
    // If not throws an error
    void assertValidCodimStructuralType(int codim, StructuralType ptype) const
    {
    	bool pass = false;

    	pass |= (ptype == Dune::PartitionType::InteriorEntity);
    	pass |= (ptype == Dune::PartitionType::GhostEntity);

    	if (codim > 0)
    	{
    		// Elements are not allowed to be boundaries
    		pass |= (ptype == Dune::PartitionType::BorderEntity);
    	}


    	if (!pass)  {
    		std::stringstream logstr;
    		logstr << "CurvilinearGridBase: Unexpected codim-structtype pair codim=" << codim;
    		logstr << " ptype=" << ptype;
    		LoggingMessage::template write<CurvGrid::LOG_MSG_PERSISTENT>( __FILE__, __LINE__, logstr.str());
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected codim-structtype pair");
    	}
    }



    /* ***************************************************************************
     * Section: Selector methods for shorthand access of specific arrays
     * ***************************************************************************/

    Local2LocalMap & selectCommMap(int codim, PartitionType ptype)
    {
    	switch (ptype)
    	{
    	case Dune::PartitionType::InteriorEntity :  return gridstorage_.boundaryInternalEntityIndexMap_[codim];  break;
    	case Dune::PartitionType::BorderEntity   :  return gridstorage_.processBoundaryIndexMap_[codim];         break;
    	case Dune::PartitionType::GhostEntity    :  return gridstorage_.ghostIndexMap_[codim];                   break;
    	default                                  :  DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected comm structural type");  break;
    	}
    }


    EntityNeighborRankVector & selectCommRankVector(int codim, PartitionType ptypesend, PartitionType ptyperecv)
    {
    	// Can only communicate over these 3 PartitionTypes
    	//assertValidCodimStructuralType(codim, ptypesend);
    	//assertValidCodimStructuralType(codim, ptyperecv);

    	switch (ptypesend)
    	{
    	case Dune::PartitionType::InteriorEntity :   // Internal -> Ghost protocol
    	{
    		assert(ptyperecv == Dune::PartitionType::GhostEntity);
    		return gridstorage_.BI2GNeighborRank_[codim];
    	} break;
    	case Dune::PartitionType::BorderEntity   :   // PB -> PB and PB -> Ghost protocols
    	{
    		if      (ptyperecv == Dune::PartitionType::BorderEntity)  { return gridstorage_.PB2PBNeighborRank_[codim]; }
    		else if (ptyperecv == Dune::PartitionType::GhostEntity)   { return gridstorage_.PB2GNeighborRank_[codim]; }
    		else { assert(0); }
    	} break;
    	case Dune::PartitionType::GhostEntity    :   // Ghost -> (Internal & PB) and Ghost -> Ghost protocols
    	{
    		if      (ptyperecv == Dune::PartitionType::InteriorEntity) { return gridstorage_.G2BIPBNeighborRank_[codim]; }
    		else if (ptyperecv == Dune::PartitionType::BorderEntity)   { return gridstorage_.G2BIPBNeighborRank_[codim]; }
    		else if (ptyperecv == Dune::PartitionType::GhostEntity)    { return gridstorage_.G2GNeighborRank_[codim]; }
    		else { assert(0); }
    	} break;
    	default: assert(0);  break;

    	}
    }



protected:


    /* ***************************************************************************
     * Section: Auxiliary Methods
     * ***************************************************************************/

    void assertStage(int expectedStage)
    {
    	if ((gridstage_ == Stage::GRID_OPERATION) && (expectedStage == Stage::GRID_CONSTRUCTION)) { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Attempted to insert entities into grid after construction"); }
    }

    // Returns a link to the set of all local indices of entities of a given codimension and specific structural type
    const LocalIndexSet & entityIndexSetSelect(int codim, PartitionType ptype, StructuralType boundaryType) const
    {
    	assertValidCodimStructuralType(codim, ptype);

    	// Check if this is a request for a boundary index set
    	if (boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary)
    	{
    		assert(codim == FACE_CODIM);  // According to convention, only faces can be boundarySegments
    		return gridstorage_.faceDomainBoundaryIndexSet_;
    	}

    	// Otherwise, this must be a request for a standard dune entity iterator
    	switch(ptype)
    	{
    	case Dune::PartitionType::InteriorEntity   : return gridstorage_.entityInternalIndexSet_[codim];          break;
    	case Dune::PartitionType::BorderEntity     : return gridstorage_.entityProcessBoundaryIndexSet_[codim];   break;
    	case Dune::PartitionType::GhostEntity      : return gridstorage_.entityGhostIndexSet_[codim];             break;
    	}
    }

    // Returns a link to the set of all local indices of entities of a given codimension, based on Dune-convention partition type
    const LocalIndexSet & entityIndexSetDuneSelect(int codim, Dune::PartitionIteratorType pitype) const
    {
    	const int DuneIPartition   = Dune::PartitionIteratorType::Interior_Partition;
    	const int DuneIBPartition  = Dune::PartitionIteratorType::InteriorBorder_Partition;
    	const int DuneGPartition   = Dune::PartitionIteratorType::Ghost_Partition;
    	const int DuneAllPartition = Dune::PartitionIteratorType::All_Partition;

    	switch(pitype)
    	{
    	case DuneIPartition     : return gridstorage_.entityDuneInteriorIndexSet_[codim];          break;
    	case DuneIBPartition    : return gridstorage_.entityDuneInteriorBorderIndexSet_[codim];    break;
    	case DuneGPartition     : return gridstorage_.entityGhostIndexSet_[codim];                 break;
    	case DuneAllPartition   : return gridstorage_.entityAllIndexSet_[codim];                   break;
    	default:
    	{
    		std::cout << "CurvilinearGridBase: Unexpected dune-pitype" << std::endl;
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected dune-pitype");         break;
    	}
    	}
    }

    EntityStorage pointData(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.point_.size());
    	const VertexStorage & thisPointData =  gridstorage_.point_[localIndex];

        EntityStorage thisPoint;
        thisPoint.geometryType.makeVertex();
        thisPoint.globalIndex  = thisPointData.globalIndex;
        thisPoint.ptype        = thisPointData.ptype;
        thisPoint.interpOrder  = 0;                    // Note: Points do not have an interpolation order
        thisPoint.physicalTag  = -1;                   // Note: Points do not have a physical tag
        thisPoint.vertexIndexSet.push_back(localIndex);   // Note: Point has only one vertex, and its index is the point index

        return thisPoint;
    }

    EntityStorage edgeData(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.edge_.size());
    	const EdgeStorage & thisEdgeData =    gridstorage_.edge_[localIndex];
        const EntityStorage & assocElement =  gridstorage_.element_[thisEdgeData.elementIndex];

        EntityStorage thisEdge;
        thisEdge.geometryType.makeLine();
        thisEdge.globalIndex  = thisEdgeData.globalIndex;
        thisEdge.ptype        = thisEdgeData.ptype;
        thisEdge.interpOrder  = assocElement.interpOrder;
        thisEdge.physicalTag  = -1;        // Note: Edges do not have a physical tag

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<InternalIndexType> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, cdim>(assocElement.geometryType, thisEdge.interpOrder, 2, thisEdgeData.subentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(unsigned int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisEdge;
    }

    EntityStorage faceData(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.face_.size());
    	const FaceStorage & thisFaceData = gridstorage_.face_[localIndex];
    	assert(thisFaceData.element1Index < gridstorage_.element_.size());
        const EntityStorage & assocElement = gridstorage_.element_[thisFaceData.element1Index];

        EntityStorage thisFace;
        thisFace.geometryType.makeTriangle();
        thisFace.globalIndex  = thisFaceData.globalIndex;
        thisFace.ptype        = thisFaceData.ptype;
        //if (thisFaceData.boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary)  { thisFace.ptype = BOUNDARY_SEGMENT_PARTITION_TYPE; }

        thisFace.interpOrder  = assocElement.interpOrder;
        thisFace.physicalTag  = thisFaceData.physicalTag;

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<InternalIndexType> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, cdim>(assocElement.geometryType, thisFace.interpOrder, 1, thisFaceData.element1SubentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(unsigned int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

        return thisFace;
    }

    EntityStorage elementData(LocalIndexType localIndex) const {
    	assert(localIndex < gridstorage_.element_.size());
    	return gridstorage_.element_[localIndex];
    }


    // Get curved geometry of an entity
    // TODO: assert mydim == element geometry type dim

    template<int codim>
    typename Codim<codim>::EntityGeometry
    entityGeometryConstructor(EntityStorage thisData) const
    {
    	assert(thisData.geometryType.dim() == cdim - codim);

        std::vector<Vertex> entityVertices;
        for (unsigned int i = 0; i < thisData.vertexIndexSet.size(); i++) {
        	LocalIndexType thisIndex = thisData.vertexIndexSet[i];
        	entityVertices.push_back(gridstorage_.point_[thisIndex].coord);
        }

        return typename Codim<codim>::EntityGeometry (thisData.geometryType, entityVertices, thisData.interpOrder);
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


private: // Private members

    // The stage of the grid determines if the grid has already been assembled or not
    int gridstage_;

    // Curvilinear Grid Constructor Class
    GridConstructorType * gridconstructor_;

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
