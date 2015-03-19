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

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
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
 *  - [FIXME] Currently nEntity counts ghost elements. Is this expected
 *  - [FIXME] Currently nEntity counts only corners. Is this expected
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
    typedef Dune::CurvilinearGridBase<ct, cdim, isCached>         GridBaseType;
    typedef Dune::CurvilinearGridStorage<ct, cdim, isCached>      GridStorageType;
    typedef Dune::CurvilinearGridConstructor<ct, cdim, isCached>  GridConstructorType;

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
    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;

    typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;

    typedef typename GridStorageType::EntityNeighborRankVector  EntityNeighborRankVector;

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


    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV      = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    // Partition type shorthands
    static const unsigned int DomainBoundaryType   = GridStorageType::PartitionType::DomainBoundary;
    static const unsigned int ProcessBoundaryType  = GridStorageType::PartitionType::ProcessBoundary;
    static const unsigned int InternalType         = GridStorageType::PartitionType::Internal;
    static const unsigned int GhostType            = GridStorageType::PartitionType::Ghost;



public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridBase(bool withGhostElements, bool verbose, bool processVerbose, MPIHelper &mpihelper ) :
        verbose_(verbose),
        processVerbose_(processVerbose),
        gridstage_(0),
        mpihelper_(mpihelper),
        gridstorage_(withGhostElements),
        gridconstructor_(verbose, processVerbose, gridstorage_, *this, mpihelper)
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
    	gridconstructor_.insertVertex(p, globalIndex);
    }

    /** \brief Insert an element into the mesh
     * \param[in] gt               geometry type of the element (hopefully tetrahedron)
     * \param[in] globalId         Index unique for the union of all faces and elements of the GMSH file
     * \param[in] vertexIndexSet   local indices of interpolatory vertices. Local with respect to the order the vertices were inserted
     * \param[in] order            interpolatory order of the element
     * \param[in] physicalTag      physical tag of the element
     *
     * */
    void insertElement(
        	Dune::GeometryType gt,
        	GlobalIndexType globalId,
        	const std::vector<LocalIndexType> & vertexIndexSet,
        	InterpolatoryOrderType order,
        	PhysicalTagType physicalTag)
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
     *  \param[in] globalId                 Index unique for the union of all faces and elements of the GMSH file
     *  \param[in] associatedElementIndex   local index of the element this face is associated to
     *  \param[in] vertexIndexSet           local indices of the interpolatory vertices of this face
     *  \param[in] order                    interpolatory order of the face
     *  \param[in] physicalTag              physical tag of the element (material property)
     *
     * */

    void insertBoundarySegment(
    		Dune::GeometryType gt,
        	GlobalIndexType globalId,
        	LocalIndexType associatedElementIndex,
        	const std::vector<LocalIndexType> & vertexIndexSet,
        	InterpolatoryOrderType order,
        	PhysicalTagType physicalTag)
    {
    	assertStage(Stage::GRID_CONSTRUCTION);
    	gridconstructor_.insertBoundarySegment(gt, globalId, associatedElementIndex, vertexIndexSet, order, physicalTag);
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

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Initializing mesh");

        gridconstructor_.generateMesh();

        // Diagnostics output
        std::stringstream log_stream;
        log_stream << "CurvilinearGridBase: Constructed Mesh ";
        log_stream << " nVertexPerMesh="             << nEntityTotal(VERTEX_CODIM);
        log_stream << " nEdgePerMesh="               << nEntityTotal(EDGE_CODIM);
        log_stream << " nFacePerMesh="               << nEntityTotal(FACE_CODIM);
        log_stream << " nElementPerMesh="            << nEntityTotal(ELEMENT_CODIM);

        log_stream << std::endl; "    *** ";
        log_stream << " nCorner="                    << nEntity(VERTEX_CODIM);
        log_stream << " nCornerInternal="            << nEntity(VERTEX_CODIM, InternalType);
        log_stream << " nCornerDomainBoundary="      << nEntity(VERTEX_CODIM, DomainBoundaryType);
        log_stream << " nCornerProcessBoundary="     << nEntity(VERTEX_CODIM, ProcessBoundaryType);
        log_stream << " nCornerGhost="               << nEntity(VERTEX_CODIM, GhostType);

        log_stream << std::endl; "    *** ";
        log_stream << " nEdge="                      << nEntity(EDGE_CODIM);
        log_stream << " nEdgeInternal="              << nEntity(EDGE_CODIM, InternalType);
        log_stream << " nEdgeDomainBoundary="        << nEntity(EDGE_CODIM, DomainBoundaryType);
        log_stream << " nEdgeProcessBoundary="       << nEntity(EDGE_CODIM, ProcessBoundaryType);
        log_stream << " nEdgeGhost="                 << nEntity(EDGE_CODIM, GhostType);

        log_stream << std::endl; "    *** ";
        log_stream << " nFace="                      << nEntity(FACE_CODIM);
        log_stream << " nFaceInternal="              << nEntity(FACE_CODIM, InternalType);
        log_stream << " nFaceDomainBoundary="        << nEntity(FACE_CODIM, DomainBoundaryType);
        log_stream << " nFaceProcessBoundary="       << nEntity(FACE_CODIM, ProcessBoundaryType);
        log_stream << " nFaceGhost="                 << nEntity(FACE_CODIM, GhostType);

        log_stream << std::endl; "    *** ";
        log_stream << " nElement="                   << nEntity(ELEMENT_CODIM);
        log_stream << " nInternalElement="           << nEntity(ELEMENT_CODIM, InternalType);
        log_stream << " nGhostElement="              << nEntity(ELEMENT_CODIM, GhostType);

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
    }



    /* ***************************************************************************
     * Section: Setting user constants
     * ***************************************************************************/

    // Sets the geometry tolerance
    void   geometryRelativeTolerance(double tolerance)  { gridstorage_.GEOMETRY_TOLERANCE = tolerance; }

    // Retrieves the geometry tolerance
    double geometryRelativeTolerance() const            { return gridstorage_.GEOMETRY_TOLERANCE; }

    bool verbose() const  { return verbose_; }

    bool processVerbose() const  { return processVerbose_; }

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
    int nEntity(int codim, StructuralType structtype) const  { return entityIndexSetSelect(codim, structtype).size(); }


    /** Get the GeometryType of entities on this process  */
    Dune::GeometryType entityGeometryType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : return Dune::GeometryType(Dune::GenericGeometry::SimplexTopology<0>::type::id, 0);   break;
    	case EDGE_CODIM    : return Dune::GeometryType(Dune::GenericGeometry::SimplexTopology<1>::type::id, 1);   break;
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
    bool findEntityLocalIndex(int codim, GlobalIndexType globalIndex, LocalIndexType & localIndex)
    {
    	if (globalIndex < 0)  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Received negative index");  }

    	Global2LocalIterator tmpIter = gridstorage_.entityIndexMap_[codim].find(globalIndex);
    	if (tmpIter != gridstorage_.entityIndexMap_[codim].end())  { localIndex = (*tmpIter).second;  return true; }
    	else  { return false; }
    }


    /** Checks if the entity with specified codim and local index exists, and if it has the specified structtype  */
    bool verifyEntity(int codim, LocalIndexType localIndex, StructuralType structtype) const
    {
    	LocalIndexSet & thisSet = entityIndexSetSelect(codim, structtype);

    	if (thisSet.find(localIndex) == thisSet.end())  { return false; }
    	else  { return true; }
    }


    IdType globalId(int codim, LocalIndexType localIndex) const
    {
    	GlobalIndexType globalIndex;
    	if (! findEntityGlobalIndex(codim, localIndex, globalIndex))  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: requested local index does not point to an entity");  }

    	IdType thisId;
    	thisId.id_ = std::pair<StructuralType, GlobalIndexType>( entityStructuralType(codim, localIndex),  globalIndex );
    	return thisId;
    }


    /** Get boundary segment index of this entity if it is a Domain Boundary Face */
    LocalIndexType boundarySegmentIndex(LocalIndexType localIndex) const
    {
    	StructuralType thisstructtype = entityStructuralType(FACE_CODIM, localIndex);
    	assert(thisstructtype == DomainBoundaryType);
    	return gridstorage_.boundarySegmentIndexMap_.at(localIndex);
    }


    /** Get Structural Type of an entity */
    StructuralType entityStructuralType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case VERTEX_CODIM  : assert(localIndex < gridstorage_.point_.size());    return gridstorage_.point_[localIndex].structuralType;    break;
    	case EDGE_CODIM    : assert(localIndex < gridstorage_.edge_.size());     return gridstorage_.edge_[localIndex].structuralType;     break;
    	case FACE_CODIM    : assert(localIndex < gridstorage_.face_.size());     return gridstorage_.face_[localIndex].structuralType;     break;
    	case ELEMENT_CODIM : assert(localIndex < gridstorage_.element_.size());  return gridstorage_.element_[localIndex].structuralType;  break;
    	default :
    	{
    		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: unexpected codim " + std::to_string(codim));
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected codimension");
    		break;
    	}
    	}
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


    /** Get vertex global coordinate */
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
    	int codim, LocalIndexType localIndex, StructuralType structSend, StructuralType structRecv
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
    		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair = (" + std::to_string(codim) + "," + std::to_string(subcodim) + ")");
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

    	switch (subcodim)
    	{
    		// Face
    		case FACE_CODIM   :  return gridstorage_.elementSubentityCodim1_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Edge
    		case EDGE_CODIM   :  return gridstorage_.elementSubentityCodim2_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Corner
    		case VERTEX_CODIM :
    		{
    			InterpolatoryOrderType interpolationOrder = gridstorage_.element_[elementLocalIndex].interpOrder;
    			InternalIndexType cornerInternalIndex = Dune::CurvilinearGeometryHelper::cornerIndex(tetrahedronGeometry, interpolationOrder, elementSubentityInternalIndex1);
    			return gridstorage_.element_[elementLocalIndex].vertexIndexSet[cornerInternalIndex];
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
        	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearPostConstructor: Unexpected neighbor subentity index =" + std::to_string(internalNeighborIndex));
        	DUNE_THROW(Dune::IOError, "CurvilinearGrid: faceNeighbor() unexpected neighbor index");  break;
        }
        }

        return rez;
    }


    /** Vertex coordinate
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


    template<int codim>
    typename Codim<codim>::EntityGeometry
    entityGeometry(LocalIndexType localIndex) const
    {
    	//std::cout << "attempting to create a codim " << codim << " entity with localIndex=" << localIndex << std::endl;
    	return entityGeometryConstructor<codim>(entityData(codim, localIndex));
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
            for (int i = 0; i < elementIndices.size(); i++) {
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
    IndexSetIterator entityIndexIterator(int codim, StructuralType structtype, LocalIndexType localIndex) const
    {
    	return entityIndexSetSelect(codim, structtype).find(localIndex);
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
    IndexSetIterator entityIndexBegin(int codim, StructuralType structtype)  { return entityIndexSetSelect(codim, structtype).begin(); }

    IndexSetIterator entityIndexEnd(int codim, StructuralType structtype)    { return entityIndexSetSelect(codim, structtype).end(); }


    IndexSetIterator entityDuneIndexBegin(int codim, Dune::PartitionIteratorType pitype) const { return entityIndexSetDuneSelect(codim, pitype).begin(); }

    IndexSetIterator entityDuneIndexEnd(int codim, Dune::PartitionIteratorType pitype)   const { return entityIndexSetDuneSelect(codim, pitype).end(); }


    /* ***************************************************************************
     * Section: Public Auxiliary Methods
     * ***************************************************************************/

    Dune::PartitionType structural2PartitionType(StructuralType structtype) const
    {
    	switch (structtype)
    	{
    	  case GridStorageType::PartitionType::Internal          : return PartitionType::InteriorEntity;
    	  case GridStorageType::PartitionType::DomainBoundary    : return PartitionType::InteriorEntity;
    	  case GridStorageType::PartitionType::ProcessBoundary   : return PartitionType::BorderEntity;
    	  case GridStorageType::PartitionType::Ghost             : return PartitionType::GhostEntity;
    	  default : DUNE_THROW(Dune::IOError, "CurvilinearGridBase: unexpected structural type for conversion");  break;
    	}
    }


    std::string PartitonTypeName(StructuralType structtype) const  { return gridstorage_.PartitonTypeName[structtype]; }

    // Checks if entities of a given codim are allowed to be of a given structural type
    // If not throws an error
    void assertValidCodimStructuralType(int codim, StructuralType structtype) const
    {
    	bool fail = false;

    	fail |= (structtype == GridStorageType::PartitionType::FrontBoundary);    // Not Implemented
    	fail |= (structtype == GridStorageType::PartitionType::Overlap);          // Not Implemented
    	fail |= (structtype == GridStorageType::PartitionType::InternalBoundary); // Not Implemented

    	if (codim == 0)
    	{
    		// Elements are not allowed to be boundaries
    		fail |= ((structtype == GridStorageType::PartitionType::ProcessBoundary));
    		fail |= ((structtype == GridStorageType::PartitionType::DomainBoundary));
    	}

    	if (fail)  {
    		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearGridBase: Unexpected codim-structtype pair");
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected codim-structtype pair");
    	}
    }



    /* ***************************************************************************
     * Section: Selector methods for shorthand access of specific arrays
     * ***************************************************************************/

    Local2LocalMap & selectCommMap(int codim, StructuralType structtype)
    {
    	switch (structtype)
    	{
    	case InternalType          :  return gridstorage_.boundaryInternalEntityIndexMap_[codim];  break;
    	case ProcessBoundaryType   :  return gridstorage_.processBoundaryIndexMap_[codim];         break;
    	case GhostType             :  return gridstorage_.ghostIndexMap_[codim];                   break;
    	default                    :  DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected comm structural type");  break;
    	}
    }


    EntityNeighborRankVector & selectCommRankVector(int codim, StructuralType structSend, StructuralType structRecv)
    {
    	// Can only communicate over these 3 PartitionTypes
    	assert((structRecv == InternalType)||(structRecv == ProcessBoundaryType)||(structRecv == GhostType));

    	switch (structSend)
    	{
    	case InternalType          :   // Internal -> Ghost protocol
    	{
    		assert(structRecv == GhostType);
    		return gridstorage_.BI2GNeighborRank_[codim];
    	} break;
    	case ProcessBoundaryType   :   // PB -> PB and PB -> Ghost protocols
    	{
    		assert((structRecv == ProcessBoundaryType)||(structRecv == GhostType));
    		if (structRecv == ProcessBoundaryType)  { return gridstorage_.PB2PBNeighborRank_[codim]; }
    		if (structRecv == GhostType)            { return gridstorage_.PB2GNeighborRank_[codim]; }
    	} break;
    	case GhostType             :   // Ghost -> (Internal & PB) and Ghost -> Ghost protocols
    	{
    		assert((structRecv == InternalType)||(structRecv == ProcessBoundaryType)||(structRecv == GhostType));
    		if (structRecv == InternalType)         { return gridstorage_.G2BIPBNeighborRank_[codim]; }
    		if (structRecv == ProcessBoundaryType)  { return gridstorage_.G2BIPBNeighborRank_[codim]; }
    		if (structRecv == GhostType)            { return gridstorage_.G2GNeighborRank_[codim]; }
    	} break;

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
    const LocalIndexSet & entityIndexSetSelect(int codim, StructuralType structtype) const
    {
    	assertValidCodimStructuralType(codim, structtype);

    	switch(structtype)
    	{
    	case InternalType          : return gridstorage_.entityInternalIndexSet_[codim];          break;
    	case DomainBoundaryType    : return gridstorage_.entityDomainBoundaryIndexSet_[codim];    break;
    	case ProcessBoundaryType   : return gridstorage_.entityProcessBoundaryIndexSet_[codim];   break;
    	case GhostType             : return gridstorage_.entityGhostIndexSet_[codim];             break;
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
        thisPoint.globalIndex     = thisPointData.globalIndex;
        thisPoint.structuralType  = thisPointData.structuralType;
        thisPoint.interpOrder     = 0;                    // Note: Points do not have an interpolation order
        thisPoint.physicalTag     = -1;                   // Note: Points do not have a physical tag
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
        thisEdge.globalIndex     = thisEdgeData.globalIndex;
        thisEdge.structuralType  = thisEdgeData.structuralType;
        thisEdge.interpOrder     = assocElement.interpOrder;
        thisEdge.physicalTag     = -1;        // Note: Edges do not have a physical tag

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<InternalIndexType> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, cdim>(assocElement.geometryType, thisEdge.interpOrder, 2, thisEdgeData.subentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisEdge.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

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
        thisFace.globalIndex     = thisFaceData.globalIndex;
        thisFace.structuralType  = thisFaceData.structuralType;
        thisFace.interpOrder     = assocElement.interpOrder;
        thisFace.physicalTag     = thisFaceData.physicalTag;

        // Get the internal element vertex indices associated with this face as a subentity
        std::vector<InternalIndexType> subentityVertexIndices =
            Dune::CurvilinearGeometryHelper::subentityInternalCoordinateSet<ct, cdim>(assocElement.geometryType, thisFace.interpOrder, 1, thisFaceData.element1SubentityIndex);

        // Calculate the localIndices's of vertices of this face by extracting them from the element vertex Ids
        for(int i = 0; i < subentityVertexIndices.size(); i++) { thisFace.vertexIndexSet.push_back(assocElement.vertexIndexSet[subentityVertexIndices[i]]); }

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

    	//std::cout << "started geom constructor " << thisData.geometryType << " " << thisData.interpOrder << std::endl;

        std::vector<Vertex> entityVertices;
        for (int i = 0; i < thisData.vertexIndexSet.size(); i++) { entityVertices.push_back(gridstorage_.point_[thisData.vertexIndexSet[i]].coord); }

        //std::cout << "  -- assembled vertices " << Dune::VectorHelper::vector2string(entityVertices) << std::endl;#

        return typename Codim<codim>::EntityGeometry (thisData.geometryType, entityVertices, thisData.interpOrder);

        //return typename Codim<codim>::EntityGeometry (thisData.geometryType, entityVertices, thisData.interpOrder);
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

    bool verbose_;
    bool processVerbose_;


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
