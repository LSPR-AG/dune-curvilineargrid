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
 *  - [TODO] Does NOT support overlapping elements at the moment
 *  - [TODO] Does NOT support non-tetrahedral meshes. Generalization to arbitrary and mixed geometry meshes is possible but will be cumbersome
 *
 * Development log
 *  - [TODO]  Consider making GridConstructor a pointer and deleting it at the end of construction phase.
 *  - [FIXME] Convert all code to unsigned ints wherever the number is definitely non-negative
 *  - [FIXME] Do Ghost elements require same properties as normal elements. For example, subentityByIndex?
 *  - [FIXME] Implement ostream operator << for IdType.
 *  - [FIXME] Currently nEntity counts ghost elements. Is this expected
 *  - [FIXME] Currently nEntity counts only corners. Is this expected
 *
 *  - [FIXME] Need to add normal and outerNormal
 *  - [FIXME] Need to wrap for Dune
 *  - [FIXME] Need to match Dune's internal subentity id convention
 *  - [FIXME] When returning Ghost elements, must check if they are defined, and throw error if not
 *  - [FIXME] Replace all DUNE_THROW output by the logging message output such that it can be read
 *
 * Usage:
 *  - [TODO] Disable all the vertex2string output for multiprocessor case - too much output
 *  - [TODO] Implement timing (and perhaps memory usage) for each of the construction operations
 *  - [TODO] Implement memory in LoggingMessage as is done in Hades
 *  - [TODO] Implement TimeSync logging message so processes do not comment on top of each other. Check Boost if already exists
 *
 *
 * Testing:
 *  - [FIXME] Write test which checks consistency of all local and global indices
 *  - [FIXME] Constructor run check if all non-owned entities have been successfully enumerated at the end
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



template <class ct, int cdim>
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
    typedef Dune::CurvilinearGridBase<ct, cdim>          GridBaseType;
    typedef Dune::CurvilinearGridStorage<ct, cdim>       GridStorageType;
    typedef Dune::CurvilinearGridConstructor<ct, cdim>   GridConstructorType;

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

    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;

    typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;

    typedef typename GridStorageType::template Codim<0>::EntityGeometry    ElementGeometry;


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
        log_stream << " nVertexPerMesh="             << nEntityTotal(3);
        log_stream << " nEdgePerMesh="               << nEntityTotal(2);
        log_stream << " nFacePerMesh="               << nEntityTotal(1);
        log_stream << " nElementPerMesh="            << nEntityTotal(0);
        log_stream << " nVertex="                    << nEntity(3);
        log_stream << " nEdge="                      << nEntity(2);
        log_stream << " nFace="                      << nEntity(1);
        log_stream << " nElement="                   << nEntity(0);
        log_stream << " nInternalElement="           << nEntity(0, InternalType);
        log_stream << " nGhostElement="              << nEntity(0, GhostType);
        log_stream << " nFaceDomainBoundary="        << nEntity(1, DomainBoundaryType);
        log_stream << " nFaceProcessBoundary="       << nEntity(1, ProcessBoundaryType);
        log_stream << " nFaceInternal="              << nEntity(1, InternalType);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_stream.str());
    }




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
    int nEntityTotal(int codim) const  { gridstorage_.nEntityTotal_[codim]; }


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
    	case 3 : return Dune::GeometryType(Dune::GenericGeometry::SimplexTopology<0>::type::id, 0);   break;
    	case 2 : return Dune::GeometryType(Dune::GenericGeometry::SimplexTopology<1>::type::id, 1);   break;
    	case 1 : return gridstorage_.face_[localIndex].geometryType;                                  break;
    	case 0 : return gridstorage_.element_[localIndex].geometryType;                               break;
    	}
    }


    /** Get physical tag based on codimension  */
    PhysicalTagType physicalTag(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case 3 : return 0;                                              break;
    	case 2 : return 0;                                              break;
    	case 1 : return gridstorage_.face_[localIndex].physicalTag;     break;
    	case 0 : return gridstorage_.element_[localIndex].physicalTag;  break;
    	}
    }


    /** Finds global index using local index and codimension of entity. Returns false if requested local index does not correspond to an entity on this process  */
    bool findEntityGlobalIndex(int codim, LocalIndexType localIndex, GlobalIndexType & globalIndex) const
    {
    	if (localIndex < 0)  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Received negative index");  }

    	switch(codim)
    	{
    	case 3 : if (localIndex >= gridstorage_.point_.size())    { return false; }  else { globalIndex = gridstorage_.point_[localIndex].globalIndex; }    break;
    	case 2 : if (localIndex >= gridstorage_.edge_.size())     { return false; }  else { globalIndex = gridstorage_.edge_[localIndex].globalIndex; }     break;
    	case 1 : if (localIndex >= gridstorage_.face_.size())     { return false; }  else { globalIndex = gridstorage_.face_[localIndex].globalIndex; }     break;
    	case 0 : if (localIndex >= gridstorage_.element_.size())  { return false; }  else { globalIndex = gridstorage_.element_[localIndex].globalIndex; }  break;
    	}
    	return true;
    }


    /** Finds local index using global index and codimension of entity. Returns false if requested global index does not correspond to an entity on this process  */
    bool findEntityLocalIndex(int codim, GlobalIndexType globalIndex, LocalIndexType & localIndex) const
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


    IdType globalId(int codim, LocalIndexType localIndex)
    {
    	GlobalIndexType globalIndex;
    	if (! findEntityGlobalIndex(codim, localIndex, globalIndex))  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: requested local index does not point to an entity");  }

    	IdType thisId;
    	thisId.id_ = std::pair<StructuralType, GlobalIndexType>( entityStructuralType(codim, localIndex),  globalIndex );
    	return thisId;
    }


    /** Get Structural Type of an entity */
    StructuralType entityStructuralType(int codim, LocalIndexType localIndex) const
    {
    	switch(codim)
    	{
    	case 3 : return gridstorage_.point_[localIndex].structuralType;                       break;
    	case 2 : return gridstorage_.edge_[localIndex].structuralType;                        break;
    	case 1 : return gridstorage_.face_[localIndex].structuralType;                        break;
    	case 0 : return gridstorage_.element_[localIndex].structuralType;                     break;
    	default : DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected subentity codimension");  break;
    	}
    }


    /** Get Interpolatory order of the entity */
    InterpolatoryOrderType entityInterpolationOrder(int codim, LocalIndexType localIndex) const
    {
    	if (codim >= 3)  { return 0; }

    	LocalIndexType localAssocElementIndex;

    	switch (codim)
    	{
    	case 0 : localAssocElementIndex = localIndex;  break;
    	case 1 : localAssocElementIndex = gridstorage_.face[localIndex].element1Index;  break;
    	case 2 : localAssocElementIndex = gridstorage_.edge[localIndex].elementIndex;  break;
    	}

    	return gridstorage_.element[localAssocElementIndex].interpOrder;
    }


    /** Get vertex global coordinate */
    std::vector<LocalIndexType> entityCornerLocalIndex(int codim, LocalIndexType localIndex) const
	{
    	std::vector<LocalIndexType> rez;

    	if (codim == 3)  { rez.push_back(localIndex);  return rez; }

    	LocalIndexType    localElementIndex;
    	InternalIndexType subentityIndex;

    	switch(codim)
    	{
    	case 0 : localElementIndex = localIndex;  break;
    	case 1 : {
    		localElementIndex = gridstorage_.face_[localIndex].element1Index;
    		subentityIndex  = gridstorage_.face_[localIndex].element1SubentityIndex;

    	} break;
    	case 2 : {
    		localElementIndex = gridstorage_.edge_[localIndex].elementIndex;
    		subentityIndex  = gridstorage_.edge_[localIndex].subentityIndex;
    	} break;
    	}


    	const EntityStorage & thisElem = gridstorage_.element_[localElementIndex];
        std::vector<LocalIndexType> elementCornerLocalIndexSet = Dune::CurvilinearGeometryHelper::entityVertexCornerSubset<ct, 3>(thisElem.geometryType, thisElem.vertexIndexSet, thisElem.interpOrder);

        std::vector<InternalIndexType> internalLinearSubentityIndices = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(thisElem.geometryType, 2, subentityIndex);

        for (int i = 0; i < internalLinearSubentityIndices.size(); i++)  {
        	rez.push_back(elementCornerLocalIndexSet[internalLinearSubentityIndices[i]]);
        }

        return rez;
	}


    /** Check if edge is a complex edge */
    bool isComplex(LocalIndexType localIndex) const
    {
    	LocalIndexType edgePBIndex = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM][localIndex];
    	return gridstorage_.processBoundaryNeighborRank_[EDGE_CODIM][edgePBIndex].size() > 1;
    }


    /** Get the neighbors of this process boundary  */

    std::vector<int> processBoundaryNeighborRankSet(int codim, LocalIndexType localIndex) const
	{
    	if ((codim <= 0)||(codim > 3)) { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected process boundary codim"); }

    	LocalIndexType entityPBIndex = gridstorage_.processBoundaryIndexMap_[codim][localIndex];
    	return gridstorage_.processBoundaryNeighborRank_[codim][entityPBIndex];
	}




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
    LocalIndexType subentityLocalIndex (LocalIndexType entityIndex, int codim, int subcodim, InternalIndexType subentityInternalIndex)
    {
    	// Stage 1) Find this subentity as an element subentity
    	// **************************************************************************
    	Dune::GeometryType tetrahedronGeometry;
    	tetrahedronGeometry.makeTetrahedron();
    	Dune::ReferenceElement<ct,cdim> & thisRefElement = Dune::ReferenceElements<ct,cdim>::general(tetrahedronGeometry);

    	if (subcodim >= codim) { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: subentityIndex(): Unexpected codim-subcodim pair"); }

    	LocalIndexType elementLocalIndex;
    	InternalIndexType elementSubentityInternalIndex1;

    	switch (codim)
    	{
    		case 0 :
    		{
    			elementLocalIndex = entityIndex;
    			elementSubentityInternalIndex1 = subentityInternalIndex;
    		} break;
    		case 1 :
    		{
    			InternalIndexType elementSubentityInternalIndex2 = gridstorage_.face_[entityIndex].element1SubentityIndex;

    			elementLocalIndex = gridstorage_.face_[entityIndex].element1Index;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    		case 2 :
    		{
    			InternalIndexType elementSubentityInternalIndex2 = gridstorage_.edge[entityIndex].subentityIndex;

    			elementLocalIndex = gridstorage_.edge_[entityIndex].element1Index;
    			elementSubentityInternalIndex1 = thisRefElement.subEntity(elementSubentityInternalIndex2, codim, subentityInternalIndex, subcodim);
    		} break;
    	}


    	// Stage 2) Find index of the element subentity
    	// **************************************************************************

    	switch (subcodim)
    	{
    		// Face
    		case 1 :  return gridstorage_.elementSubentityCodim1_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Edge
    		case 2 :  return gridstorage_.elementSubentityCodim2_[elementLocalIndex][elementSubentityInternalIndex1];  break;
    		// Corner
    		case 3 :
    		{
    			InterpolatoryOrderType interpolationOrder = gridstorage_.element_[elementLocalIndex].interpOrder;
    			InternalIndexType cornerInternalIndex = Dune::CurvilinearGeometryHelper::cornerID(tetrahedronGeometry, interpolationOrder, elementSubentityInternalIndex1);
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
     *
     *
     * Conventions of internalNeighborIndex for face types
     * * For Domain Boundary there is only one neighbor
     * * For Process Boundary 2nd neighbor is always the Ghost Element
     * * For Internal Face there is no convention on order of the neighbors
     * FIXME: Check if this contradicts Dune convention in any way
     *
     *  */
    LocalIndexType faceNeighbor(LocalIndexType localIndex, InternalIndexType internalNeighborIndex) const
    {
    	LocalIndexType rez;

        switch(internalNeighborIndex)
        {
        case 0 : rez = gridstorage_.face_[localIndex].element1Index;  break;
        case 1 : rez = gridstorage_.face_[localIndex].element2Index;  break;
        default: DUNE_THROW(Dune::IOError, "CurvilinearGrid: faceNeighbor() unexpected neighbor index");  break;
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
    	case 3 : return pointData(localIndex);     break;
    	case 2 : return edgeData(localIndex);      break;
    	case 1 : return faceData(localIndex);      break;
    	case 0 : return elementData(localIndex);   break;
    	default:
    	{
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Requested unexpected codim for entity data ");
    	} break;

    	}
    }


    template<int codim>
    typename GridStorageType::template Codim<codim>::EntityGeometry
    entityGeometry(LocalIndexType localIndex) const
    {
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
                ElementGeometry thisGeometry = entityGeometry<0>(elementIndices[i]);

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
    IndexSetIterator entityIndexIterator(int codim, LocalIndexType localIndex)
    {
    	return gridstorage_.entityAllIndexSet_[codim].find(localIndex);
    }

    // Iterator for entities of a given codimension and structural type only
    IndexSetIterator entityIndexIterator(int codim, StructuralType structtype, LocalIndexType localIndex)
    {
    	return entityIndexSetSelect(codim, structtype).find(localIndex);
    }

    // Iterator for entities of a given codimension and Dune partition type only
    IndexSetIterator entityIndexDuneIterator(int codim, Dune::PartitionIteratorType pitype, LocalIndexType localIndex)
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


    IndexSetIterator entityDuneIndexBegin(int codim, Dune::PartitionIteratorType pitype)  { return entityIndexSetDuneSelect(codim, pitype).begin(); }

    IndexSetIterator entityDuneIndexEnd(int codim, Dune::PartitionIteratorType pitype)    { return entityIndexSetDuneSelect(codim, pitype).end(); }


    /* ***************************************************************************
     * Section: Public Auxiliary Methods
     * ***************************************************************************/

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
    	default: DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected dune-pitype");         break;
    	}
    }

    EntityStorage pointData(LocalIndexType localIndex) const
    {
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
    	const FaceStorage & thisFaceData = gridstorage_.face_[localIndex];
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

    EntityStorage elementData(LocalIndexType localIndex) const { return gridstorage_.element_[localIndex]; }


    // Get curved geometry of an entity
    // TODO: assert mydim == element geometry type dim

    template<int codim>
    typename GridStorageType::template Codim<codim>::EntityGeometry
    entityGeometryConstructor(EntityStorage thisData) const
    {
        std::vector<Vertex> entityVertices;
        for (int i = 0; i < thisData.vertexIndexSet.size(); i++) { entityVertices.push_back(gridstorage_.point_[thisData.vertexIndexSet[i]].coord); }

        return Dune::CurvilinearGeometry<ct, cdim - codim, cdim> (thisData.geometryType, entityVertices, thisData.interpOrder);
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
