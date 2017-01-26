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
#include <dune/curvilineargrid/common/loggingtimer.hh>

#include <dune/curvilineargrid/curvilineargridbase/impl/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearlooseoctree.hh>

#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbasecommunication.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbasecorner.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseedge.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseentity.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseintersection.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseindexset.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseproperty.hh>
#include <dune/curvilineargrid/curvilineargridbase/methods/curvilinearoctreewrapper.hh>


/* ***************************************************************************
 * Specifications: Parallel Curvilinear Mesh Manager
 *
 * Existent Functionality:
 *     - All vertices, edges, faces, elements are accessible through their globalIndex and as subentities
 *     - All elements and faces possess physicalTag - integer associated to their material property
 *     - InterProcessBoundaries automatically generate
 *     - GhostElements automatically communicated and generated (optional)
 *     - GlobalIndex automatically generated for edges, faces and elements
 *     - Supports cuboid periodic domains with some or all dimensions periodic. Periodic ghosts are communicated as regular ghosts and not translated
 *     - [not tested] OCTree functionality allows to find element and its local coordinate corresponding to a given global coordinate
 *     - [not tested] Supports non-uniformly p-refined meshes (not tested)
 *
 *
 * Missing functionality:
 *  - [DESIGN] reader must provide globalIds for vertices and elements, they are not generated automatically.
 *  - [DESIGN] reader must provide all BoundarySegments (process boundaries), they are not communicated automatically
 *  - [TODO] Does NOT support refinement - it is not possible to dynamically refine or coarsen the mesh at the moment
 *  - [TODO] Does NOT support hanging nodes - it is not possible to load non-uniform h-refined mesh at the moment
 *  - [TODO] Does NOT support front/overlap elements at the moment
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
 *  - [TODO] Move insertVertex etc entirely to the GridConstructor. Leave GridBase only with usage methods
 *  - [TODO] There are two many methods in GridBase. Create subclasses to split file size
 *  - [TODO] Disable all the vertex2string output for multiprocessor case - too much output
 *  - [TODO] Implement TimeSync logging message so processes do not comment on top of each other. Check Boost if already exists
 *
 * Testing:
 *  - [TODO] Constructor run check if all non-owned entities have been successfully enumerated at the end
 *
 *
 * Additional Functionality - Extended Physical Tags
 *  - Solution 1: Implement Vector Physical Tags, read from GMSH. Involves changing tag reading in GMSH and tag communication in grid construction
 *  - Solution 2: Keep single integer tag. Involves implementing a tag->info mapper in the derived code.
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

namespace CurvGrid {

// *******************************************
// Forwards-declaration
// *******************************************
template<class ct, int cdim, bool isCached>
class CurvilinearGridStorage;



template <class ct, int cdim, bool isCached>
class CurvilinearGridBase {
public:

    /* public types
     * *******************************************************************/
//
	typedef ct        ctype;
	static const int dimension = cdim;
	static const int dimensionworld = cdim;
	static const int is_cached = isCached;
//
    typedef CurvilinearGridBase<ct, cdim, isCached>        GridBaseType;
    typedef CurvilinearGridStorage<ct, cdim, isCached>   GridStorageType;
//
    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;
//
    typedef typename GridStorageType::GlobalCoordinate                 GlobalCoordinate;
//    typedef typename GridStorageType::VertexStorage          VertexStorage;
//    typedef typename GridStorageType::EdgeStorage            EdgeStorage;
//    typedef typename GridStorageType::FaceStorage            FaceStorage;
//    typedef typename GridStorageType::EntityStorage          EntityStorage;
//
//    typedef typename GridStorageType::EdgeKey                EdgeKey;
//    typedef typename GridStorageType::FaceKey                FaceKey;
    typedef typename GridStorageType::IdType                 IdType;
//
//    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
//    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
//    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;
//    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
//    typedef typename GridStorageType::Global2LocalConstIterator      Global2LocalConstIterator;
//    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;
//
//    typedef typename GridStorageType::NodeType                  NodeType;
    typedef typename GridStorageType::CurvilinearLooseOctree    CurvilinearLooseOctree;
//
//    typedef typename GridStorageType::EntityNeighborRankVector  EntityNeighborRankVector;
//
    typedef Dune::CurvGrid::LoggingMessage         LoggingMessage;
    typedef Dune::CurvGrid::LoggingTimer<LoggingMessage>         LoggingTimer;
//

//
    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;
//
//    static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;
//    static const unsigned int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;
//    static const unsigned int INTERIOR_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::InteriorBoundary;
//    static const unsigned int PERIODIC_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::PeriodicBoundary;

    typedef CurvilinearGridBaseCommunication<GridBaseType> GridComm;
    typedef CurvilinearGridBaseCorner<GridBaseType> GridCorner;
    typedef CurvilinearGridBaseEdge<GridBaseType> GridEdge;
    typedef CurvilinearGridBaseEntity<GridBaseType> GridEntity;
    typedef CurvilinearGridBaseIntersection<GridBaseType> GridIntersection;
    typedef CurvilinearGridBaseIndexSet<GridBaseType> GridIndexSet;
    typedef CurvilinearGridBaseProperty<GridBaseType> GridProperty;
    typedef CurvilinearGridBaseOctreeWrapper<GridBaseType> OctreeWrapper;




public: /* public methods */

    /** Parallel constructor - USE THIS CONSTRUCTOR*/
    CurvilinearGridBase(
    		bool withGhostElements,
			bool withElementGlobalIndex,
			MPIHelper &mpihelper,
			std::vector<bool> periodicCuboidDimensions) :
		gridstorage_(mpihelper, withGhostElements, withElementGlobalIndex, periodicCuboidDimensions),
		property_(*this),
	    comm_(*this),
	    corner_(*this),
	    edge_(*this),
	    entity_(*this),
	    intersection_(*this),
	    indexset_(*this),
		octreewrapper_(*this)
    {
        std::string log_string = "Initialized CurvilinearGridBase withGhostElements=" + std::to_string(withGhostElements);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
    }

private:
    /** Copy constructor: private, undefined: disallow copy */
    CurvilinearGridBase(const CurvilinearGridBase&);

    /** Assignment operator: private, undefined: disallow assignment */
    //CurvilinearGridBase& operator=(const CurvilinearGridBase&);

public:

    GridStorageType & gridstorage() { return gridstorage_; }
    GridProperty & property() { return property_; }
    GridComm & comm() { return comm_; }
    OctreeWrapper & octreeWrapper() const { return octreewrapper_; }

    const GridCorner & corner() const { return corner_; }
    const GridEdge & edge() const { return edge_; }
    const GridEntity & entity() const { return entity_; }
    const GridIntersection & intersection() const { return intersection_; }
    const GridIndexSet & indexset() const { return indexset_; }



    /** \brief Coordinate of a requested interpolatory vertex
     *  \param[in] localIndex            local vertex index (insertion index)
     * */
    GlobalCoordinate vertex(LocalIndexType localIndex) const { return gridstorage_.point_[localIndex].coord; }


    /** Return pointer to Octree or 0 if it is not constructed. */
    const CurvilinearLooseOctree & octree() const { return *gridstorage_.octree_; }


private: // Private members

    // Curvilinear Grid Storage Class
    GridStorageType gridstorage_;

    GridProperty property_;
    GridComm comm_;
    GridCorner corner_;
    GridEdge edge_;
    GridEntity entity_;
    GridIntersection intersection_;
    GridIndexSet indexset_;
    OctreeWrapper octreewrapper_;

};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDBASE_HH
