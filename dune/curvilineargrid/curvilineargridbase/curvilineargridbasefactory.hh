// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CURVGRID_GRIDBASEFACTORY_HH
#define DUNE_CURVGRID_GRIDBASEFACTORY_HH

/** \file
 *  \author Aleksejs Fomins
 *  \brief  Implementation of Curvilinear Grid Base Factory
 */

#include <config.h>

#include <map>
#include <limits>
#include <vector>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>




namespace Dune
{


template< class GridBaseType >
class CurvilinearGridBaseFactory
{
  private:

	// Typedefs and const variables
	// -----------------------------------------------------------
	typedef typename GridBaseType::ctype      ctype;

	static const bool  iscached     = GridBaseType::is_cached;
    static const int dimension      = GridBaseType::dimension;
    static const int dimensionworld = GridBaseType::dimensionworld;

	typedef FieldVector< ctype, dimensionworld >                           VertexCoordinate;

	typedef int VertexGlobalIndex;
	typedef int VertexLocalIndex;
	typedef int ElementLocalIndex;
	typedef int ElementGlobalIndex;

  public:

    CurvilinearGridBaseFactory(
    	bool withGhostElements,
    	bool withGmshElementIndex,
    	MPIHelper &mpihelper)
    {
    	gridbase_ = new GridBaseType(withGhostElements, withGmshElementIndex, mpihelper);
    }

    ~CurvilinearGridBaseFactory ()  {}

    void insertVertex ( const VertexCoordinate &pos, const VertexGlobalIndex globalIndex )
    {
    	gridbase_->insertVertex(pos, globalIndex);
    }

    void insertElement(
      GeometryType &geometry,
      const std::vector< VertexLocalIndex > &vertexIndexSet,
      const ElementGlobalIndex globalIndex,
      const int elemOrder,
      const int physicalTag)
    {
    	gridbase_->insertElement(geometry, vertexIndexSet, globalIndex, elemOrder, physicalTag);
    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const std::vector< VertexLocalIndex > &vertexIndexSet,
        const int elemOrder,
        //const ElementLocalIndex associatedElementIndex,
        const int physicalTag,
		bool isDomainBoundary)
    {
    	// Note: associatedElementIndex no longer necessary
    	//gridbase_->insertBoundarySegment(geometry, associatedElementIndex, vertexIndexSet, elemOrder, physicalTag, isDomainBoundary);
    	gridbase_->insertBoundarySegment(geometry, vertexIndexSet, elemOrder, physicalTag, isDomainBoundary);
    }


    void insertNVertexTotal(int nVertexTotal)  { gridbase_->insertNVertexTotal(nVertexTotal); }

    void insertNElementTotal(int nElementTotal)  { gridbase_->insertNElementTotal(nElementTotal); }


    GridBaseType * createGrid()
    {
    	gridbase_->generateMesh();
    	return gridbase_;
    }

    // Variables
    // -----------------------------------------------------------
  private:

    GridBaseType * gridbase_;

  };

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDBASEFACTORY_HH
