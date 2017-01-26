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
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilineargridconstructor.hh>






namespace Dune
{

namespace CurvGrid {


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

    typedef CurvilinearGridConstructor<GridBaseType>  GridConstructor;

	typedef FieldVector< ctype, dimensionworld >                           VertexCoordinate;

	typedef int VertexGlobalIndex;
	typedef int VertexLocalIndex;
	typedef int ElementLocalIndex;
	typedef int ElementGlobalIndex;

  public:

    CurvilinearGridBaseFactory(
    	bool withGhostElements,
    	bool withGmshElementIndex,
    	MPIHelper &mpihelper,
		std::vector<bool> periodicCuboidDimensions = std::vector<bool>())
    {
    	gridbase_ = new GridBaseType(withGhostElements, withGmshElementIndex, mpihelper, periodicCuboidDimensions);
    	gridconstructor_ = new GridConstructor(*gridbase_, mpihelper);
    }

    // GridBase is constructed and deleted here
    // NOTE: FACTORY SHOULD NOT DELETE GRIDBASE, BECAUSE THAT KILLS THE GRID IF THE FACTORY IS DESTROYED BEFORE THE GRID
    ~CurvilinearGridBaseFactory ()  {
    	//if (gridbase_)  { delete gridbase_; }
    }

    void insertVertex ( const VertexCoordinate &pos, const VertexGlobalIndex globalIndex )
    {
    	assert(gridconstructor_);
    	gridconstructor_->insertVertex(pos, globalIndex);
    }

    void insertElement(
      GeometryType &geometry,
      const std::vector< VertexLocalIndex > &vertexIndexSet,
      const ElementGlobalIndex globalIndex,
      const int elemOrder,
      const int physicalTag)
    {
    	assert(gridconstructor_);
    	gridconstructor_->insertElement(geometry, vertexIndexSet, globalIndex, elemOrder, physicalTag);
    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const std::vector< VertexLocalIndex > &vertexIndexSet,
        const int elemOrder,
        const int physicalTag,
		bool isDomainBoundary)
    {
    	assert(gridconstructor_);
    	gridconstructor_->insertBoundarySegment(geometry, vertexIndexSet, elemOrder, physicalTag, isDomainBoundary);
    }

    // [TODO] This functionality is rudimentary. The mpi_comm of 2 integers is really cheap, no reason to further break facade class factory
    void insertNVertexTotal(int nVertexTotal)  { gridconstructor_->insertNVertexTotal(nVertexTotal); }
    void insertNElementTotal(int nElementTotal)  { gridconstructor_->insertNElementTotal(nElementTotal); }


    // NOTE: USER MUST DELETE THE GRIDBASE POINTER!!!
    GridBaseType * createGrid()
    {
    	assert(gridconstructor_);
    	gridconstructor_->generateMesh();

        // Free up the memory taken by the construction procedure
        delete(gridconstructor_);

    	return gridbase_;
    }

    // Variables
    // -----------------------------------------------------------
  private:

    GridBaseType * gridbase_;
    GridConstructor * gridconstructor_;

  };

} // namespace CurvGrid

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDBASEFACTORY_HH
