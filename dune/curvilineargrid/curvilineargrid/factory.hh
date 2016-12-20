// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CURVGRID_GRIDFACTORY_HH
#define DUNE_CURVGRID_GRIDFACTORY_HH

/** \file
 *  \author Aleksejs Fomins
 *  \brief  Implementation of Curvilinear Grid Factory
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

//#include <dune/grid/utility/globalindex.hh>

//#include <dune/alugrid/common/transformation.hh>
//#include <dune/alugrid/3d/alugrid.hh>
//#include <dune/alugrid/3d/gridfactory.hh>
//#include <dune/alugrid/3d/gridfactory.cc>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>




namespace Dune
{

namespace CurvGrid {


template< class GridType >
class CurvilinearGridFactory
{
  private:

	typedef typename GridType::ctype  ctype;
	typedef typename GridType::GridBaseType            GridBaseType;
	typedef typename GridBaseType::LocalIndexType      LocalIndexType;
	typedef typename GridBaseType::GlobalIndexType     GlobalIndexType;

	static const int dimension = GridType::dimension;
	typedef FieldVector< ctype, dimension >                 VertexCoordinate;


  public:

	// Real constructor
    CurvilinearGridFactory(
    		bool withGhostElements,
    		bool withGmshElementIndex,
    		MPIHelper &mpihelper,
			std::vector<bool> periodicCuboidDimensions = std::vector<bool>())
    {
    	gridbase_ = new GridBaseType(
    		withGhostElements,
    		withGmshElementIndex,
    		mpihelper,
			periodicCuboidDimensions);
    }

    // GridBase is constructed and deleted here
    // NOTE: FACTORY SHOULD NOT DELETE GRIDBASE, BECAUSE THAT KILLS THE GRID IF THE FACTORY IS DESTROYED BEFORE THE GRID
    ~CurvilinearGridFactory ()  {
    	//if (gridbase_)  { delete gridbase_; }
    }

    void insertVertex ( const VertexCoordinate &pos, const GlobalIndexType globalIndex )
    {
    	gridbase_->insertVertex(pos, globalIndex);
    }

    void insertElement(
      GeometryType &geometry,
      const std::vector< LocalIndexType > &vertexIndexSet,
      const GlobalIndexType globalIndex,
      const int elemOrder,
      const int physicalTag)
    {
    	gridbase_->insertElement(geometry, vertexIndexSet, globalIndex, elemOrder, physicalTag);
    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const std::vector< LocalIndexType > &vertexIndexSet,
        const int elemOrder,
        //const LocalIndexType associatedElementIndex,
        const int physicalTag,
		bool isDomainBoundary)
    {
    	// Note: associatedElementIndex no longer necessary
    	// gridbase_->insertBoundarySegment(geometry, associatedElementIndex, vertexIndexSet, elemOrder, physicalTag, isDomainBoundary);
    	gridbase_->insertBoundarySegment(geometry, vertexIndexSet, elemOrder, physicalTag, isDomainBoundary);
    }

    void insertNVertexTotal(int nVertexTotal)  { gridbase_->insertNVertexTotal(nVertexTotal); }

    void insertNElementTotal(int nElementTotal)  { gridbase_->insertNElementTotal(nElementTotal); }


    GridType * createGrid()
    {
    	gridbase_->generateMesh();
    	GridType * grid = new GridType(gridbase_);

    	return grid;
    }


    // Variables
    // -----------------------------------------------------------
  private:

    GridBaseType * gridbase_;



  };

} // namespace CurvGrid

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDFACTORY_HH
