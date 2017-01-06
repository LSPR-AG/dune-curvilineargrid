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
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbasefactory.hh>

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

	typedef CurvilinearGridBaseFactory<GridBaseType> BaseFactory;


  public:

	// Real constructor
    CurvilinearGridFactory(
    		bool withGhostElements,
    		bool withGmshElementIndex,
    		MPIHelper &mpihelper,
			std::vector<bool> periodicCuboidDimensions = std::vector<bool>())
    {
    	basefactory_ = new BaseFactory(withGhostElements, withGmshElementIndex, mpihelper, periodicCuboidDimensions);
    }

    void insertVertex ( const VertexCoordinate &pos, const GlobalIndexType globalIndex )
    {
    	assert(basefactory_);
    	basefactory_->insertVertex(pos, globalIndex);
    }

    void insertElement(
      GeometryType &geometry,
      const std::vector< LocalIndexType > &vertexIndexSet,
      const GlobalIndexType globalIndex,
      const int elemOrder,
      const int physicalTag)
    {
    	assert(basefactory_);
    	basefactory_->insertElement(geometry, vertexIndexSet, globalIndex, elemOrder, physicalTag);
    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const std::vector< LocalIndexType > &vertexIndexSet,
        const int elemOrder,
        const int physicalTag,
		bool isDomainBoundary)
    {
    	assert(basefactory_);
    	basefactory_->insertBoundarySegment(geometry, vertexIndexSet, elemOrder, physicalTag, isDomainBoundary);
    }

    void insertNVertexTotal(int nVertexTotal)  { basefactory_->insertNVertexTotal(nVertexTotal); }

    void insertNElementTotal(int nElementTotal)  { basefactory_->insertNElementTotal(nElementTotal); }

    // NOTE: USER MUST DELETE THE GRID POINTER
    GridType * createGrid()  {
    	assert(basefactory_);
    	grid_ = new GridType(basefactory_->createGrid());
    	delete basefactory_;
    	return grid_;
    }


  private:

    BaseFactory * basefactory_;
    GridType * grid_;

};

} // namespace CurvGrid

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDFACTORY_HH
