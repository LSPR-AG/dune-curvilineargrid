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






namespace Dune
{

namespace CurvGrid {


template< class GridType, class HostFactory >
class CurvilinearGridFactoryWrapper
{
  private:

	typedef typename GridType::ctype  ctype;
	typedef typename GridType::GridBaseType            GridBaseType;

	static const int dimension = GridType::dimension;
	typedef FieldVector< ctype, dimension >            Coordinate;

	typedef int    LocalIndexType;
	typedef int    GlobalIndexType;

  public:

	// Real constructor
	CurvilinearGridFactoryWrapper(HostFactory & hostfactory)
  	  : hostfactory_(hostfactory)
    {

    }


    ~CurvilinearGridFactoryWrapper ()  { }

    void insertVertex ( const Coordinate &pos, const GlobalIndexType globalId )
    {
    	hostfactory_.insertVertex(pos);
    	//hostfactory_.insertVertex(pos, globalId);
    }


    void insertElement(
      GeometryType &geometry,
      const std::vector< LocalIndexType > &vertexIndexSet,
      const int elemOrder,
      const int physicalTag)
    {
    	hostfactory_.insertElement(geometry, vertexIndexSet);
    }


    void insertBoundarySegment(
        GeometryType &geometry,
        const std::vector< LocalIndexType > &vertexIndexSet,
        const int elemOrder,
        //const LocalIndexType associatedElementIndex,
        const int physicalTag,
		bool isDomainBoundary)
    {
    	// NOTE: Interior Boundaries do not exist in Dune Core as of time of writing
    	// Note: associatedElementIndex no longer necessary
    	if (isDomainBoundary)  { hostfactory_.insertBoundarySegment(vertexIndexSet); }
    }


    void insertNVertexTotal(int nVertexTotal)  {
    	// GMSH reader knows the total number of vertices, so one may save time by reusing this quantity, not needing extra global communication
    }

    void insertNElementTotal(int nElementTotal)  {
    	// GMSH reader knows the total number of elements, so one may save time by reusing this quantity, not needing extra global communication
    }


    GridType * createGrid()
    {
    	GridType * grid = hostfactory_.createGrid();

    	return grid;
    }


    // Variables
    // -----------------------------------------------------------
  private:

    HostFactory & hostfactory_;



  };

} // namespace CurvGrid

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDFACTORY_HH
