// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>


using namespace Dune;
using namespace Dune::CurvGrid;

const bool isCached = true;
const int DIM0D = 0;   const int CODIM0D = 3;
const int DIM1D = 1;   const int CODIM1D = 2;
const int DIM2D = 2;   const int CODIM2D = 1;
const int DIM3D = 3;   const int CODIM3D = 0;

// example for a generic algorithm that traverses
// the entities of a given mesh in various ways
template<int codim, class GridType>
void traversal (GridType& grid)
{
  const int dim =  GridType::dimension;
  typedef typename GridType::ctype ct;
  typedef typename GridType::LeafGridView LeafGridView;

  typedef typename GridType::LocalIndexType          LocalIndexType;
  typedef typename GridType::GlobalIndexType         GlobalIndexType;
  typedef typename GridType::PhysicalTagType         PhysicalTagType;
  typedef typename GridType::InterpolatoryOrderType  InterpolatoryOrderType;


  std::string dim2name[4] = {"vertices", "edges", "faces", "elements" };
  std::cout << "Iterating over all entities of codim=" << codim << ", that is, over " << dim2name[dim - codim] << std::endl;


  // get the instance of the LeafGridView
  LeafGridView leafView = grid.leafGridView();

  typedef typename LeafGridView::template Codim<codim>::Iterator EntityLeafIterator;
  //typedef typename LeafGridView::template Codim<codim>::Geometry EntityGeometry;
  typedef typename GridType::template Codim<codim>::EntityGeometryMappingImpl  BaseGeometry;
  typedef typename BaseGeometry::GlobalCoordinate                              GlobalCoordinate;

  for (auto&& elementThis : entities(leafView, Dune::Dim<dim - codim>()))
  {
	Dune::GeometryType gt              = elementThis.type();
	LocalIndexType  localIndex         = grid.leafIndexSet().index(elementThis);
	GlobalIndexType globalIndex        = grid.template entityGlobalIndex<codim>(elementThis);
	PhysicalTagType physicalTag        = grid.template entityPhysicalTag<codim>(elementThis);
	InterpolatoryOrderType interpOrder = grid.template entityInterpolationOrder<codim>(elementThis);

	// Constructing a geometry is quite expensive, do it only once

	// NOTE: As of time of writing, Dune::Grid::Geometry does not provide the interface to access properties specific to
	// Curvilinear elements. Therefore, CurvilinearGrid provides direct access to CurvilinearGeometry, such that the user
	// has access to all its features.

	//EntityGeometry geom = it->geometry();
	BaseGeometry geom = grid.template entityBaseGeometry<codim>(elementThis);
	std::vector<GlobalCoordinate>  interpVertices = geom.vertexSet();

	std::cout << "visiting leaf " << gt
              << " localIndex=" << localIndex
              << " globalIndex=" << globalIndex
              << " physicalTag=" << physicalTag
              << " interpolationOrder=" << interpOrder
              << " consisting of interpolatory vertices " << VectorHelper::vector2string(interpVertices)
              << std::endl;
  }
}




int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dimension = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<ctype, dimension, isCached> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, argc , argv);

	// Traverse all entities of the grid and write information about each entity
	traversal<CODIM3D, GridType>(*grid);  // Elements
	traversal<CODIM2D, GridType>(*grid);  // Faces
	traversal<CODIM1D, GridType>(*grid);  // Edges
	traversal<CODIM0D, GridType>(*grid);  // Corners

	typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}

