// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/factory.hh>

#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkcommunicate.hh>
#include <dune/grid/test/checkgeometryinfather.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkadaptation.hh>
#include <dune/grid/test/checkpartition.hh>

const bool isGeometryCached = true;

using namespace Dune;
using namespace Dune::CurvGrid;


template<class GridType, int order>
struct CurvFactory
{

  static GridType * buildGrid(Dune::MPIHelper & mpihelper)
  {
    std::cout << " using curvgrid 32 with order " << order << std::endl << std::endl;

    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";
    const std::string GMSH_FILE_NAME[5] {"sphere32.msh", "sphere32ord2.msh", "sphere32ord3.msh", "sphere32ord4.msh", "sphere32ord5.msh"};

    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + "sphere2000ord3.msh"; //GMSH_FILE_NAME[order - 1];

    bool withGhostElements    = true;
    bool withGmshElementIndex = true;

    CurvilinearGridFactory< GridType > factory(withGhostElements, withGmshElementIndex, mpihelper);
    CurvilinearGmshReader< GridType >::read(factory, filename, mpihelper);

    return factory.createGrid();
  }

};



template <class GridType>
void check_grid(GridType & grid) {
  std::cout << "CurvGrid<" << GridType::dimension << ">" << std::endl;


  //std::cout << "-- Running Base GridCheck" << std::endl;
  //gridcheck(grid);

  //std::cout << "-- Checking LeafGridView" << std::endl;
  //checkIterators ( grid.leafGridView() );

  //std::cout << "-- Checking LevelGridView" << std::endl;
  //checkIterators ( grid.levelGridView(0) );

  // check communication interface
  std::cout << "-- Checking Communication" << std::endl;
  checkCommunication(grid, -1, Dune::dvverb);
  for(int l=0; l<=grid.maxLevel(); ++l)  { checkCommunication(grid,l,Dune::dvverb); }

  //std::cout << "-- Checking Geometry Lifetime" << std::endl;
  //checkGeometryLifetime( grid.leafGridView() );

  //std::cout << "-- Checking Geometry in Father" << std::endl;
  //checkGeometryInFather(grid);

  //std::cout << "-- Checking Intersection Iterator" << std::endl;
  //checkIntersectionIterator(grid);

  //std::cout << "-- Checking Adaptive Refinement" << std::endl;
  //checkAdaptRefinement(grid);

  //std::cout << "-- Checking Partition Type of LeafGridView" << std::endl;
  //checkPartitionType( grid.leafGridView() );

}


int main (int argc , char **argv) {
  try {
	// Initialize MPI, if present
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

    // Instantiation of the logging message and loggingtimer
    typedef LoggingTimer<LoggingMessage>     LoggingTimer;
    LoggingMessage::init(mpihelper);
    LoggingTimer::init(mpihelper);

	typedef Dune::CurvilinearGrid<double, 3, isGeometryCached> GridType;


    /*
	{
		GridType * grid32ord1 = CurvFactory<GridType, 1>::buildGrid(mpihelper);
		check_grid(*grid32ord1);
		delete grid32ord1;

	}
*/
	{
		GridType * grid32ord2 = CurvFactory<GridType, 2>::buildGrid(mpihelper);
		check_grid(*grid32ord2);
		delete grid32ord2;
	}
	/*
	{
		GridType * grid32ord3 = CurvFactory<GridType, 3>::buildGrid(mpihelper);
		check_grid(*grid32ord3);
		delete grid32ord3;
	}

	{
		GridType * grid32ord4 = CurvFactory<GridType, 4>::buildGrid(mpihelper);
		check_grid(*grid32ord4);
		delete grid32ord4;
	}

	{
		GridType * grid32ord5 = CurvFactory<GridType, 5>::buildGrid(mpihelper);
		check_grid(*grid32ord5);
		delete grid32ord5;
	}
	*/

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
