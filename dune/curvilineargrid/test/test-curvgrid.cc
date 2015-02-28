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

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkcommunicate.cc>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/test/checkiterators.cc>
#include <dune/grid/test/checkadaptation.cc>
#include <dune/grid/test/checkpartition.cc>


template<class ctype, int cdim, int order >
struct CurvFactory
{
  typedef Dune::CurvilinearGrid<cdim, cdim, ctype> GridType;

  static GridType * buildGrid(Dune::MPIHelper & mpihelper)
  {
    std::cout << " using curvgrid 32 with order " << order << std::endl << std::endl;

    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";
    const std::string GMSH_FILE_NAME[5] {"sphere32.msh", "sphere32ord2.msh", "sphere32ord3.msh", "sphere32ord4.msh", "sphere32ord5.msh"};

    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME[order - 1];

    bool insertBoundarySegment = true;
    bool withGhostElements = true;
    bool verbose = true;
    bool processVerbose = true;
    bool writeReaderVTKFile = false;

    Dune::CurvilinearGridFactory<ctype, cdim> factory(withGhostElements, verbose, processVerbose, mpihelper);


    int nVertexTotal;
    int nElementTotal;

    Dune::CurvilinearGmshReader< GridType >::read(factory,
                                                            filename,
                                                            mpihelper,
                                                            nVertexTotal,
                                                            nElementTotal,
                                                            verbose,
                                                            processVerbose,
                                                            writeReaderVTKFile,
                                                            insertBoundarySegment);

    return factory.createGrid(nVertexTotal, nElementTotal);
  }

};


template <class ctype, int cdim>
void check_grid(Dune::CurvilinearGrid<cdim, cdim, ctype> & grid) {
  std::cout << std::endl << "CurvGrid<" << cdim << ">" << std::endl;


  gridcheck(grid);

  checkIterators ( grid.leafGridView() );
  checkIterators ( grid.levelGridView(0) );

  // check communication interface
  checkCommunication(grid, -1, Dune::dvverb);
  for(int l=0; l<=grid.maxLevel(); ++l)  { checkCommunication(grid,l,Dune::dvverb); }

  // check geometry lifetime
  checkGeometryLifetime( grid.leafGridView() );
  // check the method geometryInFather()
  checkGeometryInFather(grid);
  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);
  // check grid adaptation interface
  checkAdaptRefinement(grid);
  checkPartitionType( grid.leafGridView() );

}


int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);
	typedef Dune::CurvilinearGrid<3, 3, double> GridType;

	{
		GridType * grid32ord1 = CurvFactory<double, 3, 1>::buildGrid(mpihelper);
		check_grid(*grid32ord1);
		delete grid32ord1;

	}

	{
		GridType * grid32ord2 = CurvFactory<double, 3, 2>::buildGrid(mpihelper);
		check_grid(*grid32ord2);
		delete grid32ord2;
	}

	{
		GridType * grid32ord3 = CurvFactory<double, 3, 3>::buildGrid(mpihelper);
		check_grid(*grid32ord3);
		delete grid32ord3;
	}

	{
		GridType * grid32ord4 = CurvFactory<double, 3, 4>::buildGrid(mpihelper);
		check_grid(*grid32ord4);
		delete grid32ord4;
	}

	{
		GridType * grid32ord5 = CurvFactory<double, 3, 5>::buildGrid(mpihelper);
		check_grid(*grid32ord5);
		delete grid32ord5;
	}


  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
