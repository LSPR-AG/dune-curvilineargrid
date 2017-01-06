/********************************************
 * Dune-CurvilinearGrid
 * Tutorial 1: Getting started
 *
 * Author: Aleksejs Fomins
 *
 * Description: This simple tutorial reads a curvilinear grid from .msh file,
 * constructs the mesh and reports the timing of reading and construction procedures
 ********************************************/

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>


const bool isCached = true;


using namespace Dune;
using namespace Dune::CurvGrid;


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, argc, argv);


    // *********************************************************
    // Do all sort of fancy things with the grid here
	// In this example we collect some information and
	// use the LoggingMessage to write it on the master process
    // *********************************************************
	std::stringstream loggingStream;
	loggingStream << "the number of grid elements = " << grid->size(0);
	loggingStream << " process rank = " << grid->comm().rank();
	loggingStream << " process size = " << grid->comm().size();
	LoggingMessage::template write<LOG_MSG_PRODUCTION>(__FILE__, __LINE__, loggingStream.str());

	// Report the parallel timing statistics
	LoggingTimer<LoggingMessage>::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
