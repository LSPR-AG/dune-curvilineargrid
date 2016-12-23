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


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, argc, argv);


    // *********************************************************
	// Check if grid recognises it is parallel
    // *********************************************************
	std::stringstream loggingStream;
	loggingStream << "Is curvilinear grid periodic = " << grid->withPeriodic();
	LoggingMessage::template write<LOG_MSG_PRODUCTION>(__FILE__, __LINE__, loggingStream.str());

	// Report the parallel timing statistics
	LoggingTimer<LoggingMessage>::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
