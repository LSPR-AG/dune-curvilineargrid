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


const bool isCached = true;


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	const int grid_file_type = 1;  // createGrid procedure provides 6 different example grids numbered 0 to 5

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, 2);


    // *********************************************************
    // Do all sort of fancy things with the grid here
	// In this example we collect some information and
	// use the LoggingMessage to write it on the master process
    // *********************************************************
	std::stringstream loggingStream;
	loggingStream << "the number of grid elements = " << grid->size(0);
	loggingStream << " process rank = " << grid->comm().rank();
	loggingStream << " process size = " << grid->comm().size();
	Dune::LoggingMessage::template write<Dune::CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, loggingStream.str());

	// Report the parallel timing statistics
	Dune::LoggingTimer<Dune::LoggingMessage>::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
