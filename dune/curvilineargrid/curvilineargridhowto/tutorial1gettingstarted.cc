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

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached, Dune::LoggingMessage> GridType;
	GridType * grid = createGrid<GridType>(mpihelper);



    // **********************************************
    // Do all sort of fancy things with the grid here
    // **********************************************
	std::cout << "the number of grid nElement=" << grid->size(0) << " rank_ = " << grid->comm().rank() << " size_ is " << grid->comm().size() << std::endl;

	typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::getInstance().reportParallel(mpihelper);

    // Delete the grid
    delete grid;


    return 0;
}
