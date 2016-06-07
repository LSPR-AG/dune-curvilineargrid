#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>

#include <dune/curvilineargrid/utility/globalboundarycontainer.hh>





/** \brief
 *
 */


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	const int grid_file_type = 5;  // createGrid procedure provides 6 different example grids numbered 0 to 5

	const bool isCached = true;
	const int ELEMENT_CODIM = 0;  // Codimension of element in 3D

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, grid_file_type);
	GridType & gridRef = (*grid);


	typedef Dune::CurvGrid::GlobalBoundaryContainer<GridType> BoundaryContainer;
	BoundaryContainer testContainer(gridRef);






	typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
