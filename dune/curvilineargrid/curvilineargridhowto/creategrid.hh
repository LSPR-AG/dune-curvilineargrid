#include <cstdlib>

#include <dune/curvilineargrid/curvilineargrid/factory.hh>
#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>




/** \brief This example procedure constructs a Curvilinear Grid from a .msh file
 *
 *   To do so it initialises the LoggingMessage and LoggingTimer classes, and then uses the CurvilinearGridFactory and CurvilinearGMSHReader to construct the grid.
 *
 *  */
template <class GridType>
GridType * createGrid(Dune::MPIHelper & mpihelper, int argc , char **argv)
{
	// Obtain path to the mesh folder from a CMake constant
    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

    // A choice of example meshes. First number is the number of elements, 2nd is the polynomial order of the mesh
    const std::vector<std::string> GMSH_FILE_NAME {
    	"sphere32.msh",
    	"sphere32ord2.msh",
    	"sphere32ord3.msh",
    	"sphere32ord4.msh",
    	"sphere32ord5.msh",
    	"sphere2000ord3.msh",
		"sphere2400ord1.msh",
		"sphere2layer.msh",
		"cube575.msh",
		"squarePeriodic2.msh"
    };


    // Check if correct arguments are provided
	if (argc <= 1) {
		std::cout << "Usage: " << std::endl;
		std::cout << "  ./tutorialfile meshindex [periodicnormals] " << std::endl;
		std::cout << "    -meshindex - index of the mesh to use in this tutorial [0.." << GMSH_FILE_NAME.size()-1 << "]" << std::endl;
		std::cout << "    -periodic normals - optional number in binary to determine which dimensions are periodic. For example, 101 will mean that X and Z are periodic, and Y is not" << std::endl;
		exit(1);
	}

    // Choice of file name
	int exampleFile = atoi(argv[1]);
    assert((exampleFile >= 0)&&(exampleFile < GMSH_FILE_NAME.size()));
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME[exampleFile];

    // Determine if periodicity is used and which one
	std::vector<bool> periodicDim;
	if (argc >2) {
		assert(std::string(argv[2]).size() == 3);

		char * periodicend;
		long int periodicint = strtol(argv[2], &periodicend, 2);

		periodicDim = {
				bool(periodicint & 4),
				bool(periodicint & 2),
				bool(periodicint & 1)
		};
	}

    // Additional constants
    bool withGhostElements    = true;      // to create Ghost elements
    bool withGmshElementIndex = true;      // to re-use Global Element Index provided by GMSH (recommended)

    // Initialize LoggingMessage and LoggingTimer
    Dune::CurvGrid::LoggingMessage::init(mpihelper);
    Dune::CurvGrid::LoggingTimer<Dune::CurvGrid::LoggingMessage>::init(mpihelper);

    // Construct the grid factory
    Dune::CurvGrid::CurvilinearGridFactory<GridType> factory(withGhostElements, withGmshElementIndex, mpihelper, periodicDim);

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvGrid::CurvilinearGmshReader< GridType >::read(factory, filename, mpihelper);

    // Create the grid
    return factory.createGrid();
}
