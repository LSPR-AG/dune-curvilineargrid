#include <dune/curvilineargrid/curvilineargrid/factory.hh>
#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>




/** \brief This example procedure constructs a Curvilinear Grid from a .msh file
 *
 *   To do so it initialises the LoggingMessage and LoggingTimer classes, and then uses the CurvilinearGridFactory and CurvilinearGMSHReader to construct the grid.
 *
 *  */
template <class GridType>
GridType * createGrid(Dune::MPIHelper & mpihelper, int exampleFile)
{
	// Obtain path to the mesh folder from a CMake constant
    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

    // A choice of example meshes. First number is the number of elements, 2nd is the polynomial order of the mesh
    const std::string GMSH_FILE_NAME[8] {
    	"sphere32.msh",
    	"sphere32ord2.msh",
    	"sphere32ord3.msh",
    	"sphere32ord4.msh",
    	"sphere32ord5.msh",
    	"sphere2000ord3.msh",
		"sphere2400ord1.msh",
		"sphere2layer.msh",
    };

    // Choice of file name
    assert((exampleFile >= 0)&&(exampleFile < 8));
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME[exampleFile];

    // Additional constants
    bool withGhostElements    = true;      // to create Ghost elements
    bool withGmshElementIndex = true;      // to re-use Global Element Index provided by GMSH (recommended)

    // Initialize LoggingMessage and LoggingTimer
    Dune::LoggingMessage::init(mpihelper);
    Dune::LoggingTimer<Dune::LoggingMessage>::init(mpihelper);

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, withGmshElementIndex, mpihelper);

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory, filename, mpihelper);

    // Create the grid
    return factory.createGrid();
}
