#include <dune/curvilineargrid/curvilineargrid/factory.hh>
#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>


template <class GridType>
GridType * createGrid(Dune::MPIHelper & mpihelper)
{
	// Obtain path to the mesh folder from a CMake constant
    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

    // A choice of 32 tetrahedra spheres with interpolation orders 1 to 5.
    const std::string GMSH_FILE_NAME[5] {"sphere32.msh", "sphere32ord2.msh", "sphere32ord3.msh", "sphere32ord4.msh", "sphere32ord5.msh"};

    // Choice of file name
    int interpOrder = 2;
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH  + "sphere2000ord3.msh"; //GMSH_FILE_NAME[interpOrder - 1]; //"bullseye-rev-400.msh";//

    // Additional constants
    bool withGhostElements    = true;      // to create Ghost elements
    bool withGmshElementIndex = true;

    const bool LOGGING_MESSAGE_VERBOSE   = true;   // If the master process should report diagnostics
    const bool LOGGING_MESSAGE_PVERBOSE  = false;  // If all processes should report diagnostics (not recommended)
    const bool LOGGING_TIMER_REALVERBOSE = false;  // If LoggingTimer should report during timing (only for debug)

    // Initialize LoggingMessage and LoggingTimer
    Dune::LoggingMessage::getInstance().init(mpihelper, LOGGING_MESSAGE_VERBOSE, LOGGING_MESSAGE_PVERBOSE);
    Dune::LoggingTimer<Dune::LoggingMessage>::getInstance().init(LOGGING_TIMER_REALVERBOSE);

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, withGmshElementIndex, mpihelper);

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory, filename, mpihelper);

    // Create the grid
    return factory.createGrid();
}
