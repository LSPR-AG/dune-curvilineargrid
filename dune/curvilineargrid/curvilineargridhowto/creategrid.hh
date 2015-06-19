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
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH  + GMSH_FILE_NAME[interpOrder - 1]; //"bullseye-rev-400.msh";//

    // Additional constants
    bool withGhostElements = true;      // to create Ghost elements

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, mpihelper);

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory, filename, mpihelper);

    // Create the grid
    return factory.createGrid();
}
