// Autotool generated header.
#include <config.h>

// Stl headers.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

/** *\brief include functionality that encapsulates MPI */
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbasefactory.hh>
#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>
#include <dune/curvilineargrid/utility/gridbasediagnostic.hh>



using namespace Dune;

using namespace Dune::CurvGrid;


// Define path to meshes
const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

// Define mesh file names
const std::string    GMSH_FILE_NAME_SPHERE32_ORD1     =    "sphere32.msh";
const std::string    GMSH_FILE_NAME_SPHERE32_ORD2     =    "sphere32ord2.msh";
const std::string    GMSH_FILE_NAME_SPHERE32_ORD3     =    "sphere32ord3.msh";
const std::string    GMSH_FILE_NAME_SPHERE32_ORD4     =    "sphere32ord4.msh";
const std::string    GMSH_FILE_NAME_SPHERE32_ORD5     =    "sphere32ord5.msh";
const std::string    GMSH_FILE_NAME_SPHERE2000_ORD3   =    "sphere2000ord3.msh";
const std::string    GMSH_FILE_NAME_BULLSEYE400_ORD1  =    "bullseye-rev-400.msh";
const std::string    GMSH_FILE_NAME_SPHEREINSPHERE100_ORD1  =    "sphere-in-sphere-in-sphere-rev-100.msh";



const bool isGeometryCached = true;



/**\brief Test program which visualizes the base functions on a dgf mesh to
 * a vtk file. */
int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);

    // Instantiation of the logging message and loggingtimer
    typedef LoggingTimer<LoggingMessage>  LoggingTimer;
    LoggingMessage::init(mpihelper);
    LoggingTimer::init(mpihelper);

    // typedef  Dune::ALUGrid<3,3,simplex,nonconforming> SimplexGridType;
    typedef CurvilinearGridBase<double, 3, isGeometryCached>  SimplexGridType;

    bool withGhostElements = true;
    bool withGmshElementIndex = true;

    /** \brief provide a grid factory object for a grid of the ALUGSimplexGrid<3,3> type */
    //Dune::GridFactory<ALUSimplexGridType> factory;
    CurvilinearGridBaseFactory<SimplexGridType> factory(withGhostElements, withGmshElementIndex, mpihelper);


    // Assemble the file name
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME_SPHERE2000_ORD3;
    CurvilinearGmshReader< SimplexGridType >::read(factory, filename, mpihelper);

    SimplexGridType * gridbase = factory.createGrid();



    // Perform diagnostics tests on the constructed grid
    CurvilinearGridBaseDiagnostic<SimplexGridType> diagnostic(mpihelper, *gridbase);

	std::vector<bool> withElements {true, true};                 // Whether to add elements to VTK: Internal / Ghost
	std::vector<bool> withFaces    {false, false, true, true};   // Whether to add faces to VTK: Internal / Ghost / ProcessBoundary / DomainBoundary
	std::vector<bool> withEdges    {false, false, false};        // Whether to add edges to VTK: Internal / Ghost / ProcessBoundary

	int  VTK_CURV_DISCRETIZATION = 7;       // 2=linear, minimal allowed discretization
	bool VTK_INTERPOLATE_DISCRETIZATION = true;
	bool VTK_EXPLODE_ELEMENTS = true;
	std::vector<bool> VTK_WRITE_CODIM {true, true, false, false};  // Use tetrahedra and triangles to discretize the mesh

    diagnostic.vtkWriteMesh(
    	withElements,
    	withFaces,
    	withEdges,
    	VTK_CURV_DISCRETIZATION,
    	VTK_INTERPOLATE_DISCRETIZATION,
    	VTK_EXPLODE_ELEMENTS,
    	VTK_WRITE_CODIM
    );

    diagnostic.vtkWriteOctree();

    diagnostic.runAnalyticTest("curvilinearMeshAnalyticTest.txt");

    LoggingTimer::reportParallel();


    // Delete the GridBase
    delete gridbase;


    /** \brief leave program peacefully */
    return(0);
}
