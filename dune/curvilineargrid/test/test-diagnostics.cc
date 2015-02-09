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
#include <dune/curvilineargrid/utility/curvilineargriddiagnostic.hh>



using namespace Dune;


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





/**\brief Test program which visualizes the base functions on a dgf mesh to
 * a vtk file. */
int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);
    int rank=mpihelper.rank();
    int size=mpihelper.size();


    // Assemble the file name
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME_SPHERE32_ORD4;

    // typedef  Dune::ALUGrid<3,3,simplex,nonconforming> SimplexGridType;
    typedef Dune::CurvilinearFakeGrid<3,3,double>  SimplexGridType;

    bool insertBoundarySegment = true;
    bool withGhostElements = true;
    bool verbose = false;
    bool processVerbose = false;

    bool writeReaderVTKFile = false;

    /** \brief provide a grid factory object for a grid of the ALUGSimplexGrid<3,3> type */
    //Dune::GridFactory<ALUSimplexGridType> factory;
    Dune::CurvilinearGridBaseFactory<SimplexGridType> factory(withGhostElements, verbose, processVerbose, mpihelper);


    /** \brief open the GMSH formatted tetrahedral mesh file into a grid factory */
    std::vector<int> boundaryId2physicalEntity;
    std::vector<int> elementIndex2PhysicalEntity;

    int nVertexTotal;
    int nElementTotal;

    Dune::CurvilinearGmshReader< SimplexGridType >::read(factory,
                                                            filename,
                                                            mpihelper,
                                                            boundaryId2physicalEntity,
                                                            elementIndex2PhysicalEntity,
                                                            nVertexTotal,
                                                            nElementTotal,
                                                            verbose,
                                                            processVerbose,
                                                            writeReaderVTKFile,
                                                            insertBoundarySegment);

    Dune::CurvilinearGridBase<double, 3> & gridbase = factory.createGrid(nVertexTotal, nElementTotal);



    // Perform diagnostics tests on the constructed grid
    Dune::CurvilinearGridDiagnostic<double, 3> diagnostic(verbose, processVerbose, mpihelper, gridbase);

	bool VTK_WRITE_ELEMENTS = true;
	bool VTK_WRITE_GHOST_ELEMENTS = true;

	bool VTK_WRITE_INTERNAL_FACE = false;
	bool VTK_WRITE_DOMAIN_BOUNDARY_FACE = true;
	bool VTK_WRITE_PROCESS_BOUNDARY_FACE = true;
	bool VTK_WRITE_GHOST_FACE = false;

	int  VTK_CURV_DISCRETIZATION = 2;       // 2=linear, minimal allowed discretization
	bool VTK_INTERPOLATE_DISCRETIZATION = true;
	bool VTK_EXPLODE_ELEMENTS = false;


	/*
    diagnostic.vtkWriteMesh(
    	VTK_WRITE_ELEMENTS,
    	VTK_WRITE_GHOST_ELEMENTS,
    	VTK_WRITE_INTERNAL_FACE,
    	VTK_WRITE_DOMAIN_BOUNDARY_FACE,
    	VTK_WRITE_PROCESS_BOUNDARY_FACE,
    	VTK_WRITE_GHOST_FACE,
    	VTK_CURV_DISCRETIZATION,
    	VTK_INTERPOLATE_DISCRETIZATION,
    	VTK_EXPLODE_ELEMENTS);

*/
    diagnostic.vtkWriteOctree();

    diagnostic.runAnalyticTest("curvilinearMeshAnalyticTest.txt");


    /** \brief leave program peacefully */
    return(0);
}
