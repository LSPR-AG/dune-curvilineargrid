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

#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/factory.hh>
#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>
#include <dune/curvilineargrid/io/file/curvilinearvtkgridwriter.hh>
#include <dune/curvilineargrid/utility/griddiagnostic.hh>



using namespace Dune;

const bool isGeometryCached     = true;      // We will be using CachedCurvilinearGeometry class
const bool withGhostElements    = true;     // We want to create a mesh with ghost elements
const bool withGmshElementIndex = true;

const bool LOGGING_MESSAGE_VERBOSE   = true;   // If the master process should report diagnostics
const bool LOGGING_MESSAGE_PVERBOSE  = false;  // If all processes should report diagnostics (not recommended)
const bool LOGGING_TIMER_REALVERBOSE = false;  // If LoggingTimer should report during timing (only for debug)

/**\brief Test program which visualizes the base functions on a dgf mesh to
 * a vtk file. */
int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper = Dune::MPIHelper::instance(argc,argv);
    int rank = mpihelper.rank();
    int size = mpihelper.size();

    //
    /***************************************************************/
    /** Instantiation of the logging message and loggingtimer      */
    /***************************************************************/
    typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
    Dune::LoggingMessage::getInstance().init(mpihelper, LOGGING_MESSAGE_VERBOSE, LOGGING_MESSAGE_PVERBOSE);
    LoggingTimerDev::getInstance().init(LOGGING_TIMER_REALVERBOSE);

    /******************************************************/
    /** Define GridType and associated factory class      */
    /******************************************************/
    typedef Dune::CurvilinearGrid<double, 3, isGeometryCached, Dune::LoggingMessage> GridType;

    //Dune::GridFactory<ALUSimplexGridType> factory;
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, withGmshElementIndex, mpihelper);

    /******************************************************/
    /* Pass filename as command line argument             */
    /******************************************************/
    assert(argc > 1);
    std::string filename = std::string(argv[1]);

    /******************************************************/
    /* Read mesh and create grid                          */
    /******************************************************/
    Dune::CurvilinearGmshReader< GridType >::read(factory, filename, mpihelper);
    GridType * grid = factory.createGrid();

    /******************************************************/
    /* Perform diagnostics tests on the constructed grid  */
    /******************************************************/
    Dune::CurvilinearGridDiagnostic<GridType> diagnostic(mpihelper, *grid);
    diagnostic.runAnalyticTest("curvilinearMeshAnalyticTest.txt");

    /******************************************************/
    /* Write grid to VTK                                  */
    /******************************************************/
	Dune::CurvilinearVTKGridWriter<GridType> gridwriter(*grid);
	gridwriter.write("basis_test");


    /******************************************************/
    /* Write OCTree to VTK                                */
    /******************************************************/
    diagnostic.vtkWriteOctree();


    /******************************************************/
    /* Report timing information                          */
    /******************************************************/
    LoggingTimerDev::reportParallel(mpihelper);



    // Delete the GridBase
    delete grid;


    /** \brief leave program peacefully */
    return(0);
}
