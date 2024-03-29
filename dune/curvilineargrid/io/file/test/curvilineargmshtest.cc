/** \file
 *  \brief Implementation of a test program that studies how mesh files in
 *         the GMSH file format can be read into Dune based codes, including
 *         boundary id's and elemental tags.
 *
 *
 *  Copyright by Benedikt Oswald and Patrick Leidenberger, 2002-2009.
 *  All rights reserved.
 *
 *  Objective: read a finite element mesh in the GMSH format from a GMSH file,
 *             assign boundary id's and elemental tags.
 *
 *  \author    Benedikt Oswald, Patrick Leidenberger
 *  \date      2009 oct 10, created, benedikt oswald
 *  \date      2009 oct 10, adapted, benedikt oswald, added gmsh reader support.
 *
 *  \warning   None.
 *  \attention Attention to nothing.
 *  \bug
 *  \todo
 */

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
#include <assert.h>

/** *\brief include functionality that encapsulates MPI */
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

//#include <dune/grid/common/mcmgmapper.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbasefactory.hh>
#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>



// Old ALUGRID
// #if HAVE_ALUGRID
// #warning "Including old AluGrid from dune/grid"
// #include <dune/grid/alugrid.hh>
// #include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
// #endif
// New ALUGRID
// #ifdef HAVE_DUNE_ALUGRID
// #include <dune/alugrid/grid.hh>
// #endif

/** \brief provide hades3d namespace */
using namespace Dune;

using namespace Dune::CurvGrid;


/**\brief Test program which visualizes the base functions on a dgf mesh to
 * a vtk file. */
int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);
    int rank=mpihelper.rank();
    int size=mpihelper.size();


    // Assemble the file name
    assert(argc > 1);  // User must provide file name
    std::string filename(argv[1]);

    // Properties of the grid
    bool withGhostElements    = true;
    bool withGmshElementIndex = true;
    bool withProcessPartition = true;
    CurvilinearGmshReaderLoadBalanceStrategy LBstrat = LoadBalanceDefault;
    //CurvilinearGmshReaderLoadBalanceStrategy LBstrat = LoadBalanceBoundary;
    bool writeCurvReaderVTK = true;

    const bool isCached = true;

    // Instantiation of the logging message
    typedef LoggingTimer<LoggingMessage>  LoggingTimer;
    LoggingMessage::init(mpihelper);
    LoggingTimer::init(mpihelper);

    // [TODO] DEBUG
    //LoggingMessage::setPVerbose(true);

    // typedef  Dune::ALUGrid<3,3,simplex,nonconforming> SimplexGridType;
    typedef CurvilinearGridBase<double, 3, isCached>  SimplexGridType;

    /** \brief provide a grid factory object for a grid of the ALUGSimplexGrid<3,3> type */
    //Dune::GridFactory<ALUSimplexGridType> factory;
    CurvilinearGridBaseFactory<SimplexGridType> factory(withGhostElements, withGmshElementIndex, mpihelper);

    CurvilinearGmshReader< SimplexGridType>::read(factory, filename, mpihelper, withGmshElementIndex, withProcessPartition, LBstrat, writeCurvReaderVTK);

    //factory.createGrid(nVertexTotal, nElementTotal);

    //LoggingTimer::report();
    LoggingTimer::reportParallel();

    /** \brief leave program peacefully */
    return(0);
}
