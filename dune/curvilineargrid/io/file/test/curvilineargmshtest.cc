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

/** *\brief include functionality that encapsulates MPI */
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/mcmgmapper.hh>

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













/** Layout function returns true for all elements with codimensions. */

//template<int dim>
//struct P0Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim; } };

//template<int dim>
//struct P1Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim-1; } };

//template<int dim>
//struct P2Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim-2; } };

//template<int dim>
//struct P3Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim-3; } };



/**\brief Test program which visualizes the base functions on a dgf mesh to
 * a vtk file. */
int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);
    int rank=mpihelper.rank();
    int size=mpihelper.size();


    // Assemble the file name
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME_SPHERE32_ORD1;


    //! Define datatypes for grid mappers and  iterators to access all elements of a specific geometry type.
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P0Layout> ElementMapper;
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P1Layout> FaceMapper;
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P2Layout> EdgeMapper;
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P3Layout> NodeMapper;



    // typedef  Dune::ALUGrid<3,3,simplex,nonconforming> SimplexGridType;
    typedef Dune::CurvilinearFakeGrid<3,3,double>  SimplexGridType;

    bool insertBoundarySegment = true;
    bool withGhostElements = true;
    bool verbose = true;
    bool processVerbose = true;

    bool writeReaderVTKFile = true;


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

    //factory.createGrid(nVertexTotal, nElementTotal);

    /** \brief leave program peacefully */
    return(0);
}
