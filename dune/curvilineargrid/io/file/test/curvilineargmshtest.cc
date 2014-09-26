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


// dune grid includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
// Old ALUGRID
#if HAVE_ALUGRID
#warning "Including old AluGrid from dune/grid"
#include <dune/grid/alugrid.hh>
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
#endif
// New ALUGRID
#ifdef HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/grid/onedgrid.hh>

// alberta related stuff
#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

// grape include
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#endif









/** \brief provide hades3d namespace */
using namespace Dune;


/************************* define specific grid file names *************************/

/** \brief straight sided grids **/

//const std::string    GMSH_FILE_NAME    =    "sphere32.msh";
//const std::string    GMSH_FILE_NAME    =    "sphere32ord2.msh";
//const std::string    GMSH_FILE_NAME    =    "sphere32ord3.msh";
//const std::string    GMSH_FILE_NAME    =    "sphere32ord4.msh";
const std::string    GMSH_FILE_NAME    =    "sphere32ord5.msh";





/** Layout function returns true for all elements with codimensions. */

template<int dim>
struct P0Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim; } };

template<int dim>
struct P1Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim-1; } };

template<int dim>
struct P2Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim-2; } };

template<int dim>
struct P3Layout { bool contains (Dune::GeometryType gt) { return gt.dim()==dim-3; } };



/**\brief Test program which visualizes the base functions on a dgf mesh to
 * a vtk file. */
int main(int argc, char** argv)
{

#ifdef HAVE_MPI
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);
    int rank=mpihelper.rank();
    int size=mpihelper.size();
#else
    int rank=0;
    int size=1;
#endif

    //! Define datatypes for grid mappers and  iterators to access all elements of a specific geometry type.
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P0Layout> ElementMapper;
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P1Layout> FaceMapper;
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P2Layout> EdgeMapper;
    //typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Dune::GridSelector::GridType,P3Layout> NodeMapper;


#ifdef HAVE_MPI
    typedef  ALUGrid<3,3,simplex,nonconforming, mpihelper> ALUSimplexGridType;
#else
    typedef  ALUGrid<3,3,simplex,nonconforming > ALUSimplexGridType;
#endif

    /** \brief provide a grid factory object for a grid of the ALUGSimplexGrid<3,3> type */
    Dune::GridFactory<ALUSimplexGridType> factory;



    /** \brief open the GMSH formatted curvilinear tetrahedral mesh file */


    /** \brief open the GMSH formatted tetrahedral mesh file into a grid factory */
    std::vector<int> boundaryId2physicalEntity;
    boundaryId2physicalEntity.clear();

    std::vector<int> elementIndex2PhysicalEntity;
    elementIndex2PhysicalEntity.clear();

    Dune::CurvilinearGmshReader< ALUSimplexGridType >::read(factory,
                                                        GMSH_FILE_NAME,
                                                        boundaryId2physicalEntity,
                                                        elementIndex2PhysicalEntity,
                                                        true,
                                                        true);


    /** \brief leave program peacefully */
    return(0);
}
