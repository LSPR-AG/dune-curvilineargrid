// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CURVGRID_GRIDFACTORY_HH
#define DUNE_CURVGRID_GRIDFACTORY_HH

/** \file
 *  \author Aleksejs Fomins
 *  \brief  Implementation of Curvilinear Grid Factory
 */

#include <config.h>

#include <map>
#include <limits>
#include <vector>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundaryprojection.hh>

//#include <dune/grid/utility/globalindex.hh>

#include <dune/alugrid/common/transformation.hh>
#include <dune/alugrid/3d/alugrid.hh>
#include <dune/alugrid/3d/gridfactory.hh>
#include <dune/alugrid/3d/gridfactory.cc>



#include <parmetis.h>






namespace Dune
{


template< class ALUGrid >
class CurvilinearGridFactory
{
  private:

	typedef double ctype;

    static const unsigned int dimension = 3;
    static const unsigned int dimensionworld = 3;

	typedef FieldVector< ctype, dimensionworld > VertexCoordinate;

	typedef unsigned int VertexGlobalId;
	typedef unsigned int VertexLocalIndex;
	typedef unsigned int ElementLocalIndex;


    // Parallel Implementation
    MPIHelper &mpihelper_;
    static const int MASTER_RANK = 0;

    bool verbose_ = true;

  public:

    CurvilinearGridFactory(MPIHelper &mpihelper) : mpihelper_(mpihelper) {}

    ~CurvilinearGridFactory ()  {}

    void insertVertex ( const VertexCoordinate &pos, const VertexGlobalId globalId )   { }

    void insertElement(
      GeometryType &geometry,
      const int globalId,
      const std::vector< VertexLocalIndex > &elementVertexId,
      const int elemOrder)
    {

    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const int globalId,
        const std::vector< VertexGlobalId > &boundaryVertexId,
        const int elemOrder,
        const ElementLocalIndex linkedElement)
    {

    }

    void createGrid()  { }


  private:

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::string tmp_str;
        for (int i = 0; i < V.size(); i++) { tmp_str += std::to_string(V[i]) + " "; }
        return tmp_str;
    }

    // Writes debug info to the command line
    // TODO: Use IFDEF to manipulate between no output, all output, or only master process output
    void print_debug(std::string s)
    {
        if (verbose_) { std::cout << "Process_" << mpihelper_.rank() << ": " << s << std::endl; }
    }



  };

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDFACTORY_HH
