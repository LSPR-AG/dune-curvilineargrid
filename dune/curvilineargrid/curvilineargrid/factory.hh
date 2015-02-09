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

//#include <dune/grid/utility/globalindex.hh>

//#include <dune/alugrid/common/transformation.hh>
//#include <dune/alugrid/3d/alugrid.hh>
//#include <dune/alugrid/3d/gridfactory.hh>
//#include <dune/alugrid/3d/gridfactory.cc>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>




namespace Dune
{


template< class ctype, int cdim >
class CurvilinearGridFactory
{
  private:

	typedef FieldVector< ctype, cdim >                 VertexCoordinate;
	typedef Dune::CurvilinearGridBase<ctype, cdim>     GridBaseType;

	typedef typename GridBaseType::LocalIndexType      LocalIndexType;
	typedef typename GridBaseType::GlobalIndexType     GlobalIndexType;

    bool verbose_;
    bool processVerbose_;

    // Parallel implementation
    MPIHelper &mpihelper_;
    int rank_;
    int size_;

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;


  public:

    CurvilinearGridFactory(
    		bool withGhostElements,
    		bool verbose,
    		bool processVerbose,
    		MPIHelper &mpihelper) :
    			gridbase_(withGhostElements, verbose, processVerbose, mpihelper),
    			verbose_(verbose),
    			processVerbose_(processVerbose),
    			mpihelper_(mpihelper)
    {
    	rank_ = mpihelper.rank();
    	size_ = mpihelper.size();
    }

    ~CurvilinearGridFactory ()  {}

    void insertVertex ( const VertexCoordinate &pos, const GlobalIndexType globalId )
    {
    	gridbase_.insertVertex(pos, globalId);
    }

    void insertElement(
      GeometryType &geometry,
      const int globalId,          // Not actually used
      const std::vector< LocalIndexType > &vertexIndexSet,
      const int elemOrder,
      const int physicalTag)
    {
    	gridbase_.insertElement(geometry, globalId, vertexIndexSet, elemOrder, physicalTag);
    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const int globalId,          // Not actually used
        const std::vector< LocalIndexType > &vertexIndexSet,
        const int elemOrder,
        const LocalIndexType associatedElementIndex,
        const int physicalTag)
    {
    	gridbase_.insertBoundarySegment(geometry, globalId, associatedElementIndex, vertexIndexSet, elemOrder, physicalTag);
    }

    Dune::CurvilinearGrid<cdim, cdim, ctype> createGrid(int nVertexTotal, int nElementTotal)
    {
    	gridbase_.generateMesh(nVertexTotal, nElementTotal);
    	return Dune::CurvilinearGrid<cdim, cdim, ctype>(gridbase_, mpihelper_);
    }


    // Variables
    // -----------------------------------------------------------
  private:

    GridBaseType gridbase_;



  };

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDFACTORY_HH
