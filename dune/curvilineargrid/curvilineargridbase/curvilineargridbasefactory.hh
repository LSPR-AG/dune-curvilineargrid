// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CURVGRID_GRIDBASEFACTORY_HH
#define DUNE_CURVGRID_GRIDBASEFACTORY_HH

/** \file
 *  \author Aleksejs Fomins
 *  \brief  Implementation of Curvilinear Grid Base Factory
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

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>




namespace Dune
{




template <int dim, int dimworld, class ct, bool isCached>
class CurvilinearFakeGrid
{
public:

	typedef  ct            ctype;
	static const int dimension = dim;
	static const int dimensionworld = dimworld;
	static const int is_cached = isCached;


	typedef Dune::CurvilinearGridStorage<ct, dim, isCached>  GridStorageType;

	CurvilinearFakeGrid() {}
};






template< class GridType >
class CurvilinearGridBaseFactory
{
  private:

	// Typedefs and const variables
	// -----------------------------------------------------------
	typedef typename GridType::ctype      ctype;

	static const bool  iscached     = GridType::is_cached;
    static const int dimension      = GridType::dimension;
    static const int dimensionworld = GridType::dimensionworld;

	typedef FieldVector< ctype, dimensionworld >                           VertexCoordinate;
	typedef Dune::CurvilinearGridBase<ctype, dimensionworld, iscached>     GridBaseType;

	typedef int VertexGlobalId;
	typedef int VertexLocalIndex;
	typedef int ElementLocalIndex;

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

    CurvilinearGridBaseFactory(
    		bool withGhostElements,
    		bool verbose,
    		bool processVerbose,
    		MPIHelper &mpihelper) :
    			verbose_(verbose),
    			processVerbose_(processVerbose),
    			mpihelper_(mpihelper)
    {
    	gridbase_ = new GridBaseType(withGhostElements, verbose, processVerbose, mpihelper);

    	rank_ = mpihelper.rank();
    	size_ = mpihelper.size();
    }

    ~CurvilinearGridBaseFactory ()  {}

    void insertVertex ( const VertexCoordinate &pos, const VertexGlobalId globalId )
    {
    	gridbase_->insertVertex(pos, globalId);
    }

    void insertElement(
      GeometryType &geometry,
      const int globalId,          // Not actually used
      const std::vector< VertexLocalIndex > &vertexIndexSet,
      const int elemOrder,
      const int physicalTag)
    {
    	gridbase_->insertElement(geometry, globalId, vertexIndexSet, elemOrder, physicalTag);
    }

    void insertBoundarySegment(
        GeometryType &geometry,
        const int globalId,          // Not actually used
        const std::vector< VertexGlobalId > &vertexIndexSet,
        const int elemOrder,
        const ElementLocalIndex associatedElementIndex,
        const int physicalTag)
    {
    	gridbase_->insertBoundarySegment(geometry, globalId, associatedElementIndex, vertexIndexSet, elemOrder, physicalTag);
    }



    void insertNVertexTotal(int nVertexTotal)  { gridbase_->insertNVertexTotal(nVertexTotal); }

    void insertNElementTotal(int nElementTotal)  { gridbase_->insertNElementTotal(nElementTotal); }


    GridBaseType * createGrid()
    {
    	gridbase_->generateMesh();
    	return gridbase_;
    }

    // Variables
    // -----------------------------------------------------------
  private:

    GridBaseType * gridbase_;



  };

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDBASEFACTORY_HH
