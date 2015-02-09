// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_CAPABILITIES_HH
#define DUNE_CURVGRID_CAPABILITIES_HH

#include <cassert>

#include <dune/common/forloop.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  // Capabilities
  // ------------

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

  	// Note: At the moment curvilinear grid only capable of dealing with tetrahedral meshes
    template< int dim, int dimworld, class ctype>
    struct hasSingleGeometryType< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
        static const bool v = true;
        static const unsigned int topologyId = GenericGeometry::SimplexTopology<3>::type::id;
    };


    template< int dim, int dimworld, class ctype, int codim >
    struct hasEntity< Dune::CurvilinearGrid< dim , dimworld, ctype>, codim >
    {
    	static const bool v = true;
    };


    template< int dim, int dimworld, class ctype>
    struct isParallel< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
    	static const bool v = true;
    };


    // FIXME: Do I need to specialize this for all codimensions to avoid, say, 4D grid requests?
    template< int dim, int dimworld, class ctype, int codim >
    struct canCommunicate< Dune::CurvilinearGrid< dim , dimworld, ctype>, codim >
    {
    	static const bool v = true;
    };


    template< int dim, int dimworld, class ctype>
    struct hasBackupRestoreFacilities< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class ctype>
    struct isLevelwiseConforming< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class ctype>
    struct isLeafwiseConforming< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class ctype>
    struct threadSafe< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class ctype>
    struct viewThreadSafe< Dune::CurvilinearGrid< dim , dimworld, ctype> >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_CAPABILITIES_HH
