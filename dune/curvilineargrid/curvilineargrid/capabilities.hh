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
    template< class HostGrid, class Allocator >
    struct hasSingleGeometryType< CurvilinearGrid< HostGrid, Allocator > >
    {
        static const bool v = true;
        static const unsigned int topologyId = GenericGeometry::SimplexTopology<3>::type::id;
    };


    template< class HostGrid, class Allocator, int codim >
    struct hasEntity< CurvilinearGrid< HostGrid, Allocator >, codim >
    {
    	static const bool v = true;
    };


    template< class HostGrid, class Allocator >
    struct isParallel< CurvilinearGrid< HostGrid, Allocator > >
    {
    	static const bool v = ???;
    };


    // FIXME: Do I need to specialize this for all codimensions to avoid, say, 4D grid requests?
    template< class HostGrid, class Allocator, int codim >
    struct canCommunicate< CurvilinearGrid< HostGrid, Allocator >, codim >
    {
    	static const bool v = true;
    };


    template< class HostGrid, class Allocator >
    struct hasBackupRestoreFacilities< CurvilinearGrid< HostGrid, Allocator > >
    {
      static const bool v = false;
    };

    template< class HostGrid, class Allocator >
    struct isLevelwiseConforming< CurvilinearGrid< HostGrid, Allocator > >
    {
      static const bool v = isLevelwiseConforming< HostGrid >::v;
    };

    template< class HostGrid, class Allocator >
    struct isLeafwiseConforming< CurvilinearGrid< HostGrid, Allocator > >
    {
      static const bool v = isLeafwiseConforming< HostGrid >::v;
    };

    template< class HostGrid, class Allocator >
    struct threadSafe< CurvilinearGrid< HostGrid, Allocator > >
    {
      static const bool v = false;
    };

    template< class HostGrid, class Allocator >
    struct viewThreadSafe< CurvilinearGrid< HostGrid, Allocator > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_CAPABILITIES_HH
