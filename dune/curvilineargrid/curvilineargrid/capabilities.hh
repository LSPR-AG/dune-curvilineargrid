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
    template< int dim, int dimworld>
    struct hasSingleGeometryType< CurvilinearGrid< dim , dimworld> >
    {
        static const bool v = true;
        static const unsigned int topologyId = GenericGeometry::SimplexTopology<3>::type::id;
    };


    template< int dim, int dimworld, int codim >
    struct hasEntity< CurvilinearGrid< dim , dimworld>, codim >
    {
    	static const bool v = true;
    };


    template< int dim, int dimworld>
    struct isParallel< CurvilinearGrid< dim , dimworld> >
    {
    	static const bool v = true;
    };


    // FIXME: Do I need to specialize this for all codimensions to avoid, say, 4D grid requests?
    template< int dim, int dimworld, int codim >
    struct canCommunicate< CurvilinearGrid< dim , dimworld>, codim >
    {
    	static const bool v = true;
    };


    template< int dim, int dimworld>
    struct hasBackupRestoreFacilities< CurvilinearGrid< dim , dimworld> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld>
    struct isLevelwiseConforming< CurvilinearGrid< dim , dimworld> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld>
    struct isLeafwiseConforming< CurvilinearGrid< dim , dimworld> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld>
    struct threadSafe< CurvilinearGrid< dim , dimworld> >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld>
    struct viewThreadSafe< CurvilinearGrid< dim , dimworld> >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_CAPABILITIES_HH
