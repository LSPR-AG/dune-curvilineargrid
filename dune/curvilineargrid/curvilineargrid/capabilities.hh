// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_CAPABILITIES_HH
#define DUNE_CURVGRID_CAPABILITIES_HH

#include <cassert>

//#include <dune/common/forloop.hh>

#include <dune/geometry/type.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
//#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune
{

  // Capabilities
  // ------------

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

  	// Note: At the moment curvilinear grid only capable of dealing with tetrahedral meshes
    template< class ctype, int cdim, bool isCached>
    struct hasSingleGeometryType< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
        static const bool v = true;
        //static const unsigned int topologyId = GenericGeometry::SimplexTopology<3>::type::id;
        // static const unsigned int topologyId = Impl::SimplexTopology< 3 >::type::id;
        static const unsigned int topologyId = GeometryTypes::simplex(3).id();
    };


    template< class ctype, int cdim, bool isCached, int codim >
    struct hasEntity< Dune::CurvilinearGrid< ctype, cdim, isCached>, codim >
    {
    	static const bool v = true;
    };


    /*
    template< class ctype, int cdim, bool isCached>
    struct isParallel< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
    	static const bool v = true;
    };
    */


    // FIXME: Do I need to specialize this for all codimensions to avoid, say, 4D grid requests?
    template< class ctype, int cdim, bool isCached, int codim >
    struct canCommunicate< Dune::CurvilinearGrid< ctype, cdim, isCached>, codim >
    {
    	static const bool v = true;
    };


    template< class ctype, int cdim, bool isCached>
    struct hasBackupRestoreFacilities< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
      static const bool v = false;
    };

    template< class ctype, int cdim, bool isCached>
    struct isLevelwiseConforming< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
      static const bool v = true;
    };

    template< class ctype, int cdim, bool isCached>
    struct isLeafwiseConforming< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
      static const bool v = true;
    };

    template< class ctype, int cdim, bool isCached>
    struct threadSafe< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
      static const bool v = false;
    };

    template< class ctype, int cdim, bool isCached>
    struct viewThreadSafe< Dune::CurvilinearGrid< ctype, cdim, isCached> >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_CAPABILITIES_HH
