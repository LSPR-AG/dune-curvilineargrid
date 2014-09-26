// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GRIDFAMILY_HH
#define DUNE_CURVGRID_GRIDFAMILY_HH

#include <dune/grid/common/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/curvilineargrid/curvilineargrid/entity.hh>
#include <dune/curvilineargrid/curvilineargrid/entityseed.hh>
#include <dune/curvilineargrid/curvilineargrid/entitypointer.hh>
#include <dune/curvilineargrid/curvilineargrid/geometry.hh>
#include <dune/curvilineargrid/curvilineargrid/gridview.hh>
#include <dune/curvilineargrid/curvilineargrid/intersection.hh>
#include <dune/curvilineargrid/curvilineargrid/intersectioniterator.hh>
#include <dune/curvilineargrid/curvilineargrid/iterator.hh>
#include <dune/curvilineargrid/curvilineargrid/idset.hh>
#include <dune/curvilineargrid/curvilineargrid/indexsets.hh>

namespace Dune
{

  /** \brief namespace containing the implementations of CurvilinearGrid
   *  \ingroup CurvGrid
   */
  namespace CurvGrid
  {

    // ExportParams
    // ------------

    template< class HG, class CF >
    class ExportParams
    {
      static const bool isCoordFunction = isCoordFunctionInterface< typename CF::Interface >::value;
      static_assert(isCoordFunction, "Invalid CoordFunction.");

    public:
      typedef HG HostGrid;
      typedef CF CoordFunction;
    };



    // GridFamily
    // ----------

    template< class HG, class CF, class Allocator >
    struct GridFamily
    {
      struct Traits
      {
        typedef CurvilinearGrid< HG, CF, Allocator > Grid;

        typedef HG HostGrid;
        typedef CF CoordFunction;

        typedef typename HostGrid::ctype ctype;

        static const int dimension = HostGrid::dimension;
        static const int dimensionworld = CoordFunction::dimRange;

        typedef Dune::Intersection< const Grid, CurvGrid::Intersection< const Grid, typename HostGrid::LeafIntersection > > LeafIntersection;
        typedef Dune::Intersection< const Grid, CurvGrid::Intersection< const Grid, typename HostGrid::LevelIntersection > > LevelIntersection;

        typedef Dune::IntersectionIterator
        < const Grid, CurvGrid::IntersectionIterator< const Grid, typename HostGrid::LeafIntersectionIterator >, CurvGrid::Intersection< const Grid, typename HostGrid::LeafIntersection > >
        LeafIntersectionIterator;
        typedef Dune::IntersectionIterator
        < const Grid, CurvGrid::IntersectionIterator< const Grid, typename HostGrid::LevelIntersectionIterator >, CurvGrid::Intersection< const Grid, typename HostGrid::LevelIntersection > >
        LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const Grid, CurvGrid::HierarchicIterator< const Grid > >
        HierarchicIterator;

        template< int codim >
        struct Codim
        {
          typedef Dune::CurvGrid::Geometry< dimension-codim, dimensionworld, const Grid > GeometryImpl;
          typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, Dune::CurvGrid::Geometry > Geometry;
          typedef typename HostGrid::template Codim< codim >::LocalGeometry LocalGeometry;

          typedef CurvGrid::EntityPointerTraits< codim, const Grid > EntityPointerTraits;
          typedef CurvGrid::EntityPointer< EntityPointerTraits > EntityPointerImpl;
          typedef Dune::EntityPointer< const Grid, EntityPointerImpl > EntityPointer;
          typedef typename EntityPointerTraits::Entity Entity;

          typedef Dune::EntitySeed< const Grid, CurvGrid::EntitySeed< codim, const Grid > > EntitySeed;

          template< PartitionIteratorType pitype >
          struct Partition
          {
            typedef CurvGrid::IteratorTraits< typename HostGrid::LeafGridView, codim, pitype, const Grid > LeafIteratorTraits;
            typedef Dune::EntityIterator< codim, const Grid, CurvGrid::Iterator< LeafIteratorTraits > > LeafIterator;

            typedef CurvGrid::IteratorTraits< typename HostGrid::LevelGridView, codim, pitype, const Grid > LevelIteratorTraits;
            typedef Dune::EntityIterator< codim, const Grid, CurvGrid::Iterator< LevelIteratorTraits > > LevelIterator;
          };

          typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
          typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        };

        typedef CurvGrid::IndexSet< const Grid, typename HostGrid::Traits::LeafIndexSet > LeafIndexSet;
        typedef CurvGrid::IndexSet< const Grid, typename HostGrid::Traits::LevelIndexSet > LevelIndexSet;

        typedef CurvGrid::IdSet< const Grid, typename HostGrid::Traits::GlobalIdSet >
        GlobalIdSet;
        typedef CurvGrid::IdSet< const Grid, typename HostGrid::Traits::LocalIdSet >
        LocalIdSet;

        typedef typename HostGrid::Traits::CollectiveCommunication CollectiveCommunication;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::GridView< CurvGrid::GridViewTraits< typename HostGrid::LeafGridView, CoordFunction, Allocator, pitype > >
          LeafGridView;
          typedef Dune::GridView< CurvGrid::GridViewTraits< typename HostGrid::LevelGridView, CoordFunction, Allocator, pitype > >
          LevelGridView;
        };
      };
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRIDFAMILY_HH
