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


    // GridFamily
    // ----------

    template <int dim, int dimworld, class ct>
    struct GridFamily
    {
      struct Traits
      {
        typedef CurvilinearGrid< dim, dimworld, ct> Grid;

        typedef typename ct ctype;

        static const int dimension       = dim;
        static const int dimensionworld  = dimworld;

        typedef CurvGrid::Intersection< const Grid, typename HostGrid::LeafIntersection >    BaseLeafIntersection;
        typedef CurvGrid::Intersection< const Grid, typename HostGrid::LevelIntersection >   BaseLevelIntersection;
        typedef CurvGrid::IntersectionIterator< const Grid, typename HostGrid::LeafIntersectionIterator >   BaseLeafIntersectionIterator;
        typedef CurvGrid::IntersectionIterator< const Grid, typename HostGrid::LevelIntersectionIterator >  BaseLevelIntersectionIterator;
        typedef CurvGrid::HierarchicIterator< const Grid >                                                  BaseHierarchicIterator;

        typedef Dune::Intersection< const Grid, BaseLeafIntersection > LeafIntersection;
        typedef Dune::Intersection< const Grid, BaseLevelIntersection > LevelIntersection;


        typedef Dune::IntersectionIterator< const Grid, BaseLeafIntersectionIterator, BaseLeafIntersection >    LeafIntersectionIterator;
        typedef Dune::IntersectionIterator< const Grid, BaseLevelIntersectionIterator, BaseLevelIntersection >  LevelIntersectionIterator;
        typedef Dune::EntityIterator< 0, const Grid, BaseHierarchicIterator>                                    HierarchicIterator;

        template< int codim >
        struct Codim
        {
          typedef Dune::CurvGrid::Geometry< dimension-codim, dimensionworld, const Grid > GeometryImpl;
          typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, GeometryImpl > Geometry;
          typedef typename HostGrid::template Codim< codim >::LocalGeometry LocalGeometry;

          typedef CurvGrid::EntityPointerTraits< codim, const Grid >     EntityPointerTraits;
          typedef CurvGrid::EntityPointer< EntityPointerTraits >         EntityPointerImpl;
          typedef Dune::EntityPointer< const Grid, EntityPointerImpl >   EntityPointer;
          typedef typename EntityPointerTraits::Entity                   Entity;

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

        typedef CurvGrid::IndexSet< const Grid, typename HostGrid::Traits::LeafIndexSet >   LeafIndexSet;
        typedef CurvGrid::IndexSet< const Grid, typename HostGrid::Traits::LevelIndexSet >  LevelIndexSet;

        typedef CurvGrid::IdSet< const Grid, typename HostGrid::Traits::GlobalIdSet >  GlobalIdSet;
        typedef CurvGrid::IdSet< const Grid, typename HostGrid::Traits::LocalIdSet >   LocalIdSet;

        typedef typename Dune::CollectiveCommunication<MPI_Comm> CollectiveCommunication;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::GridView< CurvGrid::GridViewTraits< typename HostGrid::LeafGridView, Allocator, pitype > >   LeafGridView;
          typedef Dune::GridView< CurvGrid::GridViewTraits< typename HostGrid::LevelGridView, Allocator, pitype > >  LevelGridView;
        };
      };
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRIDFAMILY_HH
