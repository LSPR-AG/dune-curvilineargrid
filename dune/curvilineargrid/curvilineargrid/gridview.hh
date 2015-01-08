// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GRIDVIEW_HH
#define DUNE_CURVGRID_GRIDVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/curvilineargrid/curvilineargrid/indexsets.hh>
#include <dune/curvilineargrid/curvilineargrid/intersection.hh>
#include <dune/curvilineargrid/curvilineargrid/intersectioniterator.hh>
#include <dune/curvilineargrid/curvilineargrid/iterator.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid, PartitionIteratorType pitype >
    class GridView;



    // GridViewTraits
    // --------------

    template< class Grid, PartitionIteratorType pitype >
    class GridViewTraits
    {
      friend class GridView< Grid, pitype >;

    public:
      typedef GridView< Grid, pitype > GridViewImp;

      typedef CurvGrid::IndexSet< const Grid > IndexSet;

      typedef Dune::Intersection< const Grid, CurvGrid::Intersection< const Grid > > Intersection;

      typedef Dune::IntersectionIterator < const Grid, CurvGrid::IntersectionIterator< const Grid>, CurvGrid::Intersection< const Grid> >
      IntersectionIterator;

      typedef typename Dune::CollectiveCommunication<MPI_Comm> CollectiveCommunication;

      template< int codim >
      struct Codim
      {
        typedef CurvGrid::IteratorTraits< codim, pitype, const Grid > IteratorTraits;
        typedef Dune::EntityIterator< codim, const Grid, CurvGrid::Iterator< IteratorTraits > > Iterator;

        typedef typename Grid::Traits::template Codim< codim >::Entity Entity;
        typedef typename Grid::Traits::template Codim< codim >::EntityPointer EntityPointer;

        typedef typename Grid::template Codim< codim >::Geometry Geometry;
        typedef typename Grid::template Codim< codim >::LocalGeometry LocalGeometry;

        template< PartitionIteratorType pit >
        struct Partition
        {
          typedef CurvGrid::IteratorTraits< codim, pit, const Grid > IteratorTraits;
          typedef Dune::EntityIterator< codim, const Grid, CurvGrid::Iterator< IteratorTraits > > Iterator;
        };
      };

      static const bool conforming = true;
    };



    // GridView
    // --------

    template< class Grid, PartitionIteratorType pitype >
    class GridView
    {
      typedef GridView< Grid, pitype > This;

    public:
      typedef GridViewTraits< Grid, pitype > Traits;

      typedef typename Traits::Grid Grid;

      typedef typename Traits::IndexSet IndexSet;

      typedef typename Traits::Intersection Intersection;

      typedef typename Traits::IntersectionIterator IntersectionIterator;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< int codim >
      struct Codim
        : public Traits::template Codim< codim >
      {};

      static const bool conforming = Traits::conforming;

      GridView ( const Grid &grid) : grid_( &grid )  { }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      const IndexSet &indexSet () const
      {
        if( !indexSet_ )  { indexSet_ = IndexSet( ); }
        return indexSet_;
      }


      //
      int size ( int codim ) const  { return hostGridView().size( codim ); }
      int size ( const GeometryType &type ) const  { return hostGridView().size( type ); }

      template< int codim >
      typename Codim< codim >::Iterator begin () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorTraits IteratorTraits;
        return CurvGrid::Iterator< IteratorTraits >( grid(), IteratorTraits::begin );
      }

      template< int codim, PartitionIteratorType pit >
      typename Codim< codim >::template Partition< pit >::Iterator begin () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorTraits IteratorTraits;
        return CurvGrid::Iterator< IteratorTraits >( grid(), IteratorTraits::begin );
      }

      template< int codim >
      typename Codim< codim >::Iterator end () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorTraits IteratorTraits;
        return CurvGrid::Iterator< IteratorTraits >( grid(), IteratorTraits::end );
      }

      template< int codim, PartitionIteratorType pit >
      typename Codim< codim >::template Partition< pit >::Iterator end () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorTraits IteratorTraits;
        return CurvGrid::Iterator< IteratorTraits >( grid(), IteratorTraits::end );
      }

      IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        typedef CurvGrid::IntersectionIterator< const Grid > IntersectionIteratorImpl;
        return IntersectionIteratorImpl(entity);
      }

      IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        typedef CurvGrid::IntersectionIterator< const Grid > IntersectionIteratorImpl;
        return IntersectionIteratorImpl(entity);
      }

      const CollectiveCommunication &comm () const
      {
        return grid().comm();
      }

      int overlapSize ( int codim ) const  { return grid().overlapSize(level_, codim); }

      int ghostSize ( int codim ) const  { return grid().ghostSize(level_, codim); }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction ) const
      {
    	  //[TODO] Call communication.hh
      }

    private:
      const Grid *grid_;
      mutable IndexSet indexSet_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRIDVIEW_HH
