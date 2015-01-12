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


    // GridView
    // --------

    template< class Grid, PartitionIteratorType pitype >
    class GridView
    {
    	typedef Dune::CurvilinearGridBase<ct, dim>       GridBaseType;
    	typedef typename GridBaseType::IndexSetIterator  IndexSetIterator;

    public:

      static const bool conforming = Traits::conforming;

      GridView (const Grid &grid, GridBaseType & gridbase)
             : grid_( &grid ), gridbase_(gridbase), indexset_(gridbase)
      { }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }

      const IndexSet &indexSet () const
      {
        if( !indexset_ )  { indexset_ = IndexSet( ); }
        return indexset_;
      }


      // Get the number of entities within this gridview
      int size ( int codim )                 const  { return indexset_.size(codim); }
      int size ( const GeometryType &type )  const  { return indexset_.size(type); }

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

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction ) const
      {
    	  //[TODO] Call communication.hh
      }

    protected:
      const Grid *grid_;
      GridBaseType & gridbase_;
      mutable IndexSet indexset_;
    };






    template< class Grid, PartitionIteratorType pitype >
    class LeafGridView : GridView<Grid, pitype >
    {
    	typedef Dune::CurvilinearGridBase<ct, dim>       GridBaseType;
    	typedef typename GridBaseType::IndexSetIterator  IndexSetIterator;

    	typedef CurvGrid::LeafIterator< codim, Grid >    LeafIterator;

    public:
    	typedef GridView<Grid, pitype >  Base;

    	using Base::grid_;
    	using Base::gridbase_;
    	using Base::indexset_;

    	LeafGridView(const Grid &grid, GridBaseType & gridbase)
             : Base (grid, gridbase)
    	{ }

        template< int codim >
        typename Codim< codim >::Iterator begin () const
        {
        	return LeafIterator(gridbase_.entityDuneIndexBegin(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Codim< codim >::template Partition< pit >::Iterator begin () const
        {
        	return LeafIterator(gridbase_.entityDuneIndexBegin(codim, pit), gridbase_, grid());
        }

        template< int codim >
        typename Codim< codim >::Iterator end () const
        {
      	  return LeafIterator(gridbase_.entityDuneIndexEnd(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Codim< codim >::template Partition< pit >::Iterator end () const
        {
      	  return LeafIterator(gridbase_.entityDuneIndexEnd(codim, pit), gridbase_, grid());
        }

        int overlapSize ( int codim ) const  { return grid().overlapSize(codim); }

        int ghostSize ( int codim ) const  { return grid().ghostSize(codim); }

    };





    template< class Grid, PartitionIteratorType pitype >
    class LevelGridView : GridView<Grid, pitype >
    {
    	typedef Dune::CurvilinearGridBase<ct, dim>       GridBaseType;
    	typedef typename GridBaseType::IndexSetIterator  IndexSetIterator;

    	typedef CurvGrid::LevelIterator< codim, Grid >    LevelIterator;

    public:
    	typedef GridView<Grid, pitype >  Base;

    	using Base::grid_;
    	using Base::gridbase_;
    	using Base::indexset_;

    	LevelGridView(const Grid &grid, GridBaseType & gridbase)
             : Base (grid, gridbase),
               level_(level)
    	{ }

        template< int codim >
        typename Codim< codim >::Iterator begin () const
        {
        	return LevelIterator(gridbase_.entityDuneIndexBegin(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Codim< codim >::template Partition< pit >::Iterator begin () const
        {
        	return LevelIterator(gridbase_.entityDuneIndexBegin(codim, pit), gridbase_, grid());
        }

        template< int codim >
        typename Codim< codim >::Iterator end () const
        {
      	  return LevelIterator(gridbase_.entityDuneIndexEnd(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Codim< codim >::template Partition< pit >::Iterator end () const
        {
      	  return LevelIterator(gridbase_.entityDuneIndexEnd(codim, pit), gridbase_, grid());
        }

        int overlapSize ( int codim ) const  { return grid().overlapSize(level_, codim); }

        int ghostSize ( int codim ) const  { return grid().ghostSize(level_, codim); }

    private:
        int level_;

    };



  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRIDVIEW_HH
