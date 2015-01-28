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

#include <dune/curvilineargrid/curvilineargrid/communication.hh>

namespace Dune
{

  namespace CurvGrid
  {


    // GridView
    // --------

    template< class Grid, PartitionIteratorType pitype >
    class GridView
    {
    	typedef typename Grid::Traits Traits;

    	typedef typename Traits::GridStorageType     GridStorageType;
    	typedef typename Traits::GridBaseType        GridBaseType;
    	typedef typename Traits::LocalIndexType      LocalIndexType;
    	typedef typename Traits::IndexSetIterator    IndexSetIterator;



    	typedef typename Traits::template Codim< 0 >::Entity  Entity;

    	typedef typename Traits::CollectiveCommunication      CollectiveCommunication;

    	typedef typename Traits::IntersectionIterator         IntersectionIterator;
    	typedef typename Traits::IntersectionIteratorImpl     IntersectionIteratorImpl;

    	typedef Dune::IteratorRange<>  GridviewIteratorRange;

    public:

      static const int dimension =  Traits::dimension;
      static const bool conforming = Traits::conforming;

      GridView (const Grid &grid, GridBaseType & gridbase)
             : grid_( &grid ), gridbase_(gridbase)
      { }

      const Grid &grid () const
      {
        assert( grid_ );
        return *grid_;
      }


	  // 1) Create IntersectionIterator with subentityIndex = 0;
	  // 2) Check the type of associated face
	  // 3) if (faceType == Ghost) then iterator++
      IntersectionIterator ibegin ( const Entity &entity ) const
      {
    	  LocalIndexType elementLocalIndex = entity.localIndex();
    	  LocalIndexType faceLocalIndex = gridbase_.subentityLocalIndex(elementLocalIndex, 0, 1, 0);

    	  IntersectionIteratorImpl iter (elementLocalIndex, 0, gridbase_);

          // Iterator must not point at a Ghost face
          // If it does, increment it, it will automatically point at the next non-ghost face
          if (gridbase_.entityStructuralType(1, faceLocalIndex) == GridStorageType::PartitionType::Ghost)  { iter.increment(); }

          return iter;
      }


      // FIXME: Replace the number 4 with SubentitySize
      IntersectionIterator iend ( const Entity &entity ) const
      {
    	  // 1) Create IntersectionIterator with subentityIndex = subentityNumber;
    	  return IntersectionIteratorImpl(entity.localIndex(), 4, gridbase_);
      }

      const CollectiveCommunication &comm () const
      {
        return grid().comm();
      }


  	  // Wrapper Communication Algorithm
  	  // 1) Loop over all codim
  	  // 1.1) Check if this codim is allowed by DataHandle
  	  // 1.2) For each InterfaceSubset, check if it is consistent with InterfaceType and CommunicationDirection
  	  // 1.3) If it is, call main communication protocol main_communicate(codim, mapSend, ranklistSend)
      template< class DataHandle, class Data, class GridViewType >
      void communicateCodim( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction,
                         GridViewType & gv ) const
      {
    	Dune::CurvGrid::Communication<Grid> communicator;
    	int level = 0; // Fake

      	if (dataHandle.contains(dimension, 0))  { communicator.communicateWrapper<DataHandle, 0>(dataHandle, interface, direction, level); }
      	if (dataHandle.contains(dimension, 1))  { communicator.communicateWrapper<DataHandle, 1>(dataHandle, interface, direction, level); }
      	if (dataHandle.contains(dimension, 2))  { communicator.communicateWrapper<DataHandle, 2>(dataHandle, interface, direction, level); }
      	if (dataHandle.contains(dimension, 3))  { communicator.communicateWrapper<DataHandle, 3>(dataHandle, interface, direction, level); }
      }

    protected:
      const Grid *grid_;
      GridBaseType & gridbase_;
    };




    template< class Grid, PartitionIteratorType pitype >
    class LeafGridView : GridView<Grid, pitype >
    {
    	typedef typename Grid::Traits Traits;

    	typedef typename Traits::GridBaseType      GridBaseType;
    	typedef typename Traits::IndexSetIterator  IndexSetIterator;

    	typedef typename Traits::LeafIndexSet  LeafIndexSet;


    public:
    	typedef LeafGridView<Grid, pitype > This;
    	typedef GridView<Grid, pitype >  Base;

    	using Base::grid_;
    	using Base::gridbase_;

    	LeafGridView(const Grid &grid, GridBaseType & gridbase)
             : Base (grid, gridbase)
    	{ }

        const IndexSet &indexSet () const
        {
          if( !indexset_ )  { indexset_ = LeafIndexSet(gridbase); }
          return indexset_;
        }


        // Get the number of entities within this gridview
        int size ( int codim )                       const  { return indexset_.size(codim); }
        int size ( const Dune::GeometryType &type )  const  { return indexset_.size(type); }


        template< int codim >
        typename Traits::template Codim< codim >::LeafIterator begin () const
        {
        	return LeafIterator(gridbase_.entityDuneIndexBegin(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Traits::template Codim< codim >::template Partition< pit >::LeafIterator begin () const
        {
        	return LeafIterator(gridbase_.entityDuneIndexBegin(codim, pit), gridbase_, grid());
        }

        template< int codim >
        typename Traits::template Codim< codim >::LeafIterator end () const
        {
      	  return LeafIterator(gridbase_.entityDuneIndexEnd(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Traits::template Codim< codim >::template Partition< pit >::LeafIterator end () const
        {
      	  return LeafIterator(gridbase_.entityDuneIndexEnd(codim, pit), gridbase_, grid());
        }

        int overlapSize ( int codim ) const  { return grid().overlapSize(codim); }

        int ghostSize ( int codim ) const  { return grid().ghostSize(codim); }




        template< class DataHandle, class Data>
        void communicate( CommDataHandleIF< DataHandle, Data > &dataHandle,
                           InterfaceType interface,
                           CommunicationDirection direction ) const
        {
        	Base::communicateCodim< DataHandle, Data, This>(dataHandle, interface, direction, *this);
        }


        template <partitions>
        Dune::PartitionSet< partitions >
        interface2partitionSet(InterfaceType interface)
        {
        	switch(interface)
        	{
            InteriorBorder_InteriorBorder_Interface=0,     //!< send/receive interior and border entities
            InteriorBorder_All_Interface=1,                //!< send interior and border, receive all entities
            All_All_Interface=4                            //!< send all and receive all entities
        	}

        }


    private:
        mutable LeafIndexSet indexset_;

    };




    template< class Grid, PartitionIteratorType pitype >
    class LevelGridView : GridView<Grid, pitype >
    {
    	typedef typename Grid::Traits Traits;

    	typedef typename Traits::GridBaseType      GridBaseType;
    	typedef typename Traits::IndexSetIterator  IndexSetIterator;

    	typedef typename Traits::LevelIndexSet  LevelIndexSet;

    public:
    	typedef GridView<Grid, pitype >  Base;

    	using Base::grid_;
    	using Base::gridbase_;

    	LevelGridView(const Grid &grid, GridBaseType & gridbase, int level)
             : Base (grid, gridbase),
               level_(level)
    	{ }


        const IndexSet &indexSet () const
        {
          if( !indexset_ )  { indexset_ = LevelIndexSet(gridbase); }
          return indexset_;
        }


        // Get the number of entities within this gridview
        int size ( int codim )                       const  { return indexset_.size(codim); }
        int size ( const Dune::GeometryType &type )  const  { return indexset_.size(type); }


        template< int codim >
        typename Traits::template Codim< codim >::LevelIterator begin () const
        {
        	return LevelIterator(gridbase_.entityDuneIndexBegin(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Traits::template Codim< codim >::template Partition< pit >::LevelIterator begin () const
        {
        	return LevelIterator(gridbase_.entityDuneIndexBegin(codim, pit), gridbase_, grid());
        }

        template< int codim >
        typename Traits::template Codim< codim >::LevelIterator end () const
        {
      	  return LevelIterator(gridbase_.entityDuneIndexEnd(codim, pitype), gridbase_, grid());
        }

        template< int codim, PartitionIteratorType pit >
        typename Traits::template Codim< codim >::template Partition< pit >::LevelIterator end () const
        {
      	  return LevelIterator(gridbase_.entityDuneIndexEnd(codim, pit), gridbase_, grid());
        }

        int overlapSize ( int codim ) const  { return grid().overlapSize(level_, codim); }

        int ghostSize ( int codim ) const  { return grid().ghostSize(level_, codim); }

    private:
        int level_;
        mutable LevelIndexSet indexset_;

    };



  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRIDVIEW_HH
