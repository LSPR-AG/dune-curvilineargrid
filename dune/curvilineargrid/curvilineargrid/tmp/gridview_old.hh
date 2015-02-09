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
    	typedef typename Traits::ctype ctype;

    	typedef typename Traits::GridStorageType     GridStorageType;
    	typedef typename Traits::GridBaseType        GridBaseType;
    	typedef typename Traits::LocalIndexType      LocalIndexType;
    	typedef typename Traits::InternalIndexType   InternalIndexType;
    	typedef typename Traits::IndexSetIterator    IndexSetIterator;

        // Codimensions of entity types for better code readability
        static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
        static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
        static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
        static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    	typedef typename Traits::template Codim< ELEMENT_CODIM >::Entity  Entity;

    	typedef typename Traits::CollectiveCommunication      CollectiveCommunication;

    	typedef typename Traits::IntersectionIterator         IntersectionIterator;
    	typedef typename Traits::IntersectionIteratorImpl     IntersectionIteratorImpl;

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
    	  InternalIndexType firstFaceSubIndex = 0;
    	  LocalIndexType elementLocalIndex = entity.localIndex();
    	  LocalIndexType faceLocalIndex = gridbase_.subentityLocalIndex(elementLocalIndex, ELEMENT_CODIM, FACE_CODIM, firstFaceSubIndex);

    	  IntersectionIteratorImpl iter (elementLocalIndex, firstFaceSubIndex, gridbase_);

          // Iterator must not point at a Ghost face
          // If it does, increment it, it will automatically point at the next non-ghost face
          if (gridbase_.entityStructuralType(FACE_CODIM, faceLocalIndex) == GridStorageType::PartitionType::Ghost)  { iter.increment(); }

          return iter;
      }


      IntersectionIterator iend ( const Entity &entity ) const
      {
    	  // 1) Create IntersectionIterator with subentityIndex = subentityNumber;
    	  InternalIndexType nSubentityFace = Dune::ReferenceElements<ctype, dimension>::general(entity.type()).size(FACE_CODIM);
    	  return IntersectionIteratorImpl(entity.localIndex(), nSubentityFace, gridbase_);
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

      	if (dataHandle.contains(dimension, ELEMENT_CODIM)) { communicator.communicateWrapper<DataHandle, ELEMENT_CODIM>(dataHandle, interface, direction, level); }
      	if (dataHandle.contains(dimension, FACE_CODIM))    { communicator.communicateWrapper<DataHandle, FACE_CODIM>(dataHandle, interface, direction, level); }
      	if (dataHandle.contains(dimension, EDGE_CODIM))    { communicator.communicateWrapper<DataHandle, EDGE_CODIM>(dataHandle, interface, direction, level); }
      	if (dataHandle.contains(dimension, VERTEX_CODIM))  { communicator.communicateWrapper<DataHandle, VERTEX_CODIM>(dataHandle, interface, direction, level); }
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
             : Base (grid, gridbase),
               indexset_(gridbase)
    	{ }

        const IndexSet &indexSet () const  { return indexset_; }

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
               indexset_(gridbase),
               level_(level)
    	{ }


        const IndexSet &indexSet () const  { return indexset_; }


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
