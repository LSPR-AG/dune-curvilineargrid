// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_ITERATOR_HH
#define DUNE_CURVGRID_ITERATOR_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/curvilineargrid/curvilineargrid/entitypointer.hh>

namespace Dune
{

  namespace CurvGrid
  {

    template< int codim, PartitionIteratorType pitype, class Grid >
    class CurvLevelIterator : public CurvEntityPointer< codim, Grid >
    {
        typedef CurvEntityPointer< codim, Grid > Base;
        typedef typename Base::GridBaseType      GridBaseType;
        typedef typename Base::IndexSetIterator  IndexSetIterator;

    public:
        CurvLevelIterator (IndexSetIterator & iter, GridBaseType & gridbase)
    	  : Base( iter, gridbase, pitype)
        {}

        void increment ()
        {
      	  // Access by reference the entity stored in entity pointer, and call its method next() to iterate the entity.
      	  this->dereference().next();
        }
    };


    // HierarchicIterator
    // Note that HierarchicIterator only available over elements (codim 0)
    // ------------------

    template<class Grid >
    class CurvHierarchicIterator  : public CurvEntityPointer< 0, Grid >
    {
        typedef CurvEntityPointer< 0, Grid > Base;
        typedef typename Base::GridBaseType      GridBaseType;
        typedef typename Base::IndexSetIterator  IndexSetIterator;

    public:
        CurvHierarchicIterator (IndexSetIterator & iter, GridBaseType & gridbase)
    	: Base( iter, gridbase, PartitionIteratorType::All_Partition)
      {}

      void increment ()
      {
    	  // Access by reference the entity stored in entity pointer, and call its method next() to iterate the entity.
    	  this->dereference().next();
      }
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_ITERATOR_HH
