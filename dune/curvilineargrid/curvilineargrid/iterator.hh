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

    // Internal Forward Declarations
    // -----------------------------

    template< int codim, class Grid >
    class LeafIterator : public EntityPointer< codim, Grid >
    {
        typedef EntityPointer< codim, Grid > Base;
        typedef typename Grid::Traits Traits;

    protected:
        typedef typename Base::EntityImpl EntityImpl;
        typedef typename Base::EntitySeed EntitySeed;

    public:
        LeafIterator ( const EntitySeed & seed, const Grid &grid)
         : Base( seed, grid)
        {}

        void increment ()
        {
      	  // Access by reference the entity stored in entity pointer, and call its method next() to iterate the entity.
      	  this->dereference().next();
        }
    };


    template< int codim, class Grid >
    class LevelIterator : public EntityPointer< codim, Grid >
    {
        typedef EntityPointer< codim, Grid > Base;
        typedef typename Grid::Traits Traits;

    protected:
        typedef typename Base::EntityImpl EntityImpl;
        typedef typename Base::EntitySeed EntitySeed;

    public:
        LevelIterator ( const EntitySeed & seed, const Grid &grid)
         : Base( seed, grid)
        {}

        void increment ()
        {
      	  // Access by reference the entity stored in entity pointer, and call its method next() to iterate the entity.
      	  this->dereference().next();
        }
    };


    // HierarchicIterator
    // ------------------

    template< int codim, class Grid >
    class HierarchicIterator  : public EntityPointer< codim, Grid >
    {
    	typedef EntityPointer< codim, Grid  > Base;

    protected:
    	typedef typename Base::EntityImpl EntityImpl;
    	typedef typename Base::EntitySeed EntitySeed;

    public:
      HierarchicIterator ( const EntitySeed & seed, const Grid &grid)
       : Base( seed, grid)
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
