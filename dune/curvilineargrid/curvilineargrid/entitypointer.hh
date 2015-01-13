// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_ENTITYPOINTER_HH
#define DUNE_CURVGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/entityseed.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // External Forward Declarations
    // -----------------------------

    template< int, int, class >
    class Entity;



    template< class Grid, int codim >
    class EntityPointer
    {
    	typedef typename Grid::Traits Traits;
    	typedef EntityPointer< codim, Grid > This;

    public:
        static const int dimension   = Traits::dimension;
        static const int codimension = Traits::codimension;

        typedef typename Traits::Entity Entity;

    protected:
        typedef typename Traits::EntitySeed EntitySeed;
        typedef typename Traits::template Codim<codim>::Entity  EntityImpl;


        typedef typename Traits::GridBaseType      GridBaseType;
        typedef typename Traits::IndexSetIterator  IndexSetIterator;

    public:

        EntityPointer ( IndexSetIterator & iter, GridBaseType & gridbase)
          : entity_(iter, gridbase)
        {}

        explicit EntityPointer ( const EntityImpl &entity )
          : entity_( entity )
        {}

        EntityPointer ( const This &other )
           : entity_( other.getEntity() )
        {}

        const This &operator= ( const This &other )
        {
      	  entity_ = other.getEntity();
      	  return *this;
        }

        template< class T >
        bool equals ( const EntityPointer< T> &other ) const
        {
          return entity_.equals(other.getEntity());
        }

        Entity &dereference () const  { return entity_; }

        int level () const { return 0; }

        Entity getEntity()  { return entity_; }

    protected:
        EntityImpl &entityImpl () const
        {
          return Grid::getRealImplementation( entity_ );
        }

    private:
      mutable Entity entity_;

    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_ENTITYPOINTER_HH
