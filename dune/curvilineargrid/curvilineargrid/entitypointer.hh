// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_ENTITYPOINTER_HH
#define DUNE_CURVGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/entity.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // External Forward Declarations
    // -----------------------------

    template< int, int, class >
    class Entity;



    template< int codim, class Grid >
    class CurvEntityPointer
    {
    	typedef typename remove_const< Grid >::type::Traits Traits;
    	typedef typename remove_const< Grid >::type::ctype  ctype;

    	typedef CurvEntityPointer< codim, Grid > This;

    public:
        static const int dimension     = remove_const< Grid >::type::dimension;
        static const int codimension   = codim;

        typedef typename Traits::template Codim<codim>::Entity      Entity;
        typedef typename Traits::template Codim<codim>::EntitySeed  EntitySeed;

    protected:

        typedef Dune::CurvGrid::CurvEntity<codim, dimension, Grid>  EntityImpl;

        typedef Dune::CurvilinearGridBase<ctype,dimension>    GridBaseType;
        typedef typename GridBaseType::IndexSetIterator       IndexSetIterator;

    public:

        CurvEntityPointer ( const IndexSetIterator & iter, const GridBaseType & gridbase, PartitionIteratorType pitype)
          : entity_(EntityImpl(iter, gridbase, pitype))
        {}

        CurvEntityPointer ( EntitySeed seed)
        {
        	IndexSetIterator iter = seed.gridBase().entityIndexDuneIterator(codim, seed.partitionType(), seed.localIndex());
        	entity_ = Entity(EntityImpl(iter, seed.gridBase(), seed.partitionType()));
        }

        CurvEntityPointer ( const EntityImpl &entity )
          : entity_( entity )
        {}

        CurvEntityPointer ( const EntityImpl && entity )
          : entity_( std::move(entity) )
        {}

        //const This &operator= ( const This &other )  { entity_ = other.getEntity(); return *this; }

        bool equals ( const CurvEntityPointer &other ) const
        {
          return entity_.equals(other.dereference());
        }

        const Entity &dereference () const  { return entity_; }

        int level () const { return 0; }

    protected:
        // EntityImpl &entityImpl () const  { return Grid::getRealImplementation( entity_ ); }

    private:
      Entity entity_;

    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_ENTITYPOINTER_HH
