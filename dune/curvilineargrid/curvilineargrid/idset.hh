// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_IDSET_HH
#define DUNE_CURVGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // CurvIdSet
    // -----

    template< class Grid >
    class CurvIdSet
      : public Dune::IdSet< Grid, CurvIdSet< Grid >, typename Dune::CurvilinearGridBase::IdType >
    {

      typedef CurvIdSet< Grid > This;
      typedef Dune::IdSet< Grid, This, typename Dune::CurvilinearGridBase::IdType > Base;

      typedef typename remove_const< Grid >::type::Traits Traits;
	  typedef typename Traits::ctype ctype;						//! coordinate type of the grid

	  static const int dimensionworld = Traits::dimensionworld;		//! dimension of the world

    public:
      typedef Dune::CurvilinearGridBase<ctype, dimensionworld>::IdType  IdType;


      using Base::subId;

      CurvIdSet ()  { }

      CurvIdSet ( const This &other )  { }

      const This &operator= ( const This &other )  { return *this; }

      template< int codim >
      IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).id();
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity::codimension >( entity );
      }

      IdType subId ( const typename Traits::template Codim< 0 >::Entity &entity, int i, unsigned int codim ) const
      {
        return Grid::getRealImplementation( entity ).subId(i, codim);
      }

    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_IDSET_HH
