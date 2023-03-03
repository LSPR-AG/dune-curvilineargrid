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
      : public Dune::IdSet< Grid, CurvIdSet< Grid >, typename Grid::GridFamily::CurvIdType >
    {

      typedef typename Grid::Traits Traits;
  	  typedef typename Grid::ctype  ctype;					//! coordinate type of the grid
  	  typedef typename Grid::GridFamily::CurvIdType  CurvIdType;


      typedef CurvIdSet< Grid > This;
      typedef Dune::IdSet< Grid, This, CurvIdType > Base;

    public:
	  typedef CurvIdType      IdType;


      using Base::subId;

      CurvIdSet ()  { }

      CurvIdSet ( const This &other )  { }

      const This &operator= ( const This &other )  { return *this; }

      template< int codim >
      IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
    	return entity.impl().id();
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity::codimension >( entity );
      }

      IdType subId ( const typename Traits::template Codim< 0 >::Entity &entity, int i, unsigned int codim ) const
      {
        return entity.impl().subId(i, codim);
      }

    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_IDSET_HH
