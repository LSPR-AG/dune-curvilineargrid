// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_INDEXSETS_HH
#define DUNE_CURVGRID_INDEXSETS_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/curvilineargrid/curvilineargrid/declaration.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // IndexSet
    // --------

    template<class Grid>
    class IndexSet
      : public Dune::IndexSet< Grid, IndexSet< Grid> >
    {
      typedef IndexSet< Grid > This;
      typedef Dune::IndexSet< Grid, This > Base;

      typedef typename Dune::CurvilinearGridBase      CurvilinearGridBase;

      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      static const int dimension = Traits::dimension;

      typedef typename Base::IndexType IndexType;

      IndexSet ()  { }

      IndexSet ( const This &other )  { }

      const This &operator= ( const This &other )  { return *this; }

      using Base::index;
      using Base::subIndex;

      template< int cc >
      IndexType index ( const typename Traits::template Codim< cc >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).index( hostIndexSet() );
      }

      template< int cc >
      IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &entity, int i, unsigned int codim ) const
      {
        return Grid::getRealImplementation( entity ).subIndex( hostIndexSet(), i, codim );
      }

      IndexType size ( GeometryType type ) const
      {
    	  if (!type.isSimplex()) { return 0; }
    	  return size(dimension - type.dim());
      }

      IndexType size ( int codim ) const
      {
    	  switch(codim)
    	  {
    	  case 0  :   return gridbase_.nElement();  break;
    	  case 1  :   return gridbase_.nFace();     break;
    	  case 2  :   return gridbase_.nEdge();     break;
    	  case 3  :   return gridbase_.nVertex();   break;
    	  }
    	  return 0;
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).isContained( hostIndexSet() );
      }

      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

      operator bool () const { return bool( hostIndexSet_ ); }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( *this );
        return *hostIndexSet_;
      }

      CurvilinearGridBase & gridbase_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INDEXSETS_HH
