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
    class CurvIndexSet
      : public Dune::IndexSet< Grid, CurvIndexSet< Grid>,
        unsigned int
        //typename remove_const< Grid >::type::LocalIndexType
        >
    {
    	typedef typename remove_const< Grid >::type::Traits Traits;
    	typedef typename remove_const< Grid >::type::ctype ctype;

    public:
      static const int dimension = remove_const< Grid >::type::dimension;

    private:

        typedef Dune::CurvilinearGridBase<ctype,dimension>    GridBaseType;
    	typedef typename GridBaseType::LocalIndexType         LocalIndexType;


    	typedef CurvIndexSet< Grid > This;
    	typedef Dune::IndexSet< Grid, This, LocalIndexType > Base;


    public:


      typedef LocalIndexType   IndexType;

      CurvIndexSet (GridBaseType & gridbase) : gridbase_(gridbase)
      { }

      CurvIndexSet ( const This &other ) : gridbase_(other.gridbase_) { }

      const This &operator= ( const This &other )  { return *this; }

      template< int codim >
      IndexType index ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).index();
      }

      template< int codim >
      IndexType subIndex ( const typename Traits::template Codim< codim >::Entity &entity, int i, unsigned int subcodim ) const
      {
        return Grid::getRealImplementation( entity ).subIndex( i, subcodim );
      }


      IndexType size ( GeometryType type ) const
      {
    	  if (!type.isSimplex()) { return 0; }
    	  else
    	  {
    		  return size(dimension - type.dim());
    	  }
      }


      IndexType size ( int codim ) const  { return gridbase_.nEntity(codim); }

      template< class EntityType >
      bool contains ( const EntityType &entity ) const
      {
    	  int localIndex = Grid::getRealImplementation(entity).index();
    	  int globalIndex;

    	  return gridbase_.findEntityGlobalIndex(EntityType::codimension, localIndex, globalIndex);
      }

      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
    	  std::vector<GeometryType>  rez;

    	  switch (codim)
    	  {
    	  case 0 :  rez.push_back( Dune::GeometryType ( Dune::GenericGeometry::SimplexTopology<3>::type::id, 3) );  break;
    	  case 1 :  rez.push_back( Dune::GeometryType ( Dune::GenericGeometry::SimplexTopology<2>::type::id, 2) );  break;
    	  case 2 :  rez.push_back( Dune::GeometryType ( Dune::GenericGeometry::SimplexTopology<1>::type::id, 1) );  break;
    	  case 3 :  rez.push_back( Dune::GeometryType ( Dune::GenericGeometry::SimplexTopology<0>::type::id, 0) );  break;
    	  }

    	  return rez;
      }

    private:

      GridBaseType & gridbase_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INDEXSETS_HH
