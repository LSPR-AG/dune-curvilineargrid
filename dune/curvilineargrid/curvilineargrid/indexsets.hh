// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_INDEXSETS_HH
#define DUNE_CURVGRID_INDEXSETS_HH

#include <vector>

//#include <dune/common/typetraits.hh>

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
        //typename Grid::LocalIndexType
        >
    {
    	typedef typename Grid::Traits Traits;
    	typedef typename Grid::ctype ctype;

    public:
      static const int dimension = Grid::dimension;

    private:

	    typedef typename Grid::GridBaseType     GridBaseType;
    	typedef typename GridBaseType::LocalIndexType         LocalIndexType;


    	typedef CurvIndexSet< Grid > This;
    	typedef Dune::IndexSet< Grid, This, LocalIndexType > Base;


    public:


      typedef LocalIndexType   IndexType;

      CurvIndexSet (GridBaseType & gridbase)
          : gridbase_(gridbase),
      	    geomtypes_
      	    {
  		      std::vector<GeometryType> (1, GeometryType ( Dune::GenericGeometry::SimplexTopology<3>::type::id, 3)),
  		      std::vector<GeometryType> (1, GeometryType ( Dune::GenericGeometry::SimplexTopology<2>::type::id, 2)),
  		      std::vector<GeometryType> (1, GeometryType ( Dune::GenericGeometry::SimplexTopology<1>::type::id, 1)),
  		      std::vector<GeometryType> (1, GeometryType ( Dune::GenericGeometry::SimplexTopology<0>::type::id, 0))
      	    }
      {

      }

      CurvIndexSet ( const This &other )
        : gridbase_(other.gridbase_),
          geomtypes_(other.geomtypes_)
      { }

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


      template< int codim >
      IndexType indexBase ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
        return Grid::getRealImplementation( entity ).indexBase();
      }

      template< int codim >
      IndexType subIndexBase ( const typename Traits::template Codim< codim >::Entity &entity, int i, unsigned int subcodim ) const
      {
        return Grid::getRealImplementation( entity ).subIndexBase( i, subcodim );
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

      std::vector< GeometryType > types ( int codim ) const  { return geomtypes_[codim]; }

      const std::vector< GeometryType > &geomTypes ( int codim ) const  { return geomtypes_[codim]; }

    private:

      GridBaseType & gridbase_;
      std::vector<std::vector< GeometryType > > geomtypes_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INDEXSETS_HH
