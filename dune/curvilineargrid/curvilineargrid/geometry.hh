// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_GEOMETRY_HH
#define DUNE_GEOGRID_GEOMETRY_HH

#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/polynomial/polynomial.hh>
#include <dune/geometry/polynomialinterpolation/curvilinearelementinterpolator.hh>
#include <dune/geometry/lagrangegeometry.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/geometrygrid/cornerstorage.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // InferHasSingleGeometryType
    // --------------------------

    template< class hasSingleGeometryType, int dim, int mydim >
    struct InferHasSingleGeometryType
    {
    private:
      static const unsigned int id = hasSingleGeometryType::topologyId;
      static const unsigned int idMask = (1u << mydim) - 1u;

    public:
      static const bool v = hasSingleGeometryType::v && ((mydim == dim) || ((id | 1u) == 1u) || ((id | 1u) == idMask));
      static const unsigned int topologyId = (v ? id & idMask : ~0u);
    };

    template< class hasSingleGeometryType, int dim >
    struct InferHasSingleGeometryType< hasSingleGeometryType, dim, 1 >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry::CubeTopology< 1 >::type::id;
    };

    template< class hasSingleGeometryType, int dim >
    struct InferHasSingleGeometryType< hasSingleGeometryType, dim, 0 >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry::CubeTopology< 0 >::type::id;
    };



    // GeometryTraits
    // --------------

    template< class Grid >
    struct GeometryTraits
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< ctype > > MatrixHelper;

      static ctype tolerance () { return 16 * std::numeric_limits< ctype >::epsilon(); }

      template< int mydim, int cdim >
      struct CornerStorage
      {
        typedef GeoGrid::CornerStorage< mydim, cdim, Grid > Type;
      };

      template< int mydim >
      struct hasSingleGeometryType
        : public InferHasSingleGeometryType< Capabilities::hasSingleGeometryType< Grid >, Traits::dimension, mydim >
      {};
    };



    // Geometry
    // --------

    template< int mydim, int cdim, class Grid >
    class Geometry
    {
      typedef Geometry< mydim, cdim, Grid > This;

      typedef typename remove_const< Grid >::type::Traits Traits;

      template< int, int, class > friend class Geometry;

    public:
      typedef typename Traits::ctype ctype;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = Traits::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef CachedLagrangeGeometry< ctype, mydimension, coorddimension, GeometryTraits< Grid > > BasicMapping;

      struct Mapping
        : public BasicMapping
      {
        template< class Vertices >
        Mapping ( const GeometryType &type, const Vertices &vertices, int order )
          : BasicMapping( type, vertices, order ),
            refCount_( 0 )
        {}


        void addReference () { ++refCount_; }
        bool removeReference () { return (--refCount_ == 0); }

      private:
        unsigned int refCount_;
      };

    public:
      typedef typename Mapping::LocalCoordinate LocalCoordinate;
      typedef typename Mapping::GlobalCoordinate GlobalCoordinate;

      typedef typename Mapping::JacobianTransposed JacobianTransposed;
      typedef typename Mapping::JacobianInverseTransposed JacobianInverseTransposed;

      typedef typename Mapping::Polynomial Polynomial;
      typedef typename Mapping::PolynomialVector PolynomialVector;

      // Types for subentity geometries

      typedef typename Mapping::SubentityCachedGeometry SubentityCachedGeometry;
      typedef typename Mapping::SubentityCachedGeometryVector SubentityCachedGeometryVector;
      typedef typename Mapping::SubLocalCoordinate SubLocalCoordinate;


      Geometry ( const Grid &grid )
        : grid_( &grid ),
          mapping_( nullptr )
      {}

      template< class Vertices >
      Geometry ( const Grid &grid, const GeometryType &type, const Vertices &vertices, int order )
        : grid_( &grid )
      {
        assert( int( type.dim() ) == mydimension );
        void *mappingStorage = grid.allocateStorage( sizeof( Mapping ) );
        mapping_ = new( mappingStorage ) Mapping( type, vertices, order );
        mapping_->addReference();
      }

      Geometry ( const This &other )
        : grid_( other.grid_ ),
          mapping_( other.mapping_ )
      {
        if( mapping_ )
          mapping_->addReference();
      }

      ~Geometry ()
      {
        if( mapping_ && mapping_->removeReference() )
          destroyMapping();
      }

      const This &operator= ( const This &other )
      {
        if( other.mapping_ )
          other.mapping_->addReference();
        if( mapping_ && mapping_->removeReference() )
          destroyMapping();
        grid_ = other.grid_;
        mapping_ = other.mapping_;
        return *this;
      }

      operator bool () const { return bool( mapping_ ); }

      bool affine () const { return mapping_->affine(); }
      GeometryType type () const { return mapping_->type(); }

      int order() const { return mapping_->order(); }

      int vertex (int i) const { return mapping_->vertex(i); }

      int vertices () const { return mapping_->vertices(); }

      int corners () const { return mapping_->corners(); }
      GlobalCoordinate corner ( const int i ) const { return mapping_->corner( i ); }
      GlobalCoordinate center () const { return mapping_->center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_->global( local ); }
      // Local returns true if the point is inside the element. Then localC is the corresponding local coordinate
      // Local returns false if the point is not inside the element. In this case local coordinate is not defined and localC is meaningless
      bool local ( const GlobalCoordinate &globalC, LocalCoordinate & localC ) const { return mapping_->local(globalC, localC); }

      // Cached Geometry classes for all dim-1 subentities
      SubentityCachedGeometryVector subentityCachedGeometries() const  { return mapping_->subentityCachedGeometries(); }

      // Normals
      GlobalCoordinate normal(const LocalCoordinate &local ) const { return mapping_->normal(local); }
      GlobalCoordinate subentityNormal(int subentityNo, const SubLocalCoordinate &local ) const { return mapping_->subentityNormal(subentityNo, local); }


      // Integration Elements
      ctype integrationElement ( const LocalCoordinate &local ) const { return mapping_->integrationElement( local ); }
      Polynomial JacobianDeterminantAnalytical() const { return mapping_->JacobianDeterminantAnalytical(); }
      PolynomialVector NormalIntegrationElementAnalytical() const { return mapping_->NormalIntegrationElementAnalytical(); }
      Polynomial IntegrationElementSquaredAnalytical() const { return mapping_->IntegrationElementSquaredAnalytical(); }

      // Explicit integrals
      ctype integrateScalar(const Polynomial & P, double tolerance) const { return mapping_->integrateScalar(P, tolerance); }
      template <typename Functor>
      ctype integrateNumerical(const Functor & f, double tolerance) const { return mapping_->integrateNumerical(f, tolerance); }
      ctype integrateAnalyticalDot(const PolynomialVector & PVec) const { return mapping_->integrateAnalyticalDot(PVec); }
      ctype volume () const { return mapping_->volume(); }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianTransposed( local ); }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const { return mapping_->jacobianInverseTransposed( local ); }

      const Grid &grid () const { return *grid_; }

    private:
      void destroyMapping ()
      {
        mapping_->~Mapping();
        grid().deallocateStorage( mapping_, sizeof( Mapping ) );
      }

      const Grid *grid_;
      Mapping* mapping_;
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_GEOMETRY_HH
