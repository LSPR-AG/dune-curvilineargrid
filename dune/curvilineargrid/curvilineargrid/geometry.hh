// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GEOMETRY_HH
#define DUNE_CURVGRID_GEOMETRY_HH

#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/grid/common/capabilities.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // Geometry
    // [TODO] Implement 2 Classes one for Geometry one for CachedGeometry. Introduce template parameter
    // --------------------------------------------------

    template<int mydim, int cdim, class Grid>
    class CurvGeometry
    {
      typedef CurvGeometry< mydim, cdim, Grid > This;

      typedef typename remove_const< Grid >::type::Traits Traits;

      //template< int, int, class > friend class Geometry;

    public:
      typedef typename Traits::ctype ctype;
      typedef typename Traits::InternalIndexType         InternalIndexType;
      typedef typename Traits::InterpolatoryOrderType    InterpolatoryOrderType;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = Traits::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef Dune::CachedCurvilinearGeometry< ctype, mydimension, coorddimension > BasicMapping;


    public:
      typedef typename BasicMapping::LocalCoordinate LocalCoordinate;
      typedef typename BasicMapping::GlobalCoordinate GlobalCoordinate;

      typedef typename BasicMapping::JacobianTransposed JacobianTransposed;
      typedef typename BasicMapping::JacobianInverseTransposed JacobianInverseTransposed;

      typedef typename BasicMapping::Polynomial Polynomial;
      typedef typename BasicMapping::PolynomialVector PolynomialVector;

      // Types for subentity geometries

      typedef typename BasicMapping::SubentityCachedGeometry SubentityCachedGeometry;
      typedef typename BasicMapping::SubentityCachedGeometryVector SubentityCachedGeometryVector;
      typedef typename BasicMapping::SubLocalCoordinate SubLocalCoordinate;


      template< class Vertices >
      CurvGeometry (const GeometryType &type, const Vertices &vertices, InterpolatoryOrderType order )
      	  : mapping_ ( type, vertices, order )
      {
          assert( int( type.dim() ) == mydimension );
      }

      CurvGeometry (BasicMapping & mapping ) : mapping_ ( mapping )  { }

      CurvGeometry ( const This &other ) : mapping_( other.mapping_ )  { }

      const This &operator= ( const This &other )
      {
    	  mapping_ = other.mapping_;
    	  return *this;
      }

      operator bool () const { return bool( mapping_ ); }

      bool           affine ()          const  { return mapping_->affine(); }
      GeometryType   type ()            const  { return mapping_->type(); }
      int order()                       const  { return mapping_->order(); }
      int vertex (InternalIndexType i)  const  { return mapping_->vertex(i); }
      int vertices ()                   const  { return mapping_->vertices(); }
      int corners ()                    const  { return mapping_->corners(); }

      GlobalCoordinate corner ( const InternalIndexType i ) const { return mapping_->corner( i ); }
      GlobalCoordinate center ()                            const { return mapping_->center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_->global( local ); }
      // Local returns true if the point is inside the element. Then localC is the corresponding local coordinate
      // Local returns false if the point is not inside the element. In this case local coordinate is not defined and localC is meaningless
      bool local ( const GlobalCoordinate &globalC, LocalCoordinate & localC ) const { return mapping_->local(globalC, localC); }

      // Integration Elements
      ctype integrationElement ( const LocalCoordinate &local )  const { return mapping_->integrationElement( local ); }
      Polynomial JacobianDeterminantAnalytical()                 const { return mapping_->JacobianDeterminantAnalytical(); }
      PolynomialVector NormalIntegrationElementAnalytical()      const { return mapping_->NormalIntegrationElementAnalytical(); }
      Polynomial IntegrationElementSquaredAnalytical()           const { return mapping_->IntegrationElementSquaredAnalytical(); }

      GlobalCoordinate subentityIntegrationNormal (InternalIndexType subIndex, LocalCoordinate localCoord)  const { return mapping_->subentityIntegrationNormal(subIndex, localCoord); }
      GlobalCoordinate subentityNormal            (InternalIndexType subIndex, LocalCoordinate localCoord)  const { return mapping_->subentityNormal(subIndex, localCoord); }
      GlobalCoordinate subentityUnitNormal        (InternalIndexType subIndex, LocalCoordinate localCoord)  const { return mapping_->subentityUnitNormal(subIndex, localCoord); }


      // Explicit integrals
      ctype integrateScalar(const Polynomial & P, double tolerance) const { return mapping_->integrateScalar(P, tolerance); }
      template <typename Functor>
      ctype integrateNumerical(const Functor & f, double tolerance) const { return mapping_->integrateNumerical(f, tolerance); }
      ctype integrateAnalyticalDot(const PolynomialVector & PVec)   const { return mapping_->integrateAnalyticalDot(PVec); }
      ctype volume () const  { return mapping_->volume(); }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local )                const { return mapping_->jacobianTransposed( local ); }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local )  const { return mapping_->jacobianInverseTransposed( local ); }

    private:

      BasicMapping mapping_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GEOMETRY_HH
