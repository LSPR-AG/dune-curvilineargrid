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
    // [TODO]  Implement 2 Classes one for Geometry one for CachedGeometry. Introduce template parameter
    // [FIXME] Need a method to obtain RELATIVE_INTEGRAL_TOLERANCE from user. Ask Dune
    // --------------------------------------------------

    template<int mydim, int cdim, class Grid>
    class CurvGeometry
    {
      typedef CurvGeometry< mydim, cdim, Grid > This;

      typedef typename remove_const< Grid >::type::Traits Traits;
      typedef typename remove_const< Grid >::type::ctype ctype;

      //template< int, int, class > friend class Geometry;

    public:
      typedef Dune::CurvilinearGridBase<ctype,cdim>            GridBaseType;
      typedef typename GridBaseType::InternalIndexType         InternalIndexType;
      typedef typename GridBaseType::InterpolatoryOrderType    InterpolatoryOrderType;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = remove_const< Grid >::type::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef Dune::CachedCurvilinearGeometry< ctype, mydimension, coorddimension > BasicMapping;


    public:
      typedef typename BasicMapping::LocalCoordinate LocalCoordinate;
      typedef typename BasicMapping::GlobalCoordinate GlobalCoordinate;

      typedef typename BasicMapping::JacobianTransposed JacobianTransposed;
      typedef typename BasicMapping::JacobianInverseTransposed JacobianInverseTransposed;

      typedef typename BasicMapping::LocalPolynomial   Polynomial;
      typedef typename BasicMapping::PolynomialVector  PolynomialVector;


      template< class Vertices >
      CurvGeometry (
          const GeometryType &type,
          const Vertices &vertices,
          InterpolatoryOrderType order,
          const GridBaseType & gridbase)
      	      : mapping_ (type, vertices, order),
      	      gridbase_(&gridbase)
      {
          assert( int( type.dim() ) == mydimension );
      }


      CurvGeometry (
          const BasicMapping & mapping,
          const GridBaseType & gridbase)
              : mapping_ ( mapping ),
                gridbase_(&gridbase)
      { }


      CurvGeometry ( const This &other )
          : mapping_( other.mapping_ ),
            gridbase_(other.gridbase_)
      { }


      const This &operator= ( const This &other )
      {
    	  mapping_ = other.mapping_;
    	  gridbase_ = other.gridbase_;
    	  return *this;
      }

      operator bool () const { return bool( mapping_ ); }

      bool           affine ()          const  { return mapping_.affine(); }
      GeometryType   type ()            const  { return mapping_.type(); }
      int order()                       const  { return mapping_.order(); }
      int vertex (InternalIndexType i)  const  { return mapping_.vertex(i); }
      int vertices ()                   const  { return mapping_.nVertex(); }
      int corners ()                    const  { return mapping_.nCorner(); }

      GlobalCoordinate corner ( const InternalIndexType i ) const { return mapping_.corner( i ); }
      GlobalCoordinate center ()                            const { return mapping_.center(); }

      GlobalCoordinate global ( const LocalCoordinate &local ) const { return mapping_.global( local ); }
      // Local returns true if the point is inside the element. Then localC is the corresponding local coordinate
      // Local returns false if the point is not inside the element. In this case local coordinate is not defined and localC is meaningless
      LocalCoordinate local (const GlobalCoordinate &globalC) const
      {
    	  LocalCoordinate localC;
    	  bool isInside = mapping_.local(globalC, localC);

    	  if (isInside)  { return localC; }

    	  // By convention must always return a local coordinate.
    	  // So return a whatever coordinate outside the element
    	  localC[0] = -100;  localC[1] = 0;  localC[2] = 0;
    	  return localC;

    	  /*
    	  if (!isInside)  {
    		  std::cout << "searching for global coordinate " << globalC << " in the element given by " << Dune::VectorHelper::vector2string(mapping_.vertexSet()) << std::endl;
    		  DUNE_THROW( IOError, "Failed to find requested global coordinate inside the entity" );
    	  }
    	  */

      }

      // Integration Elements
      ctype integrationElement ( const LocalCoordinate &local )  const { return mapping_.integrationElement( local ); }
      Polynomial JacobianDeterminantAnalytical()                 const { return mapping_.JacobianDeterminantAnalytical(); }
      PolynomialVector NormalIntegrationElementAnalytical()      const { return mapping_.NormalIntegrationElementAnalytical(); }
      Polynomial IntegrationElementSquaredAnalytical()           const { return mapping_.IntegrationElementSquaredAnalytical(); }

      GlobalCoordinate subentityIntegrationNormal (InternalIndexType subIndex, const LocalCoordinate & localCoord)  const { return mapping_.subentityIntegrationNormal(subIndex, localCoord); }
      GlobalCoordinate subentityNormal            (InternalIndexType subIndex, const LocalCoordinate & localCoord)  const { return mapping_.subentityNormal(subIndex, localCoord); }
      GlobalCoordinate subentityUnitNormal        (InternalIndexType subIndex, const LocalCoordinate & localCoord)  const { return mapping_.subentityUnitNormal(subIndex, localCoord); }


      // Explicit integrals
      ctype integrateScalar(const Polynomial & P, double tolerance) const { return mapping_.integrateScalar(P, tolerance); }
      template <typename Functor>
      ctype integrateNumerical(const Functor & f, double tolerance) const { return mapping_.integrateNumerical(f, tolerance); }
      ctype integrateAnalyticalDot(const PolynomialVector & PVec)   const { return mapping_.integrateAnalyticalDot(PVec); }
      ctype volume () const  { return mapping_.volume(gridbase_->geometryRelativeTolerance()); }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local )                const { return mapping_.jacobianTransposed( local ); }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local )  const { return mapping_.jacobianInverseTransposed( local ); }

    private:

      const GridBaseType * gridbase_;
      BasicMapping mapping_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GEOMETRY_HH
