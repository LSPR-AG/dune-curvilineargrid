// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GEOMETRY_HH
#define DUNE_CURVGRID_GEOMETRY_HH

//#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/curvilineargeometry/interpolation/polynomial.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/vectorhelper.hh>





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

      typedef typename Grid::Traits Traits;
      typedef typename Grid::ctype ctype;

      //template< int, int, class > friend class Geometry;

    public:
	  typedef typename Grid::GridBaseType     GridBaseType;
      typedef typename GridBaseType::InternalIndexType         InternalIndexType;
      typedef typename GridBaseType::InterpolatoryOrderType    InterpolatoryOrderType;

      static const int mydimension = mydim;
      static const int coorddimension = cdim;
      static const int dimension = Grid::dimension;
      static const int codimension = dimension - mydimension;

    protected:
      typedef typename Grid::template Codim<codimension>::EntityGeometryMappingImpl   BasicMapping;


    public:
      typedef typename BasicMapping::LocalCoordinate LocalCoordinate;
      typedef typename BasicMapping::GlobalCoordinate GlobalCoordinate;

      typedef typename BasicMapping::JacobianTransposed JacobianTransposed;
      typedef typename BasicMapping::JacobianInverseTransposed JacobianInverseTransposed;

      typedef typename BasicMapping::LocalPolynomial             Polynomial;
      typedef typename BasicMapping::PolynomialGlobalCoordinate  PolynomialGlobalCoordinate;


      template< class Vertices >
      CurvGeometry (
          const GeometryType &type,
          const Vertices &vertices,
          InterpolatoryOrderType order,
          GridBaseType & gridbase) :
          	  gridbase_(&gridbase),
          	  mapping_ (type, vertices, order)
      {
          assert( int( type.dim() ) == mydimension );
      }


      CurvGeometry (
          const BasicMapping & mapping,
          GridBaseType & gridbase)
              : gridbase_(&gridbase),
                mapping_ ( mapping )
      { }


      CurvGeometry ( const This &other )
          : gridbase_(other.gridbase_),
            mapping_( other.mapping_ )

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

    	  if (!isInside)  {
    		  std::cout << "searching for global coordinate " << globalC << " in the element given by " << VectorHelper::vector2string(mapping_.vertexSet()) << std::endl;
    		  DUNE_THROW( IOError, "Failed to find requested global coordinate inside the entity" );
    	  }

    	  return localC;
      }


      bool local(const GlobalCoordinate &globalC, LocalCoordinate &localC) const
      {
    	  return mapping_.local(globalC, localC);
      }



      // Integration Elements
      ctype integrationElement ( const LocalCoordinate &local )         const { return mapping_.integrationElement( local ); }
      Polynomial JacobianDeterminantAnalytical()                        const { return mapping_.JacobianDeterminantAnalytical(); }
      PolynomialGlobalCoordinate NormalIntegrationElementAnalytical()   const { return mapping_.NormalIntegrationElementAnalytical(); }
      Polynomial IntegrationElementSquaredAnalytical()                  const { return mapping_.IntegrationElementSquaredAnalytical(); }

      // Subentity geometry wrapper for Geometry
      // Return pointer to save computation time inside intersection class
      /*
      template<int subdim>
      CurvGeometry<subdim, cdim, Grid> * subentityGeometry(InternalIndexType subentityIndex)
      {
          return new CurvGeometry<subdim, cdim, Grid>(mapping_.template subentityGeometry<subdim>(subentityIndex), *gridbase_);
      }

      // Functions that return the normal of a face that is a subentity of this element
      GlobalCoordinate subentityIntegrationNormal (InternalIndexType subIndex, const LocalCoordinate & localCoord)  const { return mapping_.subentityIntegrationNormal(subIndex, localCoord); }
      GlobalCoordinate subentityNormal            (InternalIndexType subIndex, const LocalCoordinate & localCoord)  const { return mapping_.subentityNormal(subIndex, localCoord); }
      GlobalCoordinate subentityUnitNormal        (InternalIndexType subIndex, const LocalCoordinate & localCoord)  const { return mapping_.subentityUnitNormal(subIndex, localCoord); }

      // Explicit integrals
      ctype integrateScalar(const Polynomial & P, double tolerance) const { return mapping_.integrateScalar(P, tolerance); }
      template <typename Functor>
      ctype integrateNumerical(const Functor & f, double tolerance) const { return mapping_.integrateNumerical(f, tolerance); }
      ctype integrateAnalyticalDot(const PolynomialVector & PVec)   const { return mapping_.integrateAnalyticalDot(PVec); }
      */

      ctype volume () const  { return mapping_.volume(gridbase_->property().geometryRelativeTolerance()); }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local )                const { return mapping_.jacobianTransposed( local ); }
      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local )  const { return mapping_.jacobianInverseTransposed( local ); }


      const BasicMapping & basegeometry()  { return mapping_; }


    private:

      GridBaseType * gridbase_;
      BasicMapping mapping_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GEOMETRY_HH
