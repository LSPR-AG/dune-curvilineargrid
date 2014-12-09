// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_INTERSECTION_HH
#define DUNE_CURVGRID_INTERSECTION_HH

#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/curvilineargrid/curvilineargrid/entitypointer.hh>

namespace Dune
{

  namespace CurvGrid
  {

    // Intersection
    // ------------

    template< class Grid, class HostIntersection >
    class Intersection
    {
      typedef typename HostIntersection::Geometry HostGeometry;
      typedef typename HostIntersection::LocalGeometry HostLocalGeometry;

      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;

    private:
      typedef CurvGrid::IntersectionCoordVector< Grid > CoordVector;

      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;

      typedef typename Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
      typedef typename Traits::template Codim< 0 >::GeometryImpl ElementGeometryImpl;

    public:
      explicit Intersection ( const ElementGeometry &insideGeo )
        : insideGeo_( Grid::getRealImplementation( insideGeo ) ),
          hostIntersection_( 0 ),
          geo_( grid() )
      {}

      Intersection ( const Intersection &other )
        : insideGeo_( other.insideGeo_ ),
          hostIntersection_( 0 ),
          geo_( grid() )
      {}

      const Intersection &operator= ( const Intersection &other )
      {
        insideGeo_ = other.insideGeo_;
        invalidate();
        return *this;
      }

      operator bool () const { return bool( hostIntersection_ ); }

      EntityPointer inside () const
      {
        return EntityPointerImpl( insideGeo_, hostIntersection().inside() );
      }

      EntityPointer outside () const
      {
        return EntityPointerImpl( grid(), hostIntersection().outside() );
      }

      bool boundary () const { return hostIntersection().boundary(); }

      bool conforming () const { return hostIntersection().conforming(); }

      bool neighbor () const { return hostIntersection().neighbor(); }

      int boundaryId () const { return hostIntersection().boundaryId(); }

      size_t boundarySegmentIndex () const
      {
        return hostIntersection().boundarySegmentIndex();
      }

      LocalGeometry geometryInInside () const
      {
        return hostIntersection().geometryInInside();
      }

      LocalGeometry geometryInOutside () const
      {
        return hostIntersection().geometryInOutside();
      }

      Geometry geometry () const
      {
        if( !geo_ )
        {
          CoordVector coords( insideGeo_, geometryInInside() );
          geo_ = GeometryImpl( grid(), type(), coords );
        }
        return Geometry( geo_ );
      }

      GeometryType type () const { return hostIntersection().type(); }

      int indexInInside () const
      {
        return hostIntersection().indexInInside();
      }

      int indexInOutside () const
      {
        return hostIntersection().indexInOutside();
      }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        return insideGeo_.subentityIntegrationNormal( indexInInside(), geometryInInside().global( local ) );
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
    	  return insideGeo_.subentityNormal( indexInInside(), geometryInInside().global( local ) );
      }

      FieldVector< ctype, dimensionworld >
      unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
    	  return insideGeo_.subentityUnitNormal( indexInInside(), geometryInInside().global( local ) );
      }

      // TODO: This does not work in 2D
      FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
      {
        const ReferenceElement< ctype, dimension-1 > &refFace
          = ReferenceElements< ctype, dimension-1 >::general( type() );
        return unitOuterNormal( refFace.position( 0, 0 ) );
      }

      const HostIntersection &hostIntersection () const
      {
        assert( *this );
        return *hostIntersection_;
      }

      const Grid &grid () const { return insideGeo_.grid(); }

      void invalidate ()
      {
        hostIntersection_ = 0;
        geo_ = GeometryImpl( grid() );
      }

      void initialize ( const HostIntersection &hostIntersection )
      {
        assert( !(*this) );
        hostIntersection_ = &hostIntersection;
      }

    private:
      ElementGeometryImpl insideGeo_;
      const HostIntersection *hostIntersection_;
      mutable GeometryImpl geo_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
