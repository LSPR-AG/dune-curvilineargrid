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

    template< class Grid >
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







      EntityPointer inside () const
      {
        return EntityPointerImpl( insideGeo_, hostIntersection().inside() );
      }

      EntityPointer outside () const
      {
        return EntityPointerImpl( grid(), hostIntersection().outside() );
      }

      // By dune-convention, domain and process boundaries are considered boundaries
      bool boundary () const {
    	  StructuralType structtype = gridbase_.entityStructuralType(1, localIndex_);

    	  return
    		(structtype == Dune::CurvilinearGridStorage::PartitionType::DomainBoundary) ||
    		(structtype == Dune::CurvilinearGridStorage::PartitionType::ProcessBoundary) ||
    		(structtype == Dune::CurvilinearGridStorage::PartitionType::ComplexBoundary);
      }


      /** \note Non-conformal grids not implemented atm  **/
      bool conforming () const { return true; }


      // By dune-convention, everything has a neighbor, except domain boundaries, and (process boundaries in serial case)
      bool neighbor () const
      {
    	  StructuralType structtype = gridbase_.entityStructuralType(1, localIndex_);

    	  if (structtype == Dune::CurvilinearGridStorage::PartitionType::DomainBoundary)  { return false; }

    	  if (
    		(structtype == Dune::CurvilinearGridStorage::PartitionType::ProcessBoundary) ||
    		(structtype == Dune::CurvilinearGridStorage::PartitionType::ComplexBoundary))
    	  {
    		  return !gridbase_.isSerial();
    	  }

    	  return true;
      }


      size_t boundarySegmentIndex () const
      {
        return hostIntersection().boundarySegmentIndex();
      }

      // Geometry of the element that calls this intersection
      LocalGeometry geometryInInside () const
      {

      }

      // Geometry of the other element
      LocalGeometry geometryInOutside () const
      {

      }

      // Geometry of this face
      Geometry geometry () const
      {
        if(!geo_)  { geo_ = GeometryImpl( grid(), type(), coords ); }
        return Geometry( geo_ );
      }

      GeometryType type () const { return geo_.type(); }

      // Face Subentity index as viewed from calling element
      InternalIndexType indexInInside () const
      {
    	  return subentityIndex(localElementIndex_);
      }

      // Face Subentity index as viewed from other element
      InternalIndexType indexInOutside () const
      {
    	  return subentityIndex(outsideElementIndex());
      }

      // All normals as viewed from the calling element

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


      const Grid &grid () const { return insideGeo_.grid(); }


      // Additional methods
      // *****************************************************************
      LocalIndexType outsideElementIndex()
      {
    	  // 1) gridbase_.face().element1index and element2index and choose the one that is not the inside
      }

      InternalIndexType subentityIndex(LocalIndexType localElementIndex)
      {
    	  // 1) Get local index of this face
    	  // 2) Loop over all subentity faces of this element
    	  // 3) Find face matching this local index, get its subentityIndex

    	  // Implement this in gridbase
      }

    private:
      LocalIndexType     localElementIndex_;
      InternalIndexType  faceSubIndex_;
      GeometryImpl       geo_;
      GridBaseType & gridbase_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
