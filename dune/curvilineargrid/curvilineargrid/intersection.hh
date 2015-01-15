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

      typedef typename Grid::Traits Traits;

    public:
      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;


      typedef typename Traits::GridStorageType    GridStorageType;
      typedef typename Traits::GridBaseType       GridBaseType;

      typedef typename Traits::LocalIndexType             LocalIndexType;
      typedef typename Traits::InternalIndexType          InternalIndexType;
      typedef typename Traits::StructuralType             StructuralType;
      typedef typename Traits::InterpolatoryOrderType     InterpolatoryOrderType;

      typedef typename Traits::IndexSetIterator   IndexSetIterator;

    private:

      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;

      typedef typename Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
      typedef typename Traits::template Codim< 0 >::GeometryImpl ElementGeometryImpl;

      typedef Dune::FieldVector< ctype, dimension-1 >      LocalCoordinate;
      typedef Dune::FieldVector< ctype, dimensionworld >   GlobalCoordinate;


      typedef typename Traits::template BaseCodim<0>::EntityGeometry  ElementBaseGeometry;
      typedef typename Traits::template BaseCodim<1>::EntityGeometry  FaceBaseGeometry;

    public:
      Intersection (
    		  LocalIndexType localIndexInside,
    		  InternalIndexType subIndexInside,
    		  GridBaseType & gridbase)
    	: localIndexInside_(localIndexInside),
    	  subIndexInside_(subIndexInside),
    	  gridbase_(gridbase)
      {
    	  localFaceIndex_ = gridbase.subentityLocalIndex (localIndexInside, 0, 1, subIndexInside);

    	  computeOutside();
      }

      Intersection ( const Intersection &other )
        : insideGeo_( other.insideGeo_ )
      {}

      const Intersection &operator= ( const Intersection &other )
      {
        insideGeo_ = other.insideGeo_;
        return *this;
      }







      EntityPointer inside () const
      {
    	  IndexSetIterator thisIter = gridbase_.entityIndexIterator(0, localIndexInside_);
    	  return EntityPointerImpl(thisIter, gridbase_);
      }

      EntityPointer outside () const
      {
    	  IndexSetIterator thisIter = gridbase_.entityIndexIterator(0, localIndexOutside_);
    	  return EntityPointerImpl(thisIter, gridbase_);
      }

      // By dune-convention, domain and process boundaries are considered boundaries
      bool boundary () const {
    	  StructuralType structtype = gridbase_.entityStructuralType(1, localFaceIndex_);

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
    	  StructuralType structtype = gridbase_.entityStructuralType(1, localFaceIndex_);

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
      LocalGeometry geometryInInside () const  { return Geometry(type(), refCoord(subIndexInside_), 1); }

      // Geometry of the other element
      LocalGeometry geometryInOutside () const  { return Geometry(type(), refCoord(subIndexOutside_), 1); }

      // Geometry of this face
      Geometry geometry () const
      {

        if(!geo_)  {
        	InterpolatoryOrderType interporder = gridbase_.entityInterpolationOrder(0, localIndexInside_);

        	ElementBaseGeometry entityGeometryBase = gridbase_.entityGeometry(0, localIndexInside_);
        	FaceBaseGeometry faceGeometryBase = entityGeometryBase.subentityGeometry(subIndexInside_);

        	geo_ = GeometryImpl(faceGeometryBase);
        }

        return Geometry( geo_ );
      }

      GeometryType type () const { return geo_.type(); }

      // Face Subentity index as viewed from calling element
      InternalIndexType indexInInside () const  { return subIndexInside_; }

      // Face Subentity index as viewed from other element
      InternalIndexType indexInOutside () const  { return subIndexOutside_; }

      // All normals as viewed from the calling element

      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate &localCoord ) const
      {
        return insideGeo_.subentityIntegrationNormal( indexInInside(), localCoord );
      }

      GlobalCoordinate outerNormal ( const LocalCoordinate &localCoord ) const
      {
    	  return insideGeo_.subentityNormal( indexInInside(), localCoord );
      }

      GlobalCoordinate unitOuterNormal ( const LocalCoordinate &localCoord ) const
      {
    	  return insideGeo_.subentityUnitNormal( indexInInside(), localCoord );
      }

      // TODO: This does not work in 2D
      GlobalCoordinate centerUnitOuterNormal () const
      {
    	  const ReferenceElement< ctype, dimension-1 > &refFace = ReferenceElements< ctype, dimension-1 >::general( type() );
    	  return unitOuterNormal( refFace.position( 0, 0 ) );
      }


      // Auxiliary methods
      // *******************************************************

      LocalIndexType intersectionIndex() { return localFaceIndex_; }

      // Iterates over subentities of the inside entity by increasing the inside element subentity index
      // If intersection is a ghost intersection, skip it
      // If next() is called beyond the allowed size, throw error
      void next()
      {
    	  const int SUBENTITY_SIZE = 4;

    	  if (subIndexInside_ >= SUBENTITY_SIZE) { DUNE_THROW(Dune::IOError, "intersection: next() called with unexpected subentity index"); }

    	  // Increase iterator until find a non-ghost face or reach the end
    	  bool inc = true;
    	  while( inc && (subIndexInside_ < SUBENTITY_SIZE) )
    	  {
    		  subIndexInside_++;
    		  localFaceIndex_ = gridbase_.subentityLocalIndex (localIndexInside_, 0, 1, subIndexInside_);
    		  inc = (gridbase_.entityStructuralType(1, localFaceIndex_) == GridStorageType::PartitionType::Ghost);
    	  }

    	  // If this is not the end, update contents of the intersection
    	  if (subIndexInside_ != SUBENTITY_SIZE)  { computeOutside(); }
      }

      // Finds outside entity local index, and this face subentity index
      void computeOutside()
      {
    	  LocalIndexType tmpIndex1 = gridbase_.faceNeighbor(localFaceIndex_, 0);
    	  LocalIndexType tmpIndex2 = gridbase_.faceNeighbor(localFaceIndex_, 1);

    	  localIndexOutside_ = (tmpIndex1 == localIndexInside_) ? tmpIndex2 : tmpIndex1;

    	  for (InternalIndexType iFace = 0; iFace < 4; iFace++)
    	  {
    		  if (localFaceIndex_ == gridbase_.subentityLocalIndex(localIndexOutside_, 0, 1, iFace))  { subIndexOutside_ = iFace; }
    	  }
      }

      // Creates coordinates of a linear intersection as one of the faces of the reference element
      std::vector<GlobalCoordinate> refCoord(InternalIndexType subentityIndex)
	  {
    	  std::vector<InternalIndexType> referenceSubset = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(type(), 1, subentityIndex);

    	  std::vector<ctype> referenceCoord {
    		  {0.0, 0.0, 0.0},
    		  {1.0, 0.0, 0.0},
    		  {0.0, 1.0, 0.0},
    		  {0.0, 0.0, 1.0}
    	  };

    	  std::vector<GlobalCoordinate> coord;

    	  for (int i = 0; i < referenceSubset.size(); i++)  { coord.push_back(referenceCoord[referenceSubset[i]]); }

    	  return coord;
	  }




    private:
      LocalIndexType     localFaceIndex_;
      LocalIndexType     localIndexInside_;
      LocalIndexType     localIndexOutside_;
      InternalIndexType  subIndexInside_;
      InternalIndexType  subIndexOutside_;
      GeometryImpl       geo_;
      LocalGeometry      insideGeo_;
      GridBaseType &     gridbase_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
