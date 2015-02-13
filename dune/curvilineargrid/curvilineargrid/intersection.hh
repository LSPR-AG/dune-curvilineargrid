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

    // [FIXME]  Create constructor that allows passing inside-geometry, since it need not be computed several times when iterating over intersections
    // [FIXME]  Since geo_ is not a pointer, must fix (if !geo) statements


    // Intersection
    // ------------

    template< class Grid >
    class CurvIntersection
    {

      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      typedef typename remove_const< Grid >::type::ctype ctype;

      static const int dimension = remove_const< Grid >::type::dimension;
      static const int dimensionworld = remove_const< Grid >::type::dimensionworld;

      typedef Dune::CurvilinearGridStorage<ctype,dimension>    GridStorageType;
      typedef Dune::CurvilinearGridBase<ctype,dimension>       GridBaseType;

      // Codimensions of entity types for better code readability
      static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
      static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
      static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
      static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

      typedef typename GridBaseType::LocalIndexType             LocalIndexType;
      typedef typename GridBaseType::InternalIndexType          InternalIndexType;
      typedef typename GridBaseType::StructuralType             StructuralType;
      typedef typename GridBaseType::InterpolatoryOrderType     InterpolatoryOrderType;

      typedef typename GridBaseType::IndexSetIterator   IndexSetIterator;



      typedef typename Traits::template Codim< ELEMENT_CODIM >::Entity Entity;
      typedef typename Traits::template Codim< ELEMENT_CODIM >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< FACE_CODIM >::Geometry Geometry;
      typedef typename Traits::template Codim< FACE_CODIM >::LocalGeometry LocalGeometry;

      typedef typename Traits::template Codim< ELEMENT_CODIM >::Geometry ElementGeometry;

    private:

      typedef typename Dune::CurvGrid::CurvEntityPointer<ELEMENT_CODIM, Grid>  EntityPointerImpl;

      typedef typename Traits::template Codim< FACE_CODIM >::GeometryImpl GeometryImpl;
      typedef typename Traits::template Codim< ELEMENT_CODIM >::GeometryImpl ElementGeometryImpl;

      typedef Dune::FieldVector< ctype, dimension-1 >      LocalCoordinate;
      typedef Dune::FieldVector< ctype, dimensionworld >   GlobalCoordinate;


      typedef typename remove_const< Grid >::type::template BaseCodim<ELEMENT_CODIM>::EntityGeometryImpl  ElementBaseGeometry;
      typedef typename remove_const< Grid >::type::template BaseCodim<FACE_CODIM>::EntityGeometryImpl     FaceBaseGeometry;

      const InterpolatoryOrderType LINEAR_ELEMENT_ORDER = 1;

    public:
      CurvIntersection (
    		  LocalIndexType localIndexInside,    // Index of the element wrt which this intersection is calculated
    		  InternalIndexType subIndexInside,   // Internal index of the face as the element subentity
    		  GridBaseType & gridbase)
    	: localIndexInside_(localIndexInside),
    	  subIndexInside_(subIndexInside),
    	  gridbase_(&gridbase)
      {
    	  localFaceIndex_ = gridbase.subentityLocalIndex (localIndexInside, 0, 1, subIndexInside);
    	  insideGeo_ = LocalGeometry(LocalNeighborGeometry(localIndexInside_, subIndexInside_));

    	  computeOutside();
      }

      CurvIntersection ( const CurvIntersection &other )
  	    : localIndexInside_(other.localIndexInside_),
  	      subIndexInside_(other.subIndexInside_),
  	      gridbase_(other.gridbase_),
  	      insideGeo_(LocalNeighborGeometry(localIndexInside_, subIndexInside_)),
  	      localFaceIndex_(other.localFaceIndex_),
  	      localIndexOutside_(other.localIndexOutside_),
  	      subIndexOutside_(other.subIndexOutside_)
      {

      }

      const CurvIntersection &operator= ( const CurvIntersection &other )
      {
    	  localIndexInside_ = other.localIndexInside_;
    	  subIndexInside_ = other.subIndexInside_;
    	  gridbase_ = other.gridbase_;
  	      localFaceIndex_ = other.localFaceIndex_;
  	      localIndexOutside_ = other.localIndexOutside_;
  	      subIndexOutside_ = other.subIndexOutside_;

    	  insideGeo_ = LocalNeighborGeometry(localIndexInside_, subIndexInside_);
    	  return *this;
      }

      bool equals(const CurvIntersection& other) const
      {
          return (localIndexInside_ == other.localIndexInside_) && (subIndexInside_ == other.subIndexInside_);
      }


      EntityPointer inside () const
      {
    	  IndexSetIterator thisIter = gridbase_->entityIndexIterator(ELEMENT_CODIM, localIndexInside_);
    	  StructuralType structtype = gridbase_->entityStructuralType(ELEMENT_CODIM, localIndexInside_);
    	  return EntityPointer(EntityPointerImpl(thisIter, *gridbase_, All_Partition));
      }

      EntityPointer outside () const
      {
    	  IndexSetIterator thisIter = gridbase_->entityIndexIterator(ELEMENT_CODIM, localIndexOutside_);
    	  StructuralType structtype = gridbase_->entityStructuralType(ELEMENT_CODIM, localIndexOutside_);
    	  return EntityPointer(EntityPointerImpl(thisIter, *gridbase_, All_Partition));
      }

      // By dune-convention, domain and periodic boundaries are considered boundaries
      bool boundary () const {
    	  StructuralType structtype = gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_);
    	  return (structtype == GridStorageType::PartitionType::DomainBoundary);
      }


      /** \note Non-conformal grids not implemented atm  **/
      bool conforming () const { return true; }


      // By dune-convention, everything has a neighbour, except domain boundaries, and (process boundaries in serial case)
      bool neighbor () const
      {
    	  StructuralType structtype = gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_);

    	  if (structtype == GridStorageType::PartitionType::DomainBoundary)  { return false; }
    	  if (structtype == GridStorageType::PartitionType::ProcessBoundary) { return !gridbase_->isSerial(); }

    	  return true;
      }


      size_t boundarySegmentIndex () const
      {
    	  assert(boundary());  // If this entity is not a Domain Boundary should throw error
    	  return gridbase_->boundarySegmentIndex(localFaceIndex_);
      }

      // Geometry of the intersection from within the inside element
      LocalGeometry geometryInInside () const  {
    	  return LocalGeometry(insideGeo_);
      }

      // Geometry of the intersection from within the outside element
      LocalGeometry geometryInOutside () const  {
    	  return LocalGeometry(LocalNeighborGeometry(localIndexOutside_, subIndexOutside_));
      }

      // Geometry of this face
      Geometry geometry () const
      {

        if(!geo_)  {
        	InterpolatoryOrderType interporder = gridbase_->entityInterpolationOrder(ELEMENT_CODIM, localIndexInside_);

        	ElementBaseGeometry entityGeometryBase = gridbase_->entityGeometry(ELEMENT_CODIM, localIndexInside_);
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
      // If intersection is a ghost intersection, skip it, because the intersection is not allowed to be of type ghost
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
    		  localFaceIndex_ = gridbase_->subentityLocalIndex (localIndexInside_, ELEMENT_CODIM, FACE_CODIM, subIndexInside_);
    		  inc = (gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_) == GridStorageType::PartitionType::Ghost);
    	  }

    	  // If this is not the end, update contents of the intersection
    	  if (subIndexInside_ != SUBENTITY_SIZE)  { computeOutside(); }
      }

      // Finds outside entity local index, and this face subentity index
      void computeOutside()
      {
    	  LocalIndexType neighborElementIndex0 = gridbase_->faceNeighbor(localFaceIndex_, 0);
    	  LocalIndexType neighborElementIndex1 = gridbase_->faceNeighbor(localFaceIndex_, 1);

    	  localIndexOutside_ = (neighborElementIndex0 == localIndexInside_) ? neighborElementIndex1 : neighborElementIndex0;

    	  int iFace = 0;
    	  bool found_face = false;
    	  while (!found_face)
    	  {
    		  if (iFace >= 4)  { DUNE_THROW(Dune::IOError, "Intersection: Not found face as its outside-element subentity"); }

    		  if (localFaceIndex_ == gridbase_->subentityLocalIndex(localIndexOutside_, ELEMENT_CODIM, FACE_CODIM, iFace))
    		  {
    			  found_face = true;
    			  subIndexOutside_ = iFace;
    		  }
    		  iFace++;
    	  }
      }


      ElementGeometryImpl LocalNeighborGeometry(LocalIndexType neighborIndex, InternalIndexType subentityIndex)
      {
    	  Dune::GeometryType gt = gridbase_->entityGeometryType(ELEMENT_CODIM, neighborIndex);
    	  return GeometryImpl(type(), refCoord(gt, subentityIndex), LINEAR_ELEMENT_ORDER);
      }

      // Creates coordinates of a linear intersection as one of the faces of the reference element
      std::vector<GlobalCoordinate> refCoord(Dune::GeometryType elemGT, InternalIndexType subentityIndex)
	  {
    	  std::vector<GlobalCoordinate> rez;

    	  const Dune::ReferenceElement<ctype, dimension> & ref = Dune::ReferenceElements<ctype, dimension>::general(elemGT);
    	  int nCornerPerFace = ref.size(0, FACE_CODIM, VERTEX_CODIM);

    	  for (int i = 0; i < nCornerPerFace; i++)
    	  {
    		  InternalIndexType cornerInd = ref.subEntity(subentityIndex, FACE_CODIM, i, VERTEX_CODIM);
    		  GlobalCoordinate  coord = Dune::CurvilinearGeometryHelper::cornerInternalCoordinate(elemGT, cornerInd);
    		  rez.push_back(coord);
    	  }

    	  return rez;
	  }




    private:
      LocalIndexType      localFaceIndex_;
      LocalIndexType      localIndexInside_;
      LocalIndexType      localIndexOutside_;
      InternalIndexType   subIndexInside_;
      InternalIndexType   subIndexOutside_;
      GeometryImpl        geo_;
      ElementGeometryImpl insideGeo_;
      GridBaseType *      gridbase_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
