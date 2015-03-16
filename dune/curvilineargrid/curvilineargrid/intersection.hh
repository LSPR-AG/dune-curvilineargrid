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

	  typedef typename remove_const< Grid >::type::GridStorageType  GridStorageType;
	  typedef typename remove_const< Grid >::type::GridBaseType     GridBaseType;

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


      typedef typename remove_const< Grid >::type::template Codim<ELEMENT_CODIM>::EntityGeometryMappingImpl  ElementBaseGeometry;
      typedef typename remove_const< Grid >::type::template Codim<FACE_CODIM>::EntityGeometryMappingImpl     FaceBaseGeometry;

      const InterpolatoryOrderType LINEAR_ELEMENT_ORDER = 1;

    public:
      CurvIntersection (
    		  LocalIndexType localIndexInside,    // Index of the element wrt which this intersection is calculated
    		  InternalIndexType subIndexInside,   // Internal index of the face as the element subentity
    		  GridBaseType & gridbase)
    	: localIndexInside_(localIndexInside),
    	  subIndexInside_(subIndexInside),
    	  gridbase_(&gridbase),
    	  insideGeo_(gridbase.template entityGeometry<ELEMENT_CODIM>(localIndexInside), gridbase),
    	  geo_(nullptr),
    	  geoInInside_(nullptr),
    	  geoInOutside_(nullptr)
      {
    	  // If subentity index is equal to the number of subentities, then this intersection represents the
    	  // intersectioniterator::end(), and its other parameters are unphysical and irrelevant
    	  assert(subIndexInside_ <= 4);
    	  if (subIndexInside_ < 4)
    	  {
        	  localFaceIndex_ = gridbase.subentityLocalIndex (localIndexInside, 0, 1, subIndexInside);
        	  computeOutside();
    	  }
      }


      CurvIntersection ( const CurvIntersection &other )
  	    : localIndexInside_(other.localIndexInside_),
  	      subIndexInside_(other.subIndexInside_),
  	      gridbase_(other.gridbase_),
  	      insideGeo_(other.gridbase_->template entityGeometry<ELEMENT_CODIM>(other.localIndexInside_), *other.gridbase_),
  	      localFaceIndex_(other.localFaceIndex_),
  	      localIndexOutside_(other.localIndexOutside_),
  	      subIndexOutside_(other.subIndexOutside_),
  	      geo_(nullptr),
    	  geoInInside_(nullptr),
    	  geoInOutside_(nullptr)
      {

      }


      ~CurvIntersection()
      {
    	  if (geo_)          { delete geo_; };
    	  if (geoInInside_)  { delete geoInInside_; };
    	  if (geoInOutside_) { delete geoInOutside_; };
      }


      const CurvIntersection &operator= ( const CurvIntersection &other )
      {
    	  localIndexInside_ = other.localIndexInside_;
    	  subIndexInside_ = other.subIndexInside_;
    	  gridbase_ = other.gridbase_;
  	      localFaceIndex_ = other.localFaceIndex_;
  	      localIndexOutside_ = other.localIndexOutside_;
  	      subIndexOutside_ = other.subIndexOutside_;
    	  insideGeo_ = other.insideGeo_;
    	  geo_ = nullptr;
    	  geoInInside_ = nullptr;
    	  geoInOutside_ = nullptr;

    	  return *this;
      }

      bool equals(const CurvIntersection& other) const
      {
          return (localIndexInside_ == other.localIndexInside_) && (subIndexInside_ == other.subIndexInside_);
      }


      EntityPointer inside () const
      {
    	  //std::cout << "inside" << std::endl;

    	  IndexSetIterator thisIter = gridbase_->entityIndexIterator(ELEMENT_CODIM, localIndexInside_);
    	  StructuralType structtype = gridbase_->entityStructuralType(ELEMENT_CODIM, localIndexInside_);
    	  return EntityPointer(EntityPointerImpl(thisIter, *gridbase_, All_Partition));
      }

      EntityPointer outside () const
      {
    	  //std::cout << "outside" << std::endl;

    	  if (!neighbor())  { DUNE_THROW(Dune::IOError, "Intersection: entityPointer of non-existing outside entity requested"); }

    	  IndexSetIterator thisIter = gridbase_->entityIndexIterator(ELEMENT_CODIM, localIndexOutside_);
    	  StructuralType structtype = gridbase_->entityStructuralType(ELEMENT_CODIM, localIndexOutside_);
    	  return EntityPointer(EntityPointerImpl(thisIter, *gridbase_, All_Partition));
      }

      // By dune-convention, domain and periodic boundaries are considered boundaries
      bool boundary () const {
    	  //std::cout << "boundary" << std::endl;

    	  StructuralType structtype = gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_);
    	  return (structtype == GridStorageType::PartitionType::DomainBoundary);
      }


      /** \note Non-conformal grids not implemented atm  **/
      bool conforming () const { return true; }


      // By dune-convention, everything has a neighbour, except domain boundaries, and (process boundaries in serial case)
      bool neighbor () const
      {
    	  //std::cout << "neighbor" << std::endl;
    	  StructuralType structtype = gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_);

    	  if (structtype == GridStorageType::PartitionType::DomainBoundary)  { return false; }
    	  if (structtype == GridStorageType::PartitionType::ProcessBoundary) { return !gridbase_->isSerial(); }

    	  //std::cout << "finished neighbor" << std::endl;

    	  return true;
      }


      size_t boundarySegmentIndex () const
      {
    	  assert(boundary());  // If this entity is not a Domain Boundary should throw error
    	  return gridbase_->boundarySegmentIndex(localFaceIndex_);
      }

      // Geometry of the intersection from within the inside element
      LocalGeometry geometryInInside () const  {
    	  if (!geoInInside_)  { generateLocalGeometries(); }
    	  assert(geoInInside_);
    	  return *geoInInside_;
      }

      // Geometry of the intersection from within the outside element
      // NOTE!!! BY CURRENT CONVENTION OUTSIDE GEOMETRY HAS THE SAME ORIENTATION AS INSIDE GEOMETRY
      LocalGeometry geometryInOutside () const  {
    	  if (!neighbor())  { DUNE_THROW(Dune::IOError, "Intersection: geometry of non-existing outside entity requested"); }
    	  if (!geoInOutside_)  { generateLocalGeometries(); }
    	  assert(geoInOutside_);
    	  return *geoInOutside_;
      }

      // Geometry of this face
      Geometry geometry () const
      {
    	  //std::cout << "started geometry dim="<<dimension << " localIndexInside_="<<localIndexInside_ << " subIndexInside_=" << subIndexInside_ << std::endl;

        if(!geo_)  {
        	//InterpolatoryOrderType interporder = gridbase_->entityInterpolationOrder(ELEMENT_CODIM, localIndexInside_);
        	//std::cout << "constructing element geom" << std::endl;
        	ElementBaseGeometry entityGeometryBase = gridbase_->template entityGeometry<ELEMENT_CODIM>(localIndexInside_);

        	//std::cout << "constructing face geom" << std::endl;
        	FaceBaseGeometry faceGeometryBase = entityGeometryBase.template subentityGeometry<dimension - FACE_CODIM>(subIndexInside_);

        	//std::cout << "constructing grid-geom" << std::endl;
        	geo_ = new GeometryImpl(faceGeometryBase, *gridbase_);
        }

        //std::cout << "finished geometry" << std::endl;

        assert(geo_);
        return Geometry( *geo_ );
      }

      GeometryType type () const {
    	  if (!geo_)  { geometry(); }
    	  assert(geo_);
    	  return geo_->type();
      }

      // Face Subentity index as viewed from calling element
      InternalIndexType indexInInside () const  { return subIndexInside_; }

      // Face Subentity index as viewed from other element
      InternalIndexType indexInOutside () const  {
    	  if (!neighbor())  { DUNE_THROW(Dune::IOError, "Intersection: index of non-existing outside entity requested"); }
    	  return subIndexOutside_;
      }

      // All normals as viewed from the calling element

      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate &localFaceCoord ) const
      {
    	  GlobalCoordinate localElemCoord( geometryInInside().global( localFaceCoord ) );

    	  //std::cout << "integrationOuterNormal faceCoord=" << localFaceCoord << " elemCoord=" << localElemCoord << std::endl;

    	  return insideGeo_.subentityIntegrationNormal( indexInInside(), localElemCoord );
      }

      GlobalCoordinate outerNormal ( const LocalCoordinate &localFaceCoord ) const
      {
    	  GlobalCoordinate localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_.subentityNormal( indexInInside(), localElemCoord );
      }

      GlobalCoordinate unitOuterNormal ( const LocalCoordinate &localFaceCoord ) const
      {
    	  GlobalCoordinate localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_.subentityUnitNormal( indexInInside(), localElemCoord );
      }

      // TODO: This does not work in 2D
      GlobalCoordinate centerUnitOuterNormal () const
      {
    	  const ReferenceElement< ctype, dimension-1 > &refFace = ReferenceElements< ctype, dimension-1 >::general( type() );
    	  return unitOuterNormal( refFace.position( 0, 0 ) );
      }


      // Auxiliary methods
      // *******************************************************

      LocalIndexType intersectionIndex() const  { return localFaceIndex_; }

      // Iterates over subentities of the inside entity by increasing the inside element subentity index
      // If intersection is a ghost intersection, skip it, because the intersection is not allowed to be of type ghost
      // If next() is called beyond the allowed size, throw error
      void next()
      {
    	  //std::cout << "next" << std::endl;

    	  // When iterating over intersections, current intersection geometry becomes invalid
    	  if (geo_)           { delete geo_; }           geo_ = nullptr;
          if (geoInInside_)   { delete geoInInside_; }   geoInInside_ = nullptr;
          if (geoInOutside_)  { delete geoInOutside_; }  geoInOutside_ = nullptr;

    	  const int SUBENTITY_SIZE = 4;

    	  if (subIndexInside_ >= SUBENTITY_SIZE) { DUNE_THROW(Dune::IOError, "intersection: next() called with unexpected subentity index"); }

    	  // Increase iterator until find a non-ghost face or reach the end
    	  bool inc = true;
    	  subIndexInside_++;
    	  while( inc && (subIndexInside_ < SUBENTITY_SIZE) )
    	  {
    		  localFaceIndex_ = gridbase_->subentityLocalIndex (localIndexInside_, ELEMENT_CODIM, FACE_CODIM, subIndexInside_);

    		  //std::cout << "brrr localIndexInside_=" << localIndexInside_ << " subIndexInside_=" << subIndexInside_ << " localFaceIndex_=" << localFaceIndex_ << std::endl;

    		  inc = (gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_) == GridStorageType::PartitionType::Ghost);
    		  if (inc) { subIndexInside_++; }
    	  }

    	  // If this is not the end, update contents of the intersection
    	  if (subIndexInside_ != SUBENTITY_SIZE)  { computeOutside(); }

    	  //std::cout << "finished next" << std::endl;
      }

      // Finds outside entity local index, and this face subentity index
      void computeOutside()
      {
    	  // Only do something if the outside entity exists at all
    	  if (neighbor())
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
      }


      // [TODO] Currently assumes same geometry type on both sides
      void generateLocalGeometries() const
      {
    	  assert(!geoInInside_);
    	  assert(!geoInOutside_);

    	  //std::cout << "Intersection: started generating local geometries" << std::endl;

    	  Dune::GeometryType elemGT = insideGeo_.type();


    	  //std::cout << "aaa" << std::endl;

    	  // Find number of corners per face
    	  const Dune::ReferenceElement<ctype, dimension> & ref = Dune::ReferenceElements<ctype, dimension>::general(elemGT);
    	  int nCornerPerFace = ref.size(0, FACE_CODIM, VERTEX_CODIM);


    	  //std::cout << "bbb" << std::endl;

    	  // Generate GeometryInInside
    	  // *******************************************************8
    	  std::vector<InternalIndexType> insideCornerInternalIndex;
    	  std::vector<InternalIndexType> insideCornerLocalIndex;

		  // 1) Produce corner internal coordinates for inside and outside
		  // 2) Produce corner local coordinates

    	  //std::cout << "ccc" << std::endl;

    	  for (int i = 0; i < nCornerPerFace; i++)
    	  {
    		  insideCornerInternalIndex. push_back(ref.subEntity(subIndexInside_,  FACE_CODIM, i, VERTEX_CODIM));
    		  insideCornerLocalIndex. push_back(gridbase_->subentityLocalIndex(localIndexInside_,  ELEMENT_CODIM, VERTEX_CODIM, insideCornerInternalIndex[i]));
    	  }

    	  //std::cout << "ddd" << std::endl;

    	  // 3) Find local coordinates of outside in inside, flip internal coordinates of outside accordingly
    	  std::vector<GlobalCoordinate> inCoord;
    	  for (int i = 0; i < nCornerPerFace; i++)  { inCoord.push_back(Dune::CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, dimensionworld>(elemGT, insideCornerInternalIndex[i]));
    	  }


    	  //std::cout << "eee attempt to construct with coords=" << Dune::VectorHelper::vector2string(inCoord) << std::endl;

    	  Dune::GeometryType gt = type();

    	  //std::cout << "fff" << std::endl;

    	  geoInInside_  = new LocalGeometry(GeometryImpl(gt, inCoord,  LINEAR_ELEMENT_ORDER, *gridbase_));



    	  //std::cout << " ** done with inside, checking outside with neighbor=" << neighbor() << std::endl;

    	  if (neighbor())
    	  {
        	  // Generate GeometryInOutside
        	  // *******************************************************8
        	  std::vector<InternalIndexType> outsideCornerInternalIndex;
        	  std::vector<InternalIndexType> outsideCornerLocalIndex;

        	  // 1) Produce corner internal coordinates for inside and outside
        	  // 2) Produce corner local coordinates
        	  for (int i = 0; i < nCornerPerFace; i++)
        	  {
        		  outsideCornerInternalIndex.push_back(ref.subEntity(subIndexOutside_, FACE_CODIM, i, VERTEX_CODIM));
        		  outsideCornerLocalIndex.push_back(gridbase_->subentityLocalIndex(localIndexOutside_, ELEMENT_CODIM, VERTEX_CODIM, outsideCornerInternalIndex[i]));
        	  }

        	  // 3) Find local coordinates of outside in inside, flip internal coordinates of outside accordingly
        	  std::vector<InternalIndexType> outsideCornerInternalIndexNew;
        	  for (int iInCorner = 0; iInCorner < nCornerPerFace; iInCorner++)
        	  {
        		  for (int iOutCorner = 0; iOutCorner < nCornerPerFace; iOutCorner++)
            	  {
            		  if (insideCornerLocalIndex[iInCorner] == outsideCornerLocalIndex[iOutCorner])
            		  {
            			  outsideCornerInternalIndexNew.push_back(outsideCornerInternalIndex[iOutCorner]);
            		  }
            	  }
        	  }

        	  // Make sure we did not find too many or too few corners
        	  assert(outsideCornerInternalIndexNew.size() == nCornerPerFace);

        	  // 4) Fill internal coordinates for both geometries
        	  std::vector<GlobalCoordinate> outCoord;
        	  for (int i = 0; i < nCornerPerFace; i++)  { outCoord.push_back(Dune::CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, dimensionworld>(elemGT, outsideCornerInternalIndexNew[i])); }
        	  geoInOutside_ = new LocalGeometry(GeometryImpl(type(), outCoord, LINEAR_ELEMENT_ORDER, *gridbase_));
    	  }

    	  //std::cout << "Intersection: finished generating local geometries" << std::endl;
      }


    private:
      LocalIndexType      localFaceIndex_;
      LocalIndexType      localIndexInside_;
      LocalIndexType      localIndexOutside_;
      InternalIndexType   subIndexInside_;
      InternalIndexType   subIndexOutside_;
      ElementGeometryImpl insideGeo_;
      GridBaseType *      gridbase_;

      mutable GeometryImpl * geo_;
      mutable LocalGeometry * geoInInside_;
      mutable LocalGeometry * geoInOutside_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
