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

      typedef typename Dune::CurvGrid::CurvEntity<ELEMENT_CODIM, dimension, Grid>   EntityImpl;
      typedef typename Dune::CurvGrid::CurvEntityPointer<ELEMENT_CODIM, Grid>       EntityPointerImpl;

      typedef typename Traits::template Codim< FACE_CODIM >::GeometryImpl GeometryImpl;
      typedef typename Traits::template Codim< ELEMENT_CODIM >::GeometryImpl ElementGeometryImpl;

      typedef Dune::FieldVector< ctype, dimension-1 >      LocalCoordinate;
      typedef Dune::FieldVector< ctype, dimensionworld >   GlobalCoordinate;


      typedef typename remove_const< Grid >::type::template Codim<ELEMENT_CODIM>::EntityGeometryMappingImpl  ElementBaseGeometry;
      typedef typename remove_const< Grid >::type::template Codim<FACE_CODIM>::EntityGeometryMappingImpl     FaceBaseGeometry;

      const InterpolatoryOrderType LINEAR_ELEMENT_ORDER = 1;

    public:

      // Default constructor, as required by Dune. DO NOT USE
      CurvIntersection() { }

      CurvIntersection (
    		  LocalIndexType localIndexInside,    // Index of the element wrt which this intersection is calculated
    		  InternalIndexType subIndexInside,   // Internal index of the face as the element subentity
    		  GridBaseType & gridbase,
    		  bool ghostcheck = false             // If ghostcheck=true, then subIndex should be the first index that does not point to a ghost
      )
    	: localIndexInside_(localIndexInside),
    	  subIndexInside_(subIndexInside),
    	  gridbase_(&gridbase),
    	  geo_(nullptr),
    	  insideGeo_(nullptr),
    	  geoInInside_(nullptr),
    	  geoInOutside_(nullptr)
      {
    	  // If subentity index is equal to the number of subentities, then this intersection represents the
    	  // intersectioniterator::end(), and its other parameters are unphysical and irrelevant
    	  assert(subIndexInside_ <= 4);
    	  if (subIndexInside_ < 4)
    	  {
        	  localFaceIndex_ = gridbase.subentityLocalIndex (localIndexInside, 0, 1, subIndexInside);

        	  //std::cout << "Constructing Intersection elementIndex=" << localIndexInside << " internal index=" << subIndexInside << " gets face index=" << localFaceIndex_ << std::endl;

        	  StructuralType faceType = gridbase.entityStructuralType(FACE_CODIM, localFaceIndex_);

        	  if (faceType == GridStorageType::PartitionType::Ghost)
        	  {
        		  if (ghostcheck)  { next(); }
        		  else
        		  {
        			  std::cout << "Error: Intersection: Attempt to construct an intersection for a ghost face" << std::endl;
        			  DUNE_THROW(Dune::IOError, "Intersection: geometry of non-existing outside entity requested");
        		  }
        	  } else  { computeOutsideIndex(); }
    	  }
      }


      CurvIntersection ( const CurvIntersection &other )
  	    : localIndexInside_(other.localIndexInside_),
  	      subIndexInside_(other.subIndexInside_),
  	      gridbase_(other.gridbase_),
  	      localFaceIndex_(other.localFaceIndex_),
  	      localIndexOutside_(other.localIndexOutside_),
  	      subIndexOutside_(other.subIndexOutside_),
  	      geo_(nullptr),
  	      insideGeo_(nullptr),
    	  geoInInside_(nullptr),
    	  geoInOutside_(nullptr)
      {

      }


      ~CurvIntersection()
      {
    	  if (geo_)          { delete geo_; }
    	  if (insideGeo_)    { delete insideGeo_; }
    	  if (geoInInside_)  { delete geoInInside_; }
    	  if (geoInOutside_) { delete geoInOutside_; }
      }


      const CurvIntersection &operator= ( const CurvIntersection &other )
      {
    	  localIndexInside_ = other.localIndexInside_;
    	  subIndexInside_ = other.subIndexInside_;
    	  gridbase_ = other.gridbase_;
  	      localFaceIndex_ = other.localFaceIndex_;
  	      localIndexOutside_ = other.localIndexOutside_;
  	      subIndexOutside_ = other.subIndexOutside_;
    	  geo_ = nullptr;
    	  insideGeo_ = nullptr;
    	  geoInInside_ = nullptr;
    	  geoInOutside_ = nullptr;

    	  return *this;
      }


      bool equals(const CurvIntersection& other) const
      {
          return (localIndexInside_ == other.localIndexInside_) && (subIndexInside_ == other.subIndexInside_);
      }


      Entity inside () const
      {
    	  //std::cout << "inside" << std::endl;

    	  IndexSetIterator thisIter = gridbase_->entityIndexIterator(ELEMENT_CODIM, localIndexInside_);
    	  StructuralType structtype = gridbase_->entityStructuralType(ELEMENT_CODIM, localIndexInside_);
    	  return Entity(EntityImpl(thisIter, *gridbase_, All_Partition));
      }


      Entity outside () const
      {
    	  //std::cout << "outside" << std::endl;

    	  if (!neighbor())  {
    		  std::cout << "Error: Intersection: entityPointer of non-existing outside entity requested" << std::endl;
    		  DUNE_THROW(Dune::IOError, "Intersection: entityPointer of non-existing outside entity requested");
    	  }

    	  IndexSetIterator thisIter = gridbase_->entityIndexIterator(ELEMENT_CODIM, localIndexOutside_);
    	  StructuralType structtype = gridbase_->entityStructuralType(ELEMENT_CODIM, localIndexOutside_);
    	  return Entity(EntityImpl(thisIter, *gridbase_, All_Partition));
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
        if(!geo_)  {
        	// Compute the geometry of inside entity if it is not yet there, then use it to get the intersection geometry
        	computeInsideGeo();
        	geo_ = insideGeo_->template subentityGeometry<dimension - FACE_CODIM>(subIndexInside_);
        }

        assert(geo_);
        return Geometry( *geo_ );
      }

      GeometryType type () const {
    	  Dune::GeometryType gt;
    	  gt.makeSimplex(dimension - FACE_CODIM);
    	  return gt;
      }

      GeometryType typeInside() const {
    	  Dune::GeometryType gt;
    	  gt.makeSimplex(dimension - ELEMENT_CODIM);
    	  return gt;
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
    	  computeInsideGeo();
    	  GlobalCoordinate localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_->subentityIntegrationNormal( indexInInside(), localElemCoord );
      }

      GlobalCoordinate outerNormal ( const LocalCoordinate &localFaceCoord ) const
      {
    	  computeInsideGeo();
    	  GlobalCoordinate localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_->subentityNormal( indexInInside(), localElemCoord );
      }

      GlobalCoordinate unitOuterNormal ( const LocalCoordinate &localFaceCoord ) const
      {
    	  computeInsideGeo();
    	  GlobalCoordinate localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_->subentityUnitNormal( indexInInside(), localElemCoord );
      }

      // TODO: This does not work in 2D
      GlobalCoordinate centerUnitOuterNormal () const
      {
    	  const ReferenceElement< ctype, dimension-1 > &refFace = ReferenceElements< ctype, dimension-1 >::general( type() );
    	  return unitOuterNormal( refFace.position( 0, 0 ) );
      }


      // Auxiliary methods
      // *******************************************************

      //LocalIndexType intersectionIndex() const  { return localFaceIndex_; }


      /** \brief Iterates over subentities of the inside entity by increasing the inside element subentity index
       *   If intersection is a ghost intersection, skip it, because the intersection is not allowed to be of type ghost.
       *   Allow to go 1 index past the maximal index, thus representing the end-iterator
       * */
      void next()
      {
    	  // When iterating over intersections, current intersection geometry becomes invalid
    	  if (geo_)           { delete geo_; }           geo_ = nullptr;
    	  if (insideGeo_)     { delete insideGeo_; }     insideGeo_ = nullptr;
          if (geoInInside_)   { delete geoInInside_; }   geoInInside_ = nullptr;
          if (geoInOutside_)  { delete geoInOutside_; }  geoInOutside_ = nullptr;

    	  const int SUBENTITY_SIZE = 4;

    	  if (subIndexInside_ >= SUBENTITY_SIZE) {
    		  std::cout << "Error: Intersection: next() called with unexpected subentity index" << std::endl;
    		  DUNE_THROW(Dune::IOError, "intersection: next() called with unexpected subentity index");
    	  }

    	  // Increase iterator until find a non-ghost face or reach the end
    	  bool inc = true;
    	  subIndexInside_++;
    	  while( inc && (subIndexInside_ < SUBENTITY_SIZE) )
    	  {
    		  localFaceIndex_ = gridbase_->subentityLocalIndex (localIndexInside_, ELEMENT_CODIM, FACE_CODIM, subIndexInside_);

    		  inc = (gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_) == GridStorageType::PartitionType::Ghost);
    		  if (inc) { subIndexInside_++; }
    	  }

    	  // If this is not the end, update contents of the intersection
    	  if (subIndexInside_ != SUBENTITY_SIZE)  { computeOutsideIndex(); }
      }


      /** \brief Calculates inside geometry */
      void computeInsideGeo() const {
    	  if (!insideGeo_)  { insideGeo_ = new ElementGeometryImpl(gridbase_->template entityGeometry<ELEMENT_CODIM>(localIndexInside_), *gridbase_); }
      }


      /** \brief Finds outside entity local index, and this face subentity index as seen from outside */
      void computeOutsideIndex()
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
        		  if (iFace >= 4)  {
        			  std::cout << "*** when searching faceIndex=" << localFaceIndex_ << " of structural type" << gridbase_->entityStructuralType(FACE_CODIM, localFaceIndex_) <<" of outsideIndex=" << localIndexOutside_ << std::endl;
        			  std::cout << "Error: Intersection: Not found face as its outside-element subentity" << std::endl;
        			  DUNE_THROW(Dune::IOError, "Intersection: Not found face as its outside-element subentity");
        		  }

        		  if (localFaceIndex_ == gridbase_->subentityLocalIndex(localIndexOutside_, ELEMENT_CODIM, FACE_CODIM, iFace))
        		  {
        			  found_face = true;
        			  subIndexOutside_ = iFace;
        		  }
        		  iFace++;
        	  }
    	  }
      }


      /** \brief Finds outside entity local index, and this face subentity index as seen from outside
       * [TODO] Currently assumes same geometry type on both sides
       * */
      void generateLocalGeometries() const
      {
    	  assert(!geoInInside_);
    	  assert(!geoInOutside_);

    	  // Find number of corners per face
    	  Dune::GeometryType elemGT = typeInside();
    	  const Dune::ReferenceElement<ctype, dimension> & ref = Dune::ReferenceElements<ctype, dimension>::general(elemGT);
    	  int nCornerPerFace = ref.size(0, FACE_CODIM, VERTEX_CODIM);

    	  // Generate GeometryInInside
    	  // *******************************************************8
    	  std::vector<InternalIndexType> insideCornerInternalIndex;
    	  std::vector<InternalIndexType> insideCornerLocalIndex;

		  // 1) Produce corner internal coordinates for inside and outside
		  // 2) Produce corner local coordinates
    	  for (int i = 0; i < nCornerPerFace; i++)
    	  {
    		  insideCornerInternalIndex. push_back(ref.subEntity(subIndexInside_,  FACE_CODIM, i, VERTEX_CODIM));
    		  insideCornerLocalIndex. push_back(gridbase_->subentityLocalIndex(localIndexInside_,  ELEMENT_CODIM, VERTEX_CODIM, insideCornerInternalIndex[i]));
    	  }

    	  // 3) Find local coordinates of outside in inside, flip internal coordinates of outside accordingly
    	  std::vector<GlobalCoordinate> inCoord;
    	  for (int i = 0; i < nCornerPerFace; i++)  { inCoord.push_back(Dune::CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, dimensionworld>(elemGT, insideCornerInternalIndex[i]));  }

    	  geoInInside_  = new LocalGeometry(GeometryImpl(type(), inCoord,  LINEAR_ELEMENT_ORDER, *gridbase_));

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
      }


    private:
      LocalIndexType      localFaceIndex_;
      LocalIndexType      localIndexInside_;
      LocalIndexType      localIndexOutside_;
      InternalIndexType   subIndexInside_;
      InternalIndexType   subIndexOutside_;

      GridBaseType *      gridbase_;

      // Store all geometries as pointers and init only when needed to accelerate the code
      mutable GeometryImpl * geo_;
      mutable ElementGeometryImpl * insideGeo_;
      mutable LocalGeometry * geoInInside_;
      mutable LocalGeometry * geoInOutside_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
