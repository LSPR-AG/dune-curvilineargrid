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
    class CurvIntersection
    {

      typedef typename Grid::Traits Traits;

    public:
      typedef typename Grid::ctype ctype;

      static const int dimension = Grid::dimension;
      static const int dimensionworld = Grid::dimensionworld;

	  typedef typename Grid::GridStorageType  GridStorageType;
	  typedef typename Grid::GridBaseType     GridBaseType;

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
      //typedef typename Traits::template Codim< ELEMENT_CODIM >::EntityPointer EntityPointer;
      typedef typename Traits::template Codim< FACE_CODIM >::Geometry Geometry;
      typedef typename Traits::template Codim< FACE_CODIM >::LocalGeometry LocalGeometry;

      typedef typename Traits::template Codim< ELEMENT_CODIM >::Geometry ElementGeometry;

    private:

      typedef typename Dune::CurvGrid::CurvEntity<ELEMENT_CODIM, dimension, Grid>   EntityImpl;
      //typedef typename Dune::CurvGrid::CurvEntityPointer<ELEMENT_CODIM, Grid>       EntityPointerImpl;

      typedef typename Traits::template Codim< FACE_CODIM >::GeometryImpl GeometryImpl;
      typedef typename Traits::template Codim< ELEMENT_CODIM >::GeometryImpl ElementGeometryImpl;

      typedef Dune::FieldVector< ctype, dimension-1 >      LocalCoordinateFace;
      typedef Dune::FieldVector< ctype, dimension   >      LocalCoordinateElement;
      typedef Dune::FieldVector< ctype, dimensionworld >   GlobalCoordinate;


      typedef typename Grid::template Codim<ELEMENT_CODIM>::EntityGeometryMappingImpl  ElementBaseGeometry;
      typedef typename Grid::template Codim<FACE_CODIM>::EntityGeometryMappingImpl     FaceBaseGeometry;

      typedef typename ElementBaseGeometry::JacobianInverseTransposed   JacobianInverseTransposed;

      const InterpolatoryOrderType LINEAR_ELEMENT_ORDER = 1;

    public:

      // Default constructor, as required by Dune. DO NOT USE
      CurvIntersection() { }

      CurvIntersection (
    		  LocalIndexType localIndexInside,    // Index of the element wrt which this intersection is calculated
    		  InternalIndexType subIndexInside,   // Internal index of the face as the element subentity
    		  GridBaseType & gridbase,
    		  bool ghostcheck = false             // If ghostcheck=true, then subIndex should be the first index that does not point to a ghost
      ) :
    	  localIndexInside_(localIndexInside),
    	  subIndexInside_(subIndexInside),
    	  gridbase_(gridbase),
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
        	  localFaceIndex_ = gridbase_.entity().subentityLocalIndex (localIndexInside, ELEMENT_CODIM, FACE_CODIM, subIndexInside);

        	  //std::cout << "Constructing Intersection elementIndex=" << localIndexInside << " internal index=" << subIndexInside << " gets face index=" << localFaceIndex_ << std::endl;

        	  StructuralType faceType = gridbase_.entity().partitionType(FACE_CODIM, localFaceIndex_);

        	  if (faceType == Dune::PartitionType::GhostEntity)
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


      CurvIntersection ( const CurvIntersection &other ) :
    	  localFaceIndex_(other.localFaceIndex_),
  	      localIndexInside_(other.localIndexInside_),
  	      localIndexOutside_(other.localIndexOutside_),
  	      subIndexInside_(other.subIndexInside_),
  	      subIndexOutside_(other.subIndexOutside_),
  	      gridbase_(other.gridbase_),
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
    	  IndexSetIterator thisIter = gridbase_.indexset().entityIndexSetSelect(ELEMENT_CODIM).find(localIndexInside_);
    	  return Entity(EntityImpl(thisIter, gridbase_, All_Partition));
      }


      Entity outside () const
      {
    	  if (!neighbor())  {
    		  std::cout << "Error: Intersection: entityPointer of non-existing outside entity requested" << std::endl;
    		  DUNE_THROW(Dune::IOError, "Intersection: entityPointer of non-existing outside entity requested");
    	  }

    	  IndexSetIterator thisIter = gridbase_.indexset().entityIndexSetSelect(ELEMENT_CODIM).find(localIndexOutside_);
    	  return Entity(EntityImpl(thisIter, gridbase_, All_Partition));
      }


      // By dune-convention, domain and periodic boundaries are considered boundaries
      bool boundary () const {
    	  StructuralType boundaryType = gridbase_.intersection().boundaryType(localFaceIndex_);
    	  bool isDB = (boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary);
    	  bool isPeriodic = (boundaryType == GridStorageType::FaceBoundaryType::PeriodicBoundary);
    	  return (isDB || isPeriodic);
      }


      /** \note Non-conformal grids not implemented atm  **/
      bool conforming () const { return true; }


      // By dune-convention, everything has a neighbour, except domain boundaries, and (process boundaries in serial case)
      bool neighbor () const
      {
    	  StructuralType boundaryType  = gridbase_.intersection().boundaryType(localFaceIndex_);
    	  StructuralType partitionType = gridbase_.entity().partitionType(FACE_CODIM, localFaceIndex_);

    	  bool isDB = (boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary);
    	  bool isPB = (partitionType == Dune::PartitionType::BorderEntity);
    	  bool isPeriodic = (boundaryType == GridStorageType::FaceBoundaryType::PeriodicBoundary);

    	  // Check that no more than 1 one of the above is true
    	  int nSpecialProperty = int(isDB) + int(isPB) + int(isPeriodic);
    	  if (nSpecialProperty > 1) {
    		  std::stringstream logstr;
    		  logstr << "Intersection is checked for (isDB, isPB, isPeriodic) = (" << isDB << ", " << isPB << ", "
					  << isPeriodic << "). Having two properties at the same time is unexpected ";
    		  DUNE_THROW(Dune::IOError, logstr.str() );
    	  }

    	  if			(isDB)			{ return false; }	// Pure domain boundaries should not have outer neighbors
    	  else if	(isPB)			{ return gridbase_.property().withGhostElements(); }  // Boundaries that can have ghost neighbor must check if ghosts are present
    	  else if	(isPeriodic)	{ return gridbase_.property().withGhostElements() || gridbase_.intersection().checkOuterNeighbor(localFaceIndex_); }		// Periodic boundaries may be interior or ghost
    	  else {
    		  // Interior case - must always have neighbor
    		  assert(gridbase_.intersection().checkOuterNeighbor(localFaceIndex_));
    		  return true;
    	  }
      }


      size_t boundarySegmentIndex () const
      {
    	  assert(boundary());  // If this entity is not a Domain Boundary should throw error
    	  return gridbase_.intersection().boundarySegmentIndex(localFaceIndex_);
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
        	geo_ = new GeometryImpl(
        		insideGeo_->basegeometry().template subentityGeometry<dimension - FACE_CODIM>(subIndexInside_),
        		gridbase_
        		);
        }

        assert(geo_);
        return Geometry( *geo_ );
      }

      GeometryType type () const {
    	  Dune::GeometryType gt=Dune::GeometryTypes::simplex(dimension - FACE_CODIM);
    	  return gt;
      }

      GeometryType typeInside() const {
    	  Dune::GeometryType gt=Dune::GeometryTypes::simplex(dimension - FACE_CODIM);
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

      GlobalCoordinate outerNormal            (const LocalCoordinateFace &localFaceCoord) const { return outerNormal(localFaceCoord, nullptr); }
      GlobalCoordinate unitOuterNormal        (const LocalCoordinateFace &localFaceCoord) const { return unitOuterNormal(localFaceCoord, nullptr); }
      GlobalCoordinate integrationOuterNormal (const LocalCoordinateFace &localFaceCoord) const { return integrationOuterNormal(localFaceCoord, nullptr); }

      GlobalCoordinate outerNormal ( const LocalCoordinateFace &localFaceCoord, JacobianInverseTransposed * thisjit ) const
      {
    	  computeInsideGeo();
    	  LocalCoordinateElement localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_->basegeometry().subentityNormal( indexInInside(), localElemCoord, thisjit );
      }

      GlobalCoordinate unitOuterNormal ( const LocalCoordinateFace &localFaceCoord, JacobianInverseTransposed * thisjit ) const
      {
    	  computeInsideGeo();
    	  LocalCoordinateElement localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_->basegeometry().subentityUnitNormal( indexInInside(), localElemCoord, thisjit );
      }

      GlobalCoordinate integrationOuterNormal ( const LocalCoordinateFace &localFaceCoord, JacobianInverseTransposed * thisjit) const
      {
    	  computeInsideGeo();
    	  LocalCoordinateElement localElemCoord( geometryInInside().global( localFaceCoord ) );
    	  return insideGeo_->basegeometry().subentityIntegrationNormal( indexInInside(), localElemCoord, thisjit);
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
    		  localFaceIndex_ = gridbase_.entity().subentityLocalIndex (localIndexInside_, ELEMENT_CODIM, FACE_CODIM, subIndexInside_);

    		  inc = (gridbase_.entity().partitionType(FACE_CODIM, localFaceIndex_) == Dune::PartitionType::GhostEntity);
    		  if (inc) { subIndexInside_++; }
    	  }

    	  // If this is not the end, update contents of the intersection
    	  if (subIndexInside_ != SUBENTITY_SIZE)  { computeOutsideIndex(); }
      }


      /** \brief Calculates inside geometry */
      void computeInsideGeo() const {
    	  if (!insideGeo_)  { insideGeo_ = new ElementGeometryImpl(gridbase_.entity().template geometry<ELEMENT_CODIM>(localIndexInside_), gridbase_); }
      }


      /** \brief Finds outside entity local index, and this face subentity index as seen from outside */
      void computeOutsideIndex()
      {
    	  // Only do something if the outside entity exists at all
    	  if (neighbor())
    	  {
        	  LocalIndexType neighborElementIndex0 = gridbase_.intersection().neighborElement(localFaceIndex_, 0);
        	  LocalIndexType neighborElementIndex1 = gridbase_.intersection().neighborElement(localFaceIndex_, 1);

        	  // Determine the outside element local index
        	  // NOTE: If the outside element is not interior, then the inside elem is always 0th
        	  // But if both neighbors are interior, then the order is fixed but non-trivial. Find order by comparison
        	  bool thisInsideElemIs0 = neighborElementIndex0 == localIndexInside_;
        	  bool thisInsideElemIs1 = neighborElementIndex1 == localIndexInside_;
        	  assert(thisInsideElemIs0 != thisInsideElemIs1);  // Exactly one of them should be the inside element
        	  localIndexOutside_ = thisInsideElemIs0 ? neighborElementIndex1 : neighborElementIndex0;
        	  InternalIndexType neighborOrderOutside = thisInsideElemIs0 ? 1 : 0;

        	  // IMPORTANT NOTE: IN CASE OF PERIODIC FACES, THE INTERSECTION AS SEEN FROM INSIDE AND OUTSIDE ARE DIFFERENT ENTITIES
        	  subIndexOutside_ = gridbase_.intersection().subIndexInNeighbor(localFaceIndex_, neighborOrderOutside);


//        	  int iFace = 0;
//        	  bool found_face = false;
//        	  while (!found_face)
//        	  {
//        		  if (iFace >= 4)  {
//        			  std::stringstream logstr;
//        			  logstr << "Error: Intersection: Not found face as its outside-element subentity" << std::endl;
//        			  logstr << "*** when searching faceIndex=" << localFaceIndex_;
//        			  logstr << " of structural type" << gridbase_.entityPartitionType(FACE_CODIM, localFaceIndex_);
//        			  logstr <<" of outsideIndex=" << localIndexOutside_;
//
//        			  DUNE_THROW(Dune::IOError, logstr.str());
//        		  }
//
//        		  if (localFaceIndex_ == gridbase_.entity().subentityLocalIndex(localIndexOutside_, ELEMENT_CODIM, FACE_CODIM, iFace))
//        		  {
//        			  found_face = true;
//        			  subIndexOutside_ = iFace;
//        		  }
//        		  iFace++;
//        	  }
    	  }
      }


      /** \brief Finds outside entity local index, and this face subentity index as seen from outside
       * [TODO] Currently assumes same geometry type on both sides
       *
       * Theory:
       * 1) Let L0, L1 and L2 are the internal indices of tetrahedron vertices, that exactly correspond to the face vertices p0, p1 and p2
       * 2) Let R0, R1 and R2 are the reference tetrahedron vertices, exactly corresponding to L0, L1, L2
       * 3) Assume an affine map M(u,v) = R0 + u(R1-R0) + v(R2-R0)
       * 4) Then it can be shown that M(triangle local coordinate) = tetrahedron local coordinate of the same point
       *
       * Algorithm:
       * 1) Find L0, L1, L2 by direct comparison of local indices of vertices of triangle and tetrahedron
       * 2) Find R0, R1, R2 from the reference element
       * 3) Store geometry of a linear triangle T(R0, R1, R2)
       * 4) Then T.global(triangle_local) = tetrahedron_local
       *
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
    	  std::vector<InternalIndexType> insideCornerIndexInParent;
    	  std::vector<InternalIndexType> insideCornerLocalIndex;

		  // 1) Produce corner internal coordinates for inside and outside
		  // 2) Produce corner local coordinates
    	  for (int i = 0; i < nCornerPerFace; i++)
    	  {
    		  insideCornerIndexInParent. push_back(ref.subEntity(subIndexInside_,  FACE_CODIM, i, VERTEX_CODIM));
    		  insideCornerLocalIndex. push_back(gridbase_.entity().subentityLocalIndex(localIndexInside_,  ELEMENT_CODIM, VERTEX_CODIM, insideCornerIndexInParent[i]));
    	  }

    	  // 3) Find local coordinates of outside in inside, flip internal coordinates of outside accordingly
    	  std::vector<GlobalCoordinate> inCoord;
    	  for (int i = 0; i < nCornerPerFace; i++)  { inCoord.push_back(CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, dimensionworld>(elemGT, insideCornerIndexInParent[i]));  }

    	  geoInInside_  = new LocalGeometry(GeometryImpl(type(), inCoord,  LINEAR_ELEMENT_ORDER, gridbase_));

    	  if (neighbor())
    	  {
        	  // Generate GeometryInOutside
        	  // *******************************************************8
        	  std::vector<InternalIndexType> outsideCornerIndexInParent;
        	  std::vector<InternalIndexType> outsideCornerLocalIndex;

        	  // 1) Produce corner internal coordinates for inside and outside
        	  // 2) Produce corner local coordinates
        	  for (int i = 0; i < nCornerPerFace; i++)
        	  {
        		  outsideCornerIndexInParent.push_back(ref.subEntity(subIndexOutside_, FACE_CODIM, i, VERTEX_CODIM));
        		  outsideCornerLocalIndex.push_back(gridbase_.entity().subentityLocalIndex(localIndexOutside_, ELEMENT_CODIM, VERTEX_CODIM, outsideCornerIndexInParent[i]));
        	  }

        	  std::vector<InternalIndexType> outsideCornerIndexInParentNew;
        	  std::vector<GlobalCoordinate> outCoord;

        	  StructuralType boundaryType  = gridbase_.intersection().boundaryType(localFaceIndex_);
        	  if (boundaryType != GridStorageType::FaceBoundaryType::PeriodicBoundary) {

        		  // If the two neighboring elements are not periodic, they share a common face.
        		  // Can determine the correct rotation of the geometryInOutside by matching the corner local indices
            	  for (int iInCorner = 0; iInCorner < nCornerPerFace; iInCorner++)
            	  {
            		  for (int iOutCorner = 0; iOutCorner < nCornerPerFace; iOutCorner++)
                	  {
                		  if (insideCornerLocalIndex[iInCorner] == outsideCornerLocalIndex[iOutCorner])
                		  {
                			  outsideCornerIndexInParentNew.push_back(outsideCornerIndexInParent[iOutCorner]);
                		  }
                	  }
            	  }

            	  // Make sure we did not find too many or too few corners
            	  assert(outsideCornerIndexInParentNew.size() == nCornerPerFace);

        	  } else {

        		  // In case of periodic boundaries, they are disjoint, and thus additional information is required to determine the correct rotation
        		  // A permutation index is constructed during the periodic face construction, so we will use it
        		  unsigned int permutationIndexInner = gridbase_.intersection().periodicPermutationInner(localFaceIndex_);
        		  unsigned int permutationIndexOuter = gridbase_.intersection().periodicPermutationOuter(localFaceIndex_);

        		  outsideCornerIndexInParentNew = CurvilinearGeometryHelper::permuteVec(outsideCornerIndexInParent, permutationIndexInner);
        	  }

        	  // 4) Fill internal coordinates for both geometries
        	  for (int i = 0; i < nCornerPerFace; i++)  { outCoord.push_back(CurvilinearGeometryHelper::cornerInternalCoordinate<ctype, dimensionworld>(elemGT, outsideCornerIndexInParentNew[i])); }

        	  geoInOutside_ = new LocalGeometry(GeometryImpl(type(), outCoord, LINEAR_ELEMENT_ORDER, gridbase_));
    	  }
      }


    private:
      LocalIndexType      localFaceIndex_;
      LocalIndexType      localIndexInside_;
      LocalIndexType      localIndexOutside_;
      InternalIndexType   subIndexInside_;
      InternalIndexType   subIndexOutside_;

      GridBaseType & gridbase_;

      // Store all geometries as pointers and init only when needed to accelerate the code
      mutable GeometryImpl * geo_;
      mutable ElementGeometryImpl * insideGeo_;
      mutable LocalGeometry * geoInInside_;
      mutable LocalGeometry * geoInOutside_;
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_INTERSECTION_HH
