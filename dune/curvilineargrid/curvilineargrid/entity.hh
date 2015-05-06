// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVILINEARGRID_ENTITY_HH
#define DUNE_CURVILINEARGRID_ENTITY_HH

#include <dune/common/iteratorrange.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>


namespace Dune
{

  namespace CurvGrid
  {


  // Forwards-Declaration
  template<int codim, class Grid>
  class CurvEntitySeed;

  template<class Grid >
  class CurvHierarchicIterator;

  template<class Grid>
  class CurvIntersectionIterator;







  template<int codim, int dim, class Grid>
  class CurvEntityBase
  {
  public:
	  typedef typename remove_const< Grid >::type::Traits Traits;

	  typedef typename remove_const< Grid >::type::ctype ctype;

	  typedef typename Traits::template Codim< codim >::EntitySeed  EntitySeed;
	  typedef Dune::CurvGrid::CurvEntitySeed<codim, Grid>           EntitySeedImpl;

  public:

	  static const int codimension     = codim;				                         //! codimensioon of the entity
	  static const int dimension       = remove_const< Grid >::type::dimension;		 //! dimension of the grid
	  static const int mydimension     = dimension - codimension;		                 //! dimension of the entity
	  static const int dimensionworld  = remove_const< Grid >::type::dimensionworld;   //! dimension of the world


	  typedef typename Traits::template Codim< codimension >::Geometry      Geometry;	    //! type of corresponding geometry
	  typedef typename Traits::template Codim< codimension >::GeometryImpl  GeometryImpl;

	  typedef typename remove_const< Grid >::type::GridStorageType  GridStorageType;
	  typedef typename remove_const< Grid >::type::GridBaseType     GridBaseType;
	  typedef typename GridStorageType::IdType                  IdType;
	  typedef typename GridBaseType::LocalIndexType             LocalIndexType;
	  typedef typename GridBaseType::StructuralType             StructuralType;

	  typedef typename GridBaseType::IndexSetIterator           IndexSetIterator;

  public:
		/** \name Construction, Initialization and Destruction
		*  \{ */

	  //! Fake constructor
	  CurvEntityBase ()
        : gridbase_(nullptr),
          geometry_(nullptr)
      { }

	  // Constructor that allows the iterator class to iterate over the entities
	  CurvEntityBase (
	    const IndexSetIterator & iter,
	    GridBaseType & gridbase,
	    Dune::PartitionIteratorType pitype
	    )
	  	  :
	  		gridbaseIndexIterator_(iter),
	  	    gridbase_(&gridbase),
	  	    pitype_(pitype),
	        geometry_(nullptr)
	  {  }


	  // Constructor for when the iterator is irrelevant. Creates default iterator over all-partition
	  CurvEntityBase (
	    LocalIndexType localIndex,
	    GridBaseType & gridbase,
	    Dune::PartitionIteratorType pitype
	    )
	  	  :
	  		gridbaseIndexIterator_(gridbase.entityIndexIterator(codim, localIndex)),
	  	    gridbase_(&gridbase),
	  	    pitype_(pitype),
	        geometry_(nullptr)
	  {  }


	  //! Copy constructor from an existing entity.
	  CurvEntityBase(const CurvEntityBase& other) :
		  gridbaseIndexIterator_ (other.gridbaseIndexIterator_),
		  gridbase_ (other.gridbase_),
		  pitype_ (other.pitype_),
          geometry_(nullptr)
	  {
	  }

	  //! Move constructor from an existing entity.
	  //CurvEntityBase(CurvEntityBase&& other) : realEntity(std::move(other.realEntity)) {}


	  //! Copy assignment operator from an existing entity.
	  CurvEntityBase& operator=(const CurvEntityBase& other)
	  {
	      gridbaseIndexIterator_ = other.gridbaseIndexIterator_;
	      gridbase_ = other.gridbase_;
	      pitype_ = other.pitype_;
	      geometry_ = nullptr;
	      return *this;
	  }

	  ~CurvEntityBase()
	  {
		  if (geometry_)  { delete geometry_; }
	  }


	  //! Move assignment operator from an existing entity.
	  //CurvEntityBase& operator=(CurvEntityBase&& other)  { realEntity = std::move(other.realEntity);  return *this; }

	  /** \} */


	  /** \brief compare two entities
	   *
	   *  Only the codimension and local index of the entity matter, all other constructions are auxiliary
	   *
	   * */
	  bool equals ( const CurvEntityBase &other) const
	  {
		  return (*gridbaseIndexIterator_ == *other.gridbaseIndexIterator_);
	  }


      /** \name Methods Shared by Entities of All Codimensions
	   *  \{ */

	  /** \brief Return the name of the reference element. The type can be used to access the Dune::ReferenceElement. */
	  GeometryType type () const { return gridbase_->entityGeometryType(codim, *gridbaseIndexIterator_); }

	  /** \brief  Returns the (refinement) level of this entity */
	  int level () const { return gridbase_->entityLevel(codim, *gridbaseIndexIterator_); }

	  /** \brief obtain the partition type of this entity */
	  PartitionType partitionType () const  {
	  	return gridbase_->entityPartitionType(codim, *gridbaseIndexIterator_);
      }


	  /** \brief obtain geometric realization of the entity */
	  Geometry geometry () const {
		  if (!geometry_) {
			  geometry_ = new Geometry(GeometryImpl(gridbase_->template entityGeometry<codim>(*gridbaseIndexIterator_), *gridbase_));
		  }
		  return *geometry_;
	  }


	  /** \brief obtain the entity's index */
	  int index () const  {
		  LocalIndexType thisIndex = indexBase();

		  if (codim == dimension)     { return gridbase_->cornerUniqueLocalIndex(thisIndex); }
		  else                        { return thisIndex; }
	  }


	  /** \brief obtain the index of a subentity from a host IndexSet
	   *
	   *  \param[in]  i         number of the subentity
	   *  \param[in]  codim        codimension of the subentity
	   */
	  int subIndex (int internalIndex, unsigned int subcodim ) const  {
	      LocalIndexType thisIndex = subIndexBase(internalIndex, subcodim);

		  if (subcodim == dimension)  { return gridbase_->cornerUniqueLocalIndex(thisIndex); }
		  else                        { return thisIndex; }
	  }


	  /** \brief obtain the entity's index */
	  int indexBase () const  { return *gridbaseIndexIterator_; }

	  int subIndexBase(int internalIndex, unsigned int subcodim) const  {
		  return gridbase_->subentityLocalIndex(*gridbaseIndexIterator_, codimension, subcodim, internalIndex);
	  }



	  /** \brief obtain the entity's id from a host IdSet */
	  IdType id () const  { return gridbase_->globalId(codim, *gridbaseIndexIterator_); }


	  IdType subId ( int internalIndex, unsigned int subcodim ) const
	  {
	      int subentityLocalIndex = subIndex(internalIndex, subcodim);
	      return gridbase_->globalId(subcodim, subentityLocalIndex);
	  }


	  /** \} */


	  /** \brief Return the entity seed which contains sufficient information to generate the entity again and uses as little memory as possible */
	  EntitySeed seed () const  { return EntitySeed(EntitySeedImpl(*gridbaseIndexIterator_, pitype_)); }


	  /** \brief moves to the next entity within the base storage. Additional functionality used by iterators */
	  void next()  {
		  gridbaseIndexIterator_++;
		  if (geometry_)  { delete geometry_;    geometry_ = nullptr; }
	  }

  protected:
	    IndexSetIterator gridbaseIndexIterator_;
	    mutable Geometry * geometry_;
	    GridBaseType * gridbase_;
	    Dune::PartitionIteratorType pitype_;
  };



  /** Generic Entity class valid for all codim */
  template<int codim, int dim, class Grid>
  class CurvEntity : public CurvEntityBase<codim, dim, Grid>
  {
	    typedef CurvEntityBase<codim, dim, Grid>   Base;

	    typedef typename Base::IndexSetIterator  IndexSetIterator;
	    typedef typename Base::GridBaseType      GridBaseType;

	    typedef typename GridBaseType::LocalIndexType             LocalIndexType;

  public:

	  //! Fake constructor
	  CurvEntity () : Base()  { }

      CurvEntity (
	      const IndexSetIterator & iter,
	      GridBaseType & gridbase,
	      Dune::PartitionIteratorType pitype)
	  	      : Base(iter, gridbase, pitype)
	  {}

	  // Constructor for when the iterator is irrelevant. Creates default iterator over all-partition
	  CurvEntity (
          LocalIndexType localIndex,
	      GridBaseType & gridbase,
	      Dune::PartitionIteratorType pitype)
	  	      :  Base(localIndex, gridbase, pitype)
	  {  }

	  // Constructor from a seed
	  template <class Seed>
	  CurvEntity(Seed & seed, GridBaseType & gridbase)
	     : Base(seed.localIndex(), gridbase, seed.partitionIteratorType())
	  {

	  }

  };



  /** Specialisation of the entity class for elements (codim=0) */
  template<int dim, class Grid>
  class CurvEntity <0, dim, Grid> : public CurvEntityBase<0, dim, Grid>
  {
	  typedef typename remove_const< Grid >::type::Traits  Traits;
	  typedef typename remove_const< Grid >::type::ctype   ctype;						//! coordinate type of the grid

	  typedef  CurvEntity <0, dim, Grid>  This;

  public:

	  /** \name Attributes
	   *  \{ */

	  static const int codimension    = 0;		                                //! codimensioon of the entity
	  static const int dimension      = remove_const< Grid >::type::dimension;	//! dimension of the grid
	  static const int mydimension    = dimension;                              //! dimension of the entity
	  static const int dimensionworld = remove_const< Grid >::type::dimensionworld;                                    //! dimension of the world
	  /** \} */


	  typedef typename Traits::template Codim< 0 >::Entity   Entity;
	  /** \brief The geometry type of this entity */
	  typedef typename Traits::template Codim< 0 >::Geometry Geometry;

	  typedef typename Traits::template Codim< codimension >::GeometryImpl GeometryImpl;

	  //! \brief The corresponding entity seed (for storage of entities)
	  typedef typename Traits::template Codim< 0 >::EntitySeed EntitySeed;

	  /** \brief The geometry type of this entity when the geometry is expressed embedded in the father element. */
	  typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

	  /** \brief The HierarchicIterator type*/
	  typedef typename Traits::HierarchicIterator               HierarchicIterator;
	  typedef Dune::CurvGrid::CurvHierarchicIterator<Grid>   HierarchicIteratorImpl;


	  typedef Dune::CurvGrid::CurvIntersectionIterator<Grid>   IntersectionIteratorImpl;


	  typedef typename remove_const< Grid >::type::GridStorageType  GridStorageType;
	  typedef typename remove_const< Grid >::type::GridBaseType     GridBaseType;

	  typedef typename GridBaseType::GlobalIndexType           GlobalIndexType;
	  typedef typename GridBaseType::LocalIndexType            LocalIndexType;
	  typedef typename GridBaseType::InternalIndexType         InternalIndexType;
	  typedef typename GridBaseType::StructuralType            StructuralType;
	  typedef typename GridBaseType::PhysicalTagType           PhysicalTagType;
	  typedef typename GridBaseType::InterpolatoryOrderType    InterpolatoryOrderType;

      // Codimensions of entity types for better code readability
      static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
      static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
      static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
      static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;


	  typedef CurvEntityBase<0, dim, Grid>             Base;
	  typedef typename Base::IndexSetIterator          IndexSetIterator;

	  using Base::gridbaseIndexIterator_;
	  using Base::pitype_;
	  using Base::gridbase_;

  public:

	  //! Fake constructor
	  CurvEntity () : Base()  { }

	  CurvEntity (
	      const IndexSetIterator & iter,
	      GridBaseType & gridbase,
	      Dune::PartitionIteratorType pitype)
	  	      : Base(iter, gridbase, pitype)
	  {}


	  // Constructor for when the iterator is irrelevant. Creates default iterator over all-partition
	  CurvEntity (
          LocalIndexType localIndex,
	      GridBaseType & gridbase,
	      Dune::PartitionIteratorType pitype)
	  	      :  Base(localIndex, gridbase, pitype)
	  {  }

	  // Constructor from a seed
	  template <class Seed>
	  CurvEntity(Seed & seed, GridBaseType & gridbase)
	     : Base(seed.localIndex(), gridbase, seed.partitionIteratorType())
	  {

	  }


   /**\brief Number of subentities with codimension <tt>codim</tt>.
     *
     * Strictly speaking this method is redundant, because the same information can be obtained
     * from the corresponding reference element. It is here for efficiency reasons only.
     */
    unsigned int subEntities(unsigned int codim) const
    {
    	LocalIndexType entityLocalIndex = *gridbaseIndexIterator_;
    	Dune::GeometryType gt = gridbase_->entityGeometryType(ELEMENT_CODIM, entityLocalIndex);

        return Dune::ReferenceElements<double, dimension>::general(gt).size(codim);
    }

    template<int codim>
    int count () const  { return subEntities(codim); }


    /** \brief Obtain a pointer to a subentity
     *
     *  \tparam  codim  codimension of the desired subentity
     *
     *  \param[in]  i  number of the subentity (in generic numbering)
     *
     *  \returns an EntityPointer to the specified subentity
     *
     *  \note The subentities are numbered 0, ..., count< codim >-1
     */
    template< int subcodim >
    typename Traits::template Codim< subcodim >::Entity
    subEntity ( int i ) const
    {
    	int subLocalIndex = gridbase_->subentityLocalIndex(*gridbaseIndexIterator_, 0, subcodim, i);
    	IndexSetIterator subIterator = gridbase_->entityIndexDuneIterator(subcodim, All_Partition, subLocalIndex);
    	return typename Traits::template Codim< subcodim >::Entity(CurvEntity<subcodim, dimension, Grid>(subIterator, *gridbase_, pitype_));
    }


    typename Traits::LevelIntersectionIterator ilevelbegin () const
    {
  	    InternalIndexType firstFaceSubIndex = 0;
  	    LocalIndexType elementLocalIndex = *gridbaseIndexIterator_;
  	    LocalIndexType faceLocalIndex = gridbase_->subentityLocalIndex(elementLocalIndex, ELEMENT_CODIM, FACE_CODIM, firstFaceSubIndex);

  	    IntersectionIteratorImpl iter (elementLocalIndex, firstFaceSubIndex, *gridbase_, true);

        return iter;
    }

    typename Traits::LevelIntersectionIterator ilevelend () const
    {
    	// Generate a face subentity index which is +1 to total number of faces, such that it is just 1 above all faces
    	LocalIndexType entityLocalIndex = *gridbaseIndexIterator_;
    	Dune::GeometryType gt = gridbase_->entityGeometryType(ELEMENT_CODIM, entityLocalIndex);
  	    InternalIndexType nSubentityFace = Dune::ReferenceElements<ctype, dimension>::general(gt).size(FACE_CODIM);
  	    return IntersectionIteratorImpl(entityLocalIndex, nSubentityFace, *gridbase_);
    }

    typename Traits::LeafIntersectionIterator ileafbegin () const
    {
    	return ilevelbegin();
    }

    typename Traits::LeafIntersectionIterator ileafend () const
    {
    	return ilevelend();
    }


    /**\brief Inter-level access to father entity on the next-coarser grid.
       The given entity resulted directly from a subdivision of its father
       entity. For the macro elements dereferencing the EntityPointer is undefined.

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    Entity father () const
    {
    	DUNE_THROW(NotImplemented, "CurvilinearGrid-Element: method father() not implemented, since there is no refinement");
    	return Entity(This());
    }

    /**\brief Return true if entity has a father entity which can be accessed
       using the father() method.
     */
    bool hasFather () const  { return false; }

    //! Returns true if the entity is contained in the leaf grid
    // NOTE: All elements are leafs since there is no refinement
    bool isLeaf () const  { return true; }

    /** @brief Returns true if element is of regular type in red/green type refinement.
       In bisection or hanging node refinement this is always true.
     */
    bool isRegular() const {
    	DUNE_THROW(NotImplemented, "CurvilinearGrid-Element: method isRegular() not implemented, since there is no refinement");
    	return false;
    }

    /** \brief Provides information how this element has been subdivided from its father element. */
    LocalGeometry geometryInFather () const {
    	DUNE_THROW(NotImplemented, "CurvilinearGrid-Element: method geometryInFather() not implemented, since there is no refinement");
    	return Geometry(GeometryImpl(gridbase_->template entityGeometry<0>(*gridbaseIndexIterator_), *gridbase_));
    }

    /**\brief Inter-level access to elements that resulted from (recursive) subdivision of this element.
     *
     * Since no refinement implemented, reuse LevelIterator iterating only over the base element itself
     *
     */
    HierarchicIterator hbegin (int maxLevel) const
    {
    	return HierarchicIterator(HierarchicIteratorImpl(
    			gridbase_->entityDuneIndexBegin(0, PartitionIteratorType::All_Partition), *gridbase_)
    	);
    }

    /** \brief Returns iterator to one past the last son element */
    HierarchicIterator hend (int maxLevel) const
    {
    	return ++hbegin(maxLevel);
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return true; }

    /**\brief Returns true, if entity might disappear during the next call to adapt().
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    bool mightVanish () const { return false; }

    /**\brief Returns true, if entity has intersections with the domain (or periodic) boundary
     */
    bool hasBoundaryIntersections () const
    {
    	for (InternalIndexType i = 0; i < 4; i++)
    	{
    		LocalIndexType thisFaceIndex = gridbase_->subentityLocalIndex(*gridbaseIndexIterator_, ELEMENT_CODIM, FACE_CODIM, i);
    		StructuralType thisBoundaryType = gridbase_->faceBoundaryType(thisFaceIndex);

    		if (thisBoundaryType == GridStorageType::FaceBoundaryType::DomainBoundary)  { return true; }
    	}

    	return false;
    }


  };





  } // Namespace CurvGrid

} // Namespace Dune

#endif // DUNE_CURVILINEARGRID_ENTITY_HH
