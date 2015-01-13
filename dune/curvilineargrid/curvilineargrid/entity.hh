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






  /********
   *
   *
   * [TODO] Move all the typedefs to the traits
   *
   * [TODO] Must condition everywhere to return ghost properties when this entity is ghost
   * [TODO] Implement entity seed
   *
   * [TODO] It is currently not allowed to ask GhostElement for its subentityindex. Is that necessary?
   *
   */






namespace Dune
{

  namespace CurvGrid
  {



  template<int codim, int dim, class GridImp>
  class EntityBase
  {
  public:
	  typedef typename remove_const< GridImp >::type::Traits Traits;

	  typedef typename Traits::ctype ctype;

  protected:
	  typedef typename Traits::template Codim< codim >::EntitySeed EntitySeed;

	  typedef typename Traits::GridStorageType      GridStorageType;
	  typedef typename Traits::GridBaseType         GridBaseType;

	  typedef typename Traits::StructuralType       StructuralType;

	  typedef typename Traits::IndexSetIterator     IndexSetIterator;

  public:
		/** \name Construction, Initialization and Destruction
		*  \{ */

	  EntityBase (
	    IndexSetIterator & iter,
	    GridBaseType & gridbase
	    )
	  	  :
	  		gridbaseIndexIterator_(iter),
	  	    gridbase_(gridbase)
	  {  }


	  //! Copy constructor from an existing entity.
	  EntityBase(const EntityBase& other) { }

	  /** \} */

	  //! Copy assignment operator from an existing entity.
	  EntityBase& operator=(const EntityBase& other)  { }

	  //! Move assignment operator from an existing entity.
	  EntityBase& operator=(EntityBase&& other)
	  {
		  realEntity = std::move(other.realEntity);
		  return *this;
	  }

	  //! Move constructor from an existing entity.
	  EntityBase(EntityBase&& other) : realEntity(std::move(other.realEntity))
	  {}


	  /** \brief compare two entities */
	  bool equals ( const EntityBase &other) const  { }

	  /** \brief Return the entity seed which contains sufficient information to generate the entity again and uses as little memory as possible */
	  EntitySeed seed () const  { return EntitySeed(*gridbaseIndexIterator_, pitype_, gridbase_); }


	  /** \brief moves to the next entity within the base storage. Additional functionality used by iterators */
	  void next()  { gridbaseIndexIterator_++; }

  protected:
	    IndexSetIterator gridbaseIndexIterator_;
	    GridBaseType & gridbase_;
  };







  template<int codim, int dim, class GridImp>
  class Entity : EntityBase<codim, dim, GridImp>
  {
	    typedef typename remove_const< GridImp >::type::Traits Traits;

  public:
	    /** \name Attributes
	     *  \{ */

	    static const int codimension = codim;				//! codimensioon of the entity
	    static const int dimension = Traits::dimension;			//! dimension of the grid
	    static const int mydimension = dimension - codimension;		//! dimension of the entity
	    static const int dimensionworld = Traits::dimensionworld;		//! dimension of the world

	    /** \} */

	    /** \name Types Required by DUNE
	     *  \{ */
	    typedef typename Traits::ctype ctype;						//! coordinate type of the grid
	    typedef typename Traits::template Codim< codimension >::Geometry Geometry;	//! type of corresponding geometry
	    /** \} */

	    typedef typename Traits::template Codim< codimension >::GeometryImpl GeometryImpl;

	    typedef typename Traits::GridBaseStorage      GridBaseStorage;
	    typedef typename Traits::GridBaseType         GridBaseType;
	    typedef typename Traits::IdType               IdType;


	    typedef EntityBase<codim, dim, GridImp>   Base;

	    using Base::gridbaseIndexIterator_;
	    using Base::gridbase_;

  public:

	    Entity (int localEntityIndex, GridBaseType & gridbase)
	  	  : Base(localEntityIndex, gridbase)
	    {}

	    /** \name Methods Shared by Entities of All Codimensions
	    *  \{ */

	    /** \brief Return the name of the reference element. The type can be used to access the Dune::ReferenceElement. */
	    GeometryType type () const { return gridbase_.entityGeometryType(codim, *gridbaseIndexIterator_); }

	    /** \brief  Returns the (refinement) level of this entity */
	    int level () const { return gridbase_.entityLevel(codim, *gridbaseIndexIterator_); }

	    /** \brief obtain the partition type of this entity */
	    PartitionType partitionType () const  {
	    	StructuralType stype = gridbase_.entityStructuralType(codim, *gridbaseIndexIterator_);

	    	switch (stype)
	    	{
	    	case GridBaseStorage::PartitionType::Internal          : return PartitionType::InteriorEntity;
	    	case GridBaseStorage::PartitionType::DomainBoundary    : return PartitionType::InteriorEntity;
	    	case GridBaseStorage::PartitionType::ProcessBoundary   : return PartitionType::BorderEntity;
	    	case GridBaseStorage::PartitionType::ComplexBoundary   : return PartitionType::BorderEntity;
	    	case GridBaseStorage::PartitionType::Ghost             : return PartitionType::GhostEntity;
	    	default : return 0;
	    	}
	    }


	    /** \brief obtain geometric realization of the entity */
	    Geometry geometry () const { return gridbase_.entityGeometry<codim>(*gridbaseIndexIterator_); }

	    /** \} */

	    /** \brief Additional method guaranteed to return local entity index as specified in CurvilinearGridBase */
	    int localIndex () const  { return *gridbaseIndexIterator_; }

	    /** \brief obtain the entity's index */
	    int index () const  { return *gridbaseIndexIterator_; }

	    /** \brief obtain the index of a subentity from a host IndexSet
	     *
	     *  \param[in]  i         number of the subentity
	     *  \param[in]  codim        codimension of the subentity
	     */
	    int subIndex (int internalIndex, unsigned int subcodim ) const  {
	    	return gridbase_.subentityIndex(*gridbaseIndexIterator_, codimension, subcodim, internalIndex);
	    }


	    /** \brief obtain the entity's id from a host IdSet */
	    IdType id () const  { return gridbase_.globalId(codim, *gridbaseIndexIterator_); }


	    IdType subId ( int internalIndex, unsigned int subcodim ) const
	    {
	    	int subentityLocalIndex = subIndex(internalIndex, subcodim);
	    	return gridbase_.globalId(subcodim, subentityLocalIndex);
	    }
  };






  template<int dim, class GridImp>
  class Entity <0,dim,GridImp>
  {
	  typedef typename remove_const< GridImp >::type::Traits Traits;

  public:

	  /** \brief The geometry type of this entity */
	  typedef typename Traits::template Codim< 0 >::Geometry Geometry;

	  //! \brief The corresponding entity seed (for storage of entities)
	  typedef typename Traits::template Codim< 0 >::EntitySeed EntitySeed;

	  /** \brief The geometry type of this entity when the geometry is expressed embedded in the father element. */
	  typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

	  /** \brief The HierarchicIterator type*/
	  typedef typename Traits::HierarchicIterator HierarchicIterator;

	  /** \name Attributes
	   *  \{ */
	  static const int codimension    = 0;		//! codimensioon of the entity
	  static const int dimension      = dim;	//! dimension of the grid
	  static const int mydimension    = dim;	//! dimension of the entity
	  static const int dimensionworld = dim;	//! dimension of the world
	  /** \} */

	  typedef typename Traits::GridStorageType      GridStorageType;
	  typedef typename Traits::GridBaseType         GridBaseType;

	  typedef typename GridBaseType::GlobalIndexType           GlobalIndexType;
	  typedef typename GridBaseType::LocalIndexType            LocalIndexType;
	  typedef typename GridBaseType::InternalIndexType         InternalIndexType;
	  typedef typename GridBaseType::StructuralType            StructuralType;
	  typedef typename GridBaseType::PhysicalTagType           PhysicalTagType;
	  typedef typename GridBaseType::InterpolatoryOrderType    InterpolatoryOrderType;


	  typedef EntityBase<0, dim, GridImp>                  Base;
	  using Base::gridbaseIndexIterator_;
	  using Base::gridbase_;

  public:


      Entity (int localEntityIndex, GridBaseType & gridbase)
  	  	  : Base(localEntityIndex, gridbase)
      {}


   /**\brief Number of subentities with codimension <tt>codim</tt>.
     *
     * Strictly speaking this method is redundant, because the same information can be obtained
     * from the corresponding reference element. It is here for efficiency reasons only.
     */
    unsigned int subEntities(unsigned int subcodim) const
    {
      return Dune::ReferenceElements<double, dim>::general().size(0, subcodim, subcodim);
    }

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
    typename Codim< subcodim >::Entity
    subEntity ( int i ) const
    {
    	int subentityLocalIndex = gridbase_.subentityIndex(*gridbaseIndexIterator_, 0, subcodim, i);
    	return Entity<dim - subcodim, dim, GridImp>(subentityLocalIndex, gridbase_);
    }


    /**\brief Inter-level access to father entity on the next-coarser grid.
       The given entity resulted directly from a subdivision of its father
       entity. For the macro elements dereferencing the EntityPointer is undefined.

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    Entity<0, dim, GridImp> father () const
    {
    	DUNE_THROW(NotImplemented, "CurvilinearGrid-Element: method father() not implemented, since there is no refinement");
    	return *this;
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

    /** \brief Provides information how this element has been subdivided from its
     *         father element.
     *
     *  The returned LocalGeometry is a model of
     *  Dune::Geometry<dimension,dimension,...>, mapping the reference element of
     *  the given entity to the reference element of its father.
     *
     *  This information is sufficient to interpolate all degrees of freedom in
     *  the conforming case.
     *  Nonconforming may require access to neighbors of the father and
     *  calculations with local coordinates.
     *  The on-the-fly case is somewhat inefficient since degrees of freedom may be
     *  visited several times.
     *  If we store interpolation matrices, this is tolerable.
     *  We assume that on-the-fly implementation of interpolation is only done for
     *  simple discretizations.
     *
     *  \note For ghost entities, this method is not guaranteed to be implemented.
     *
     *  \note Previously, the geometry was encapsulated in the entity object and
     *        a const reference was returned.
     *
     *  \note The returned geometry object is guaranteed to remain valid until the
     *        grid is modified (or deleted).
     */
    LocalGeometry geometryInFather () const {
    	DUNE_THROW(NotImplemented, "CurvilinearGrid-Element: method geometryInFather() not implemented, since there is no refinement");
    	return gridbase_.entityGeometry<0>(*gridbaseIndexIterator_);
    }

    /**\brief Inter-level access to elements that resulted from (recursive)
       subdivision of this element.

       \param[in] maxlevel Iterator does not stop at elements with level greater than maxlevel.
       \return Iterator to the first son (level is not greater than maxlevel)

       \note If the partitionType of the Entity is GhostEntity,
           it is not guaranteed that this method is working
           or implemented in general.
           For some grids it might be available, though.
     */
    HierarchicIterator hbegin (int maxLevel) const
    {
    	//[FIXME] Add iterator link when implemented
    }

    /** \brief Returns iterator to one past the last son element

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    HierarchicIterator hend (int maxLevel) const
    {
    	//[FIXME] Add iterator link when implemented
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return true; }

    /**\brief Returns true, if entity might disappear during the next call to adapt().
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    bool mightVanish () const { return false; }

    /**\brief Returns true, if entity has intersections with boundary
     */
    bool hasBoundaryIntersections () const
    {
    	for (InternalIndexType i = 0; i < 4; i++)
    	{
    		LocalIndexType thisFaceIndex = gridbase_.subentityIndex(*gridbaseIndexIterator_, 0, 1, i);
    		StructuralType thisStructType = entityStructuralType<1>(thisFaceIndex);

    		if (
    		  (thisStructType == GridStorageType::PartitionType::ProcessBoundary) ||
    		  (thisStructType == GridStorageType::PartitionType::ComplexBoundary)
    		)  { return true; }
    	}

    	return false;
    }


  };





  } // Namespace CurvGrid

} // Namespace Dune

#endif // DUNE_CURVILINEARGRID_ENTITY_HH
