// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVILINEARGRID_ENTITY_HH
#define DUNE_CURVILINEARGRID_ENTITY_HH

#include <dune/common/iteratorrange.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/dimension.hh>

#include "grid.hh"
#include "entitypointer.hh"
#include "rangegenerators.hh"







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

  template<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>
  class Entity
  {
	    typedef typename remove_const< GridImp >::type::Traits Traits;

  public:
	    /** \name Attributes
	     *  \{ */

	    static const int codimension = cd;				//! codimensioon of the entity
	    static const int dimension = Traits::dimension;			//! dimension of the grid
	    static const int mydimension = dimension - codimension;		//! dimension of the entity
	    static const int dimensionworld = Traits::dimensionworld;		//! dimension of the world

	    /** \} */

	    /** \name Types Required by DUNE
	     *  \{ */
	    typedef typename Traits::ctype ctype;						//! coordinate type of the grid
	    typedef typename Traits::template Codim< codimension >::Geometry Geometry;	//! type of corresponding geometry
	    /** \} */

	    typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;			//! type of corresponding entity seed
	    typedef typename Traits::template Codim< codim >::GeometryImpl GeometryImpl;


	    typedef Dune::CurvilinearGridStorage::IdType          IdType;

  public:

		/** \name Construction, Initialization and Destruction
		*  \{ */

	    Entity (
	    	int localEntityIndex,
	    	int structType,
	    	Dune::CurvilinearGridBase<ctype> & gridbase)
  	  	  	  : localEntityIndex_(localEntityIndex), structType_(structType), gridbase_(gridbase)
		{ }


		//! Copy constructor from an existing entity.
		Entity(const Entity& other) { }

		/** \} */

		//! Copy assignment operator from an existing entity.
		Entity& operator=(const Entity& other)  { }

		//! Move assignment operator from an existing entity.
		Entity& operator=(Entity&& other)
		{
			realEntity = std::move(other.realEntity);
			return *this;
		}

		//! Move constructor from an existing entity.
		Entity(Entity&& other) : realEntity(std::move(other.realEntity))
		{}


		/** \brief compare two entities */
		bool equals ( const EntityBase &other) const  { }


	    /** \name Methods Shared by Entities of All Codimensions
	    *  \{ */

	    /** \brief Return the name of the reference element. The type can be used to access the Dune::ReferenceElement. */
	    GeometryType type () const { return gridbase_.entityGeometryType<cd>(localEntityIndex_); }

	    /** \brief  Returns the (refinement) level of this entity */
	    int level () const { return gridbase_.entityLevel<cd>(localEntityIndex_); }

	    /** \brief obtain the partition type of this entity */
	    PartitionType partitionType () const  {
	    	switch (structType_)
	    	{
	    	case Dune::CurvilinearGridStorage::EntityStructuralType::DomainBoundary   : return PartitionType::BorderEntity;
	    	case Dune::CurvilinearGridStorage::EntityStructuralType::ProcessBoundary  : return PartitionType::BorderEntity;
	    	case Dune::CurvilinearGridStorage::EntityStructuralType::InternalFace     : return PartitionType::InteriorEntity;
	    	case Dune::CurvilinearGridStorage::EntityStructuralType::InternalElement  : return PartitionType::InteriorEntity;
	    	case Dune::CurvilinearGridStorage::EntityStructuralType::GhostElement     : return PartitionType::GhostEntity;
	    	default : return 0;
	    	}
	    }


	    /** \brief obtain geometric realization of the entity */
	    Geometry geometry () const { return gridbase_.entityGeometry<cd>(localEntityIndex_); }

	    /** \brief Return the entity seed which contains sufficient information to generate the entity again and uses as little memory as possible */
	    EntitySeed seed () const  {  }

	    /** \} */


	    /** \brief obtain the entity's index from a host IndexSet */
	    int index () const  { return localEntityIndex_; }

	    /** \brief obtain the index of a subentity from a host IndexSet
	     *
	     *  \param[in]  i         number of the subentity
	     *  \param[in]  cd        codimension of the subentity
	     */
	    int subIndex (int internalIndex, unsigned int subcodim ) const  {
	    	return gridbase_.subentityIndex(localEntityIndex_, codimension, subcodim, internalIndex);
	    }


	    /** \brief obtain the entity's id from a host IdSet */
	    IdType id () const  { }









  private:
     int localEntityIndex_;
     int structType_;
     Dune::CurvilinearGridBase<ctype> & gridbase_;




  };






  template<int dim, class GridImp, template<int,int,class> class EntityImp>
  class Entity <0,dim,GridImp,EntityImp>
  {

  public:
	  typedef typename remove_const< Grid >::type::Traits Traits;

	  /** \brief The geometry type of this entity */
	  typedef typename Traits::template Codim< 0 >::Geometry Geometry;

	  //! \brief The corresponding entity seed (for storage of entities)
	  typedef typename Traits::template Codim< 0 >::EntitySeed EntitySeed;

	  /** \brief The geometry type of this entity when the geometry is expressed embedded in the father element. */
	  typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

	  /** \brief EntityPointer types of the different codimensions */
	  template <int cd>
	  struct Codim
	  {
		  typedef typename GridImp::template Codim<cd>::Entity Entity;
	  };

	  /** \brief The HierarchicIterator type*/
	  typedef typename Traits::HierarchicIterator HierarchicIterator;

	  /** \name Attributes
	   *  \{ */
	  static const int codimension = 0;		//! codimensioon of the entity
	  static const int dimension = dim;		//! dimension of the grid
	  static const int mydimension = dim;	//! dimension of the entity
	  static const int dimensionworld = dim;	//! dimension of the world
	  /** \} */


  public:


      Entity ()  {}




   /**\brief Number of subentities with codimension <tt>codim</tt>.
     *
     * Strictly speaking this method is redundant, because the same information can be obtained
     * from the corresponding reference element. It is here for efficiency reasons only.
     */
    unsigned int subEntities(unsigned int codim) const
    {
      return realEntity.subEntities(codim);
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
    template< int codim >
    typename Codim< codim >::Entity
    subEntity ( int i ) const
    {
      return realEntity.template subEntity< codim >( i );
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
      return realEntity.father();
    }

    /**\brief Return true if entity has a father entity which can be accessed
       using the father() method.
     */
    bool hasFather () const
    {
      return realEntity.hasFather();
    }

    //! Returns true if the entity is contained in the leaf grid
    bool isLeaf () const
    {
      return realEntity.isLeaf();
    }

    /** @brief Returns true if element is of regular type in red/green type refinement.
       In bisection or hanging node refinement this is always true.
     */
    bool isRegular() const { return realEntity.isRegular(); }

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
    LocalGeometry geometryInFather () const { return realEntity.geometryInFather(); }

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
      return realEntity.hbegin(maxLevel);
    }

    /** \brief Returns iterator to one past the last son element

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    HierarchicIterator hend (int maxLevel) const
    {
      return realEntity.hend(maxLevel);
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return realEntity.isNew(); }

    /**\brief Returns true, if entity might disappear during the next call to adapt().
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    bool mightVanish () const { return realEntity.mightVanish(); }

    /**\brief Returns true, if entity has intersections with boundary
     */
    bool hasBoundaryIntersections () const { return realEntity.hasBoundaryIntersections(); }


  };





  } // Namespace CurvGrid

} // Namespace Dune

#endif // DUNE_CURVILINEARGRID_ENTITY_HH
