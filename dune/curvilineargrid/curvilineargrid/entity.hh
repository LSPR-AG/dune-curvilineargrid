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



  template<int codim, int dim, class GridImp>
  class CurvEntityBase
  {
  public:
	  typedef typename remove_const< GridImp >::type::Traits Traits;

	  typedef typename Traits::ctype ctype;

  protected:
	  typedef typename Traits::template Codim< codim >::EntitySeed EntitySeed;

	  typedef Dune::CurvilinearGridBase<ctype,dim>          GridBaseType;
	  typedef typename GridBaseType::IndexSetIterator       IndexSetIterator;

  public:
		/** \name Construction, Initialization and Destruction
		*  \{ */

	  CurvEntityBase (
	    IndexSetIterator & iter,
	    GridBaseType & gridbase,
	    Dune::PartitionType pitype
	    )
	  	  :
	  		gridbaseIndexIterator_(iter),
	  	    gridbase_(gridbase),
	  	    pitype_(pitype)
	  {  }


	  //! Copy constructor from an existing entity.
	  CurvEntityBase(const CurvEntityBase& other)
	  {
		  gridbaseIndexIterator_ = other.gridbaseIndexIterator_;
		  gridbase_ = other.gridbase_;
		  pitype_ = other.pitype_;
	  }

	  //! Move constructor from an existing entity.
	  //CurvEntityBase(CurvEntityBase&& other) : realEntity(std::move(other.realEntity)) {}


	  //! Copy assignment operator from an existing entity.
	  //CurvEntityBase& operator=(const CurvEntityBase& other)  { }

	  //! Move assignment operator from an existing entity.
	  //CurvEntityBase& operator=(CurvEntityBase&& other)  { realEntity = std::move(other.realEntity);  return *this; }



	  /** \} */


	  /** \brief compare two entities */
	  bool equals ( const CurvEntityBase &other) const
	  {
		  return ((gridbaseIndexIterator_ == other.gridbaseIndexIterator_) && (gridbase_ == other.gridbase_));
	  }


	  /** \brief Return the entity seed which contains sufficient information to generate the entity again and uses as little memory as possible */
	  EntitySeed seed () const  { return EntitySeed(*gridbaseIndexIterator_, pitype_, gridbase_); }


	  /** \brief moves to the next entity within the base storage. Additional functionality used by iterators */
	  void next()  { gridbaseIndexIterator_++; }

  protected:
	    IndexSetIterator gridbaseIndexIterator_;
	    GridBaseType & gridbase_;
	    Dune::PartitionType pitype_;
  };







  template<int codim, int dim, class GridImp>
  class CurvEntity : CurvEntityBase<codim, dim, GridImp>
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

	    typedef Dune::CurvilinearGridStorage<ctype,dim>           GridBaseStorage;
	    typedef Dune::CurvilinearGridBase<ctype,dim>              GridBaseType;
	    typedef typename GridBaseStorage::IdType                  IdType;


	    typedef CurvEntityBase<codim, dim, GridImp>   Base;

	    using Base::gridbaseIndexIterator_;
	    using Base::gridbase_;

  public:

	    CurvEntity (int localEntityIndex, GridBaseType & gridbase, Dune::PartitionType pitype)
	  	  : Base(localEntityIndex, gridbase, pitype)
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
	    	case GridBaseStorage::PartitionType::Ghost             : return PartitionType::GhostEntity;
	    	default : return 0;
	    	}
	    }


	    /** \brief obtain geometric realization of the entity */
	    Geometry geometry () const { return Geometry(GeometryImpl(gridbase_.entityGeometry<codim>(*gridbaseIndexIterator_))); }

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
  class CurvEntity <0, dim, GridImp>
  {
	  typedef typename remove_const< GridImp >::type::Traits Traits;
	  typedef typename Traits::ctype ctype;						//! coordinate type of the grid

  public:

	  /** \name Attributes
	   *  \{ */
	  static const int codimension    = 0;		//! codimensioon of the entity
	  static const int dimension      = dim;	//! dimension of the grid
	  static const int mydimension    = dim;	//! dimension of the entity
	  static const int dimensionworld = dim;	//! dimension of the world
	  /** \} */


	  /** \brief The geometry type of this entity */
	  typedef typename Traits::template Codim< 0 >::Geometry Geometry;

	  typedef typename Traits::template Codim< codimension >::GeometryImpl GeometryImpl;

	  //! \brief The corresponding entity seed (for storage of entities)
	  typedef typename Traits::template Codim< 0 >::EntitySeed EntitySeed;

	  /** \brief The geometry type of this entity when the geometry is expressed embedded in the father element. */
	  typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

	  /** \brief The HierarchicIterator type*/
	  typedef typename Traits::HierarchicIterator               HierarchicIterator;
	  typedef Dune::CurvGrid::CurvHierarchicIterator<GridImp>   HierarchicIteratorImpl;

	  typedef Dune::CurvilinearGridStorage<ctype,dim>           GridBaseStorage;
	  typedef Dune::CurvilinearGridBase<ctype,dim>              GridBaseType;

	  typedef typename GridBaseType::GlobalIndexType           GlobalIndexType;
	  typedef typename GridBaseType::LocalIndexType            LocalIndexType;
	  typedef typename GridBaseType::InternalIndexType         InternalIndexType;
	  typedef typename GridBaseType::StructuralType            StructuralType;
	  typedef typename GridBaseType::PhysicalTagType           PhysicalTagType;
	  typedef typename GridBaseType::InterpolatoryOrderType    InterpolatoryOrderType;


	  typedef CurvEntityBase<0, dim, GridImp>                  Base;
	  using Base::gridbaseIndexIterator_;
	  using Base::pitype_;
	  using Base::gridbase_;

  public:


      CurvEntity (int localEntityIndex, GridBaseType & gridbase, Dune::PartitionType pitype)
  	  	  : Base(localEntityIndex, gridbase, pitype)
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
    typename Traits::Codim< subcodim >::Entity
    subEntity ( int i ) const
    {
    	int subentityLocalIndex = gridbase_.subentityIndex(*gridbaseIndexIterator_, 0, subcodim, i);
    	return Traits::template Codim< subcodim >::Entity(CurvEntity<dim - subcodim, dim, GridImp>(subentityLocalIndex, gridbase_, pitype_));
    }


    /**\brief Inter-level access to father entity on the next-coarser grid.
       The given entity resulted directly from a subdivision of its father
       entity. For the macro elements dereferencing the EntityPointer is undefined.

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    CurvEntity<0, dim, GridImp> father () const
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

    /** \brief Provides information how this element has been subdivided from its father element. */
    LocalGeometry geometryInFather () const {
    	DUNE_THROW(NotImplemented, "CurvilinearGrid-Element: method geometryInFather() not implemented, since there is no refinement");
    	return Geometry(GeometryImpl(gridbase_.entityGeometry<0>(*gridbaseIndexIterator_)));
    }

    /**\brief Inter-level access to elements that resulted from (recursive) subdivision of this element.
     *
     * Since no refinement implemented, reuse LevelIterator iterating only over the base element itself
     *
     */
    HierarchicIterator hbegin (int maxLevel) const
    {

    	return HierarchicIterator(HierarchicIteratorImpl(
    			gridbase_.entityDuneIndexBegin(0, PartitionIteratorType::All_Partition), gridbase_)
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

    /**\brief Returns true, if entity has intersections with boundary
     */
    bool hasBoundaryIntersections () const
    {
    	for (InternalIndexType i = 0; i < 4; i++)
    	{
    		LocalIndexType thisFaceIndex = gridbase_.subentityIndex(*gridbaseIndexIterator_, 0, 1, i);
    		StructuralType thisStructType = entityStructuralType<1>(thisFaceIndex);

    		if (thisStructType == GridBaseStorage::PartitionType::ProcessBoundary)  { return true; }
    	}

    	return false;
    }


  };





  } // Namespace CurvGrid

} // Namespace Dune

#endif // DUNE_CURVILINEARGRID_ENTITY_HH
