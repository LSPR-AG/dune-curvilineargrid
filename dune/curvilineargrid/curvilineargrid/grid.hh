// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GRID_HH
#define DUNE_CURVGRID_GRID_HH

#include <dune/common/nullptr.hh>

#include <dune/grid/common/grid.hh>

#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>

namespace Dune
{

	// Forwards-Declaration
	// ****************************************************************************************
	template <int dim, int dimworld, class ct>                          class CurvilinearGrid;
	template< int mydim, int cdim, class GridImp >                      class CurvGeometry;
	template< int codim, int dim, class GridImp>                        class CurvEntity;
	template< int codim, class GridImp >                                class CurvEntityPointer;
	template< int codim, class GridImp >                                class CurvEntitySeed;
	template< int codim, PartitionIteratorType pitype, class GridImp >  class CurvLevelIterator;
	template< class GridImp >                                           class CurvIntersectionIterator;
	template< class GridImp >                                           class CurvIntersection;
	template< class GridImp >                                           class CurvHierarchicIterator;
	template< class GridImp >                                           class CurvIndexSet;
	template< class GridImp >                                           class CurvIdSet;


	// GridFamily
	// ****************************************************************************************
	template<int dim, int dimworld, class ct>
	struct CurvGridFamily
	{
#if HAVE_MPI
	    typedef CollectiveCommunication<MPI_Comm> CCType;
#else
	    typedef CollectiveCommunication<No_Comm> CCType;
#endif

	    typedef Dune::CurvilinearGridStorage<ct, dimworld>::IdType CurvIdType;

	    typedef GridTraits<dim,                                     // dimension of the grid
	        dimworld,                                               // dimension of the world space
	        Dune::CurvilinearGrid<dim, dimworld, ct>,
	        CurvGeometry,
	        CurvEntity,
	        CurvEntityPointer,
	        CurvLevelIterator,                                      // type used for the level iterator
	        CurvIntersection,              // leaf  intersection
	        CurvIntersection,              // level intersection
	        CurvIntersectionIterator,              // leaf  intersection iter
	        CurvIntersectionIterator,              // level intersection iter
	        CurvHierarchicIterator,
	        CurvLevelIterator,                                      // type used for the leaf(!) iterator
	        CurvIndexSet,                  // level index set
	        CurvIndexSet,                  // leaf index set
	        CurvIdSet,
	        CurvIdType,
	        CurvIdSet,
	        CurvIdType,
	        CCType,
	        DefaultLevelGridViewTraits,
	        DefaultLeafGridViewTraits,
	        CurvEntitySeed>
	    Traits;
	  };


















  template <int dim, int dimworld, class ct>
  class CurvilinearGrid: public GridDefaultImplementation < dim, dimworld, ct, CurvGrid::GridFamily< dim, dimworld, ct > >
      /** \endcond */
  {
    typedef CurvilinearGrid<dim, dimworld, ct> Grid;
    typedef GridDefaultImplementation < dim, dimworld, ct, CurvGrid::GridFamily< dim, dimworld, ct > > Base;

  public:

    typedef CurvGrid::GridFamily< dim, dimworld, ct > GridFamily;
    typedef typename GridFamily::Traits Traits;                               //! type of the grid traits

    // Curvilinear Grid Implementation
    typedef Dune::CurvilinearGridBase<ct, dimworld> CurvGridBase;
    typedef Dune::CurvilinearGridStorage<ct, dimworld> CurvGridStorage;

    static const int   VERTEX_CODIM   = CurvGridStorage::VERTEX_CODIM;
    static const int   EDGE_CODIM     = CurvGridStorage::EDGE_CODIM;
    static const int   FACE_CODIM     = CurvGridStorage::FACE_CODIM;
    static const int   ELEMENT_CODIM  = CurvGridStorage::ELEMENT_CODIM;







    /** \brief traits structure containing types for a codimension
     *
     *  \tparam codim  codimension
     *
     *  \nosubgrouping
     */
    template< int codim >
    struct Codim;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    //! iterator over the grid hierarchy
    typedef typename Traits::HierarchicIterator             HierarchicIterator;
    //! iterator over intersections with other entities on the leaf level
    typedef typename Traits::LeafIntersectionIterator       LeafIntersectionIterator;
    //! iterator over intersections with other entities on the same level
    typedef typename Traits::LevelIntersectionIterator      LevelIntersectionIterator;

    /** \} */

    /** \name Grid View Types
     *  \{ */

    /** \brief Types for GridView */
    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename GridFamily::Traits::template Partition< pitype >::LevelGridView    LevelGridView;
      typedef typename GridFamily::Traits::template Partition< pitype >::LeafGridView     LeafGridView;
    };

    /** \brief View types for All_Partition */
    typedef typename Partition< All_Partition >::LevelGridView     LevelGridView;
    typedef typename Partition< All_Partition >::LeafGridView      LeafGridView;

    /** \} */


    /** \name Index and Id Set Types
     *  \{ */

    typedef typename Traits::LeafIndexSet   LeafIndexSet;
    typedef typename Traits::LevelIndexSet  LevelIndexSet;
    typedef typename Traits::GlobalIdSet    GlobalIdSet;
    typedef typename Traits::LocalIdSet     LocalIdSet;

    /** \} */

    /** \name Miscellaneous Types
     * \{ */

    //! type of vector coordinates (e.g., double)
    typedef typename Traits::ctype ctype;

    //! communicator with all other processes having some part of the grid
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    /** \} */

    /** \name Construction and Destruction
     *  \{ */

    /** \brief constructor
     *
     *  [FIXME] Must initialize levelIndexSets_
     */
    CurvilinearGrid (GridBaseType & gridbase, MPIHelper &mpihelper)
      : gridbase_(gridbase),
        mpihelper_(mpihelper)
    {}


    /** \brief destructor */
    ~CurvilinearGrid ()  { }

    /** \} */

    /** \name Size Methods
     *  \{ */

    /** \brief obtain maximal grid level
     *
     *  Grid levels are numbered 0, ..., L, where L is the value returned by
     *  this method.
     *
     *  \returns maximal grid level
     */
    int maxLevel () const  { return 0; }


    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of entities of codimension \em codim on grid level
     *           \em level.
     */
    int size ( int level, int codim ) const
    {
      return (level == 0) ? size(codim) : 0;
    }


    /** \brief obtain number of leaf entities
     *
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of leaf entities of codimension \em codim
     */
    int size ( int codim ) const
    {
      return gridbase_.nEntity(codim);
    }

    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of entities with a geometry of type \em type on grid
     *           level \em level.
     */
    int size ( int level, GeometryType type ) const
    {
    	return (level == 0) ? size (type) : 0;
    }

    /** \brief obtain number of leaf entities
     *
     *  \returns number of leaf entities with a geometry of type \em type
     */
    int size ( GeometryType type ) const
    {
    	return type.isSimplex() ?  size ( dim - type.dim() ) : 0;
    }

    /** \brief returns the number of boundary segments within the macro grid
     *
     *  \returns number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const
    {
    	return gridbase_.nEntity(FACE_CODIM, DomainBoundaryType);
    }
    /** \} */

    const GlobalIdSet &globalIdSet () const
    {
        if( !globalIdSet_ )  { globalIdSet_ = GlobalIdSet(); }
        assert( globalIdSet_ );
        return globalIdSet_;
    }

    const LocalIdSet &localIdSet () const
    {
        return globalIdSet ();
    }

    const LevelIndexSet &levelIndexSet ( int level ) const
    {
        if( (level < 0) || (level > maxLevel()) )
        {
          DUNE_THROW( GridError, "LevelIndexSet for nonexisting level " << level << " requested." );
        }

        return leafIndexSet();
    }

    const LeafIndexSet &leafIndexSet () const
    {
      if( !leafIndexSet_ )  { leafIndexSet_ = LeafIndexSet(gridbase_); }
      assert( leafIndexSet_ );
      return leafIndexSet_;
    }

    // [TODO] Not Implemented yet
    void globalRefine ( int refCount )
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: globalRefine() called but refinement has not been implemented" << std::endl;
    }

    bool mark ( int refCount, const typename Codim< 0 >::Entity &entity )
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: mark() called but refinement has not been implemented" << std::endl;
    	return false;
    }

    int getMark ( const typename Codim< 0 >::Entity &entity ) const
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: getMark() called but refinement has not been implemented" << std::endl;
    	return 0;
    }

    bool preAdapt ()
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: preAdapt() called but refinement has not been implemented" << std::endl;
    	return false;
    }

    bool adapt ()
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: adapt() called but refinement has not been implemented" << std::endl;
    	return false;
    }

    void postAdapt ()
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: postAdapt() called but refinement has not been implemented" << std::endl;
    }

    /** \name Parallel Data Distribution and Communication Methods
     *  \{ */

    /** \brief obtain size of overlap region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int overlapSize ( int codim ) const  { return 0; }

    /** \brief obtain size of ghost region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int ghostSize( int codim ) const  { return gridbase_.nEntity(codim, Dune::CurvilinearGridStorage<ct, dimworld>::EntityStructuralType::GhostElement); }

    /** \brief obtain size of overlap region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     *
     *  \note There is only one level at the moment
     */
    int overlapSize ( int level, int codim ) const  { return (level == 0) ? overlapSize (codim) : 0; }

    /** \brief obtain size of ghost region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     *
     *  \note There is only one level at the moment
     */
    int ghostSize ( int level, int codim ) const  { return (level == 0) ? ghostSize(codim) : 0;  }

    /** \brief communicate information on a grid level
     *
     *  \param      dataHandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication or BackwardCommunication)
     *  \param[in]  level       grid level to communicate
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction,
                       int level ) const
    {
    	if (level == 0)  { communicate (dataHandle, interface,  direction); }
    }

    /** \brief communicate information on leaf entities
     *
     *  \param      dataHandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication, BackwardCommunication)
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
    	Dune::CurvGrid::Communication<Grid> communicator;
    	if (dataHandle.contains(dimension, ELEMENT_CODIM)) { communicator.communicateWrapper<DataHandle, ELEMENT_CODIM>(dataHandle, interface, direction); }
    	if (dataHandle.contains(dimension, FACE_CODIM))    { communicator.communicateWrapper<DataHandle, FACE_CODIM>(dataHandle, interface, direction); }
    	if (dataHandle.contains(dimension, EDGE_CODIM))    { communicator.communicateWrapper<DataHandle, EDGE_CODIM>(dataHandle, interface, direction); }
    	if (dataHandle.contains(dimension, VERTEX_CODIM))  { communicator.communicateWrapper<DataHandle, VERTEX_CODIM>(dataHandle, interface, direction); }
    }

    /** \brief obtain CollectiveCommunication object
     *
     *  The CollectiveCommunication object should be used to globally
     *  communicate information between all processes sharing this grid.
     *
     *  \note The CollectiveCommunication object returned is identical to the
     *        one returned by the host grid.
     */
    const CollectiveCommunication &comm () const  { return mpihelper_.getCollectiveCommunication(); }


    // data handle interface different between geo and interface

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \returns \b true, if the grid has changed.
     */
    bool loadBalance ()
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: loadBalance() called but refinement has not been implemented. CurvGrid gets balanced only at initialization" << std::endl;
    	return false;
    }

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  The data handle is used to communicate the data associated with
     *  entities that move from one process to another.
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \param  datahandle  communication data handle (user defined)
     *
     *  \returns \b true, if the grid has changed.
     */

    template< class DataHandle, class Data >
    bool loadBalance ( CommDataHandleIF< DataHandle, Data > &datahandle )
    {
    	//DUNE_THROW( NotImplemented, "Refinement not implemented in CurvGrid at the moment" );
    	std::cout << "::: loadBalance() called but refinement has not been implemented. CurvGrid gets balanced only at initialization" << std::endl;
    	return false;
    }


    /** \brief obtain EntityPointer from EntitySeed. */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::EntityPointer
    entityPointer ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( seed );
    }

    /** \} */


    const CurvGridBase & curvGridBase()  { return gridbase_; }



  private:
    mutable CurvGridBase & gridbase_;
    mutable MPIHelper & mpihelper_;

    //mutable std::vector< LevelIndexSet *, typename Allocator::template rebind< LevelIndexSet * >::other > levelIndexSets_;
    mutable LeafIndexSet leafIndexSet_;
    mutable GlobalIdSet globalIdSet_;
  };



  // CurvilinearGrid::Codim
  // -------------------

  template<int dim, int dimworld, class ct>
  template< int codim >
  struct CurvilinearGrid< dim, dimworld, ct>::Codim
    : public Base::template Codim< codim >
  {
    /** \name Entity and Entity Pointer Types
     *  \{ */

    /** \brief type of entity
     *
     *  The entity is a model of Dune::Entity.
     */
    typedef typename Traits::template Codim< codim >::Entity Entity;

    /** \brief type of entity pointer
     *
     *  The entity pointer is a model of Dune::EntityPointer.
     */
    typedef typename Traits::template Codim< codim >::EntityPointer EntityPointer;

    /** \} */

    /** \name Geometry Types
     *  \{ */

    /** \brief type of world geometry
     *
     *  Models the geomtry mapping of the entity, i.e., the mapping from the
     *  reference element into world coordinates.
     *
     *  The geometry is a model of Dune::Geometry, implemented through the
     *  generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::Geometry Geometry;

    /** \brief type of local geometry
     *
     *  Models the geomtry mapping into the reference element of dimension
     *  \em dimension.
     *
     *  The local geometry is a model of Dune::Geometry, implemented through
     *  the generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::LocalGeometry LocalGeometry;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Traits::template Codim< codim > ::template Partition< pitype >::LeafIterator   LeafIterator;
      typedef typename Traits::template Codim< codim > ::template Partition< pitype >::LevelIterator  LevelIterator;
    };

    /** \brief type of level iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of a
     *  grid level.
     *
     *  The level iterator is a model of Dune::LevelIterator.
     */
    typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

    /** \brief type of leaf iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of the
     *  leaf grid.
     *
     *  The leaf iterator is a model of Dune::LeafIterator.
     */
    typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

    /** \} */
  };

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRID_HH
