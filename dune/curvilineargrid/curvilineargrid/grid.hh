// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GRID_HH
#define DUNE_CURVGRID_GRID_HH

#include <dune/common/nullptr.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>

#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/communication.hh>
#include <dune/curvilineargrid/curvilineargrid/entity.hh>
#include <dune/curvilineargrid/curvilineargrid/entitypointer.hh>
#include <dune/curvilineargrid/curvilineargrid/entityseed.hh>
#include <dune/curvilineargrid/curvilineargrid/geometry.hh>
#include <dune/curvilineargrid/curvilineargrid/idset.hh>
#include <dune/curvilineargrid/curvilineargrid/indexsets.hh>
#include <dune/curvilineargrid/curvilineargrid/intersection.hh>
#include <dune/curvilineargrid/curvilineargrid/intersectioniterator.hh>
#include <dune/curvilineargrid/curvilineargrid/iterator.hh>


namespace Dune
{

    namespace CurvGrid {


	// Forwards-Declaration
	// ****************************************************************************************
	template< class ct,  int cdim, bool isCached, class LogMsg>         class CurvilinearGrid;
	template< int mydim, int cdim, class GridImp >                      class CurvGeometry;
	template< int codim, int dim,  class GridImp>                       class CurvEntity;
	template< int codim, class GridImp >                                class CurvEntityPointer;
	template< int codim, class GridImp >                                class CurvEntitySeed;
	template< int codim, PartitionIteratorType pitype, class GridImp >  class CurvLevelIterator;
	template< class GridImp >                                           class CurvIntersectionIterator;
	template< class GridImp >                                           class CurvIntersection;
	template< class GridImp >                                           class CurvHierarchicIterator;
	template< class GridImp >                                           class CurvIndexSet;
	template< class GridImp >                                           class CurvIdSet;

    }

	// GridFamily
	// ****************************************************************************************
    template <class ct, int cdim, bool isCached, class LogMsg>
	struct CurvGridFamily
	{
#if HAVE_MPI
	    typedef CollectiveCommunication<MPI_Comm> CCType;
#else
	    typedef CollectiveCommunication<No_Comm> CCType;
#endif

	    typedef ct  ctype;

	    typedef Dune::CurvilinearGrid<ct, cdim, isCached, LogMsg>               GridType;
	    typedef typename Dune::CurvilinearGridBase<ct, cdim, isCached, LogMsg>  GridBaseType;
	    typedef typename GridBaseType::GridStorageType                          GridStorageType;

	    typedef typename GridStorageType::LocalIndexType                 LocalIndexType;
	    typedef typename GridStorageType::IdType                         CurvIdType;

	    typedef GridTraits<cdim,                                // dimension of the grid
	        cdim,                                               // dimension of the world space
	        GridType,
	        CurvGrid::CurvGeometry,
	        CurvGrid::CurvEntity,
	        CurvGrid::CurvEntityPointer,
	        CurvGrid::CurvLevelIterator,             // type used for the level iterator
	        CurvGrid::CurvIntersection,              // leaf  intersection
	        CurvGrid::CurvIntersection,              // level intersection
	        CurvGrid::CurvIntersectionIterator,      // leaf  intersection iter
	        CurvGrid::CurvIntersectionIterator,      // level intersection iter
	        CurvGrid::CurvHierarchicIterator,
	        CurvGrid::CurvLevelIterator,             // type used for the leaf(!) iterator
	        CurvGrid::CurvIndexSet<const GridType>,  // level index set
	        CurvGrid::CurvIndexSet<const GridType>,  // leaf index set
	        CurvGrid::CurvIdSet<const GridType>,
	        CurvIdType,
	        CurvGrid::CurvIdSet<const GridType>,
	        CurvIdType,
	        CCType,
	        DefaultLevelGridViewTraits,
	        DefaultLeafGridViewTraits,
	        CurvGrid::CurvEntitySeed>
	    Traits;
	  };





  template <class ct, int cdim, bool isCached, class LogMsg>
  class CurvilinearGrid : public GridDefaultImplementation < cdim, cdim, ct, Dune::CurvGridFamily< ct, cdim, isCached, LogMsg > >
      /** \endcond */
  {
    typedef CurvilinearGrid<ct, cdim, isCached, LogMsg> Grid;
    typedef GridDefaultImplementation < cdim, cdim, ct, Dune::CurvGridFamily< ct, cdim, isCached, LogMsg > > Base;

    template< int, int, class >                    friend class Dune::CurvGrid::CurvEntityBase;
    template< int, class >                         friend class Dune::CurvGrid::CurvEntity;
    template< int, class >                         friend class Dune::CurvGrid::CurvEntityPointer;
    template< int, class >                         friend class Dune::CurvGrid::CurvEntitySeed;
    template< class >                              friend class Dune::CurvGrid::CurvIntersection;
    template< class >                              friend class Dune::CurvGrid::CurvIntersectionIterator;
    template< int, PartitionIteratorType , class>  friend class Dune::CurvGrid::CurvLevelIterator;
    template< class >                              friend class Dune::CurvGrid::CurvHierarchicIterator;
    template< class >                              friend class Dune::CurvGrid::CurvIdSet;
    template< class >                              friend class Dune::CurvGrid::CurvIndexSet;
    template< int, int, class >                    friend class Dune::CurvGrid::CurvGeometry;

  public:

    typedef Dune::CurvGridFamily< ct, cdim, isCached, LogMsg > GridFamily;
    typedef typename GridFamily::Traits Traits;                               //! type of the grid traits

    static const int dimension      = cdim;
    static const int dimensionworld = cdim;
    static const bool is_cached     = isCached;
    typedef ct                        ctype;
    typedef LogMsg                    LoggingMessage;


    // Curvilinear Grid Implementation
    // ************************************************************************************

    typedef Dune::CurvilinearGridBase<ct, cdim, isCached, LoggingMessage>  GridBaseType;
    typedef typename GridBaseType::GridStorageType                         GridStorageType;

    typedef typename GridStorageType::LocalIndexType          LocalIndexType;
    typedef typename GridStorageType::GlobalIndexType         GlobalIndexType;
    typedef typename GridStorageType::PhysicalTagType         PhysicalTagType;
    typedef typename GridStorageType::StructuralType          StructuralType;
    typedef typename GridStorageType::InterpolatoryOrderType  InterpolatoryOrderType;

    // Common Partition types
    static const unsigned int DomainBoundaryType   = GridBaseType::DomainBoundaryType;
    static const unsigned int ProcessBoundaryType  = GridBaseType::ProcessBoundaryType;
    static const unsigned int InternalType         = GridBaseType::InternalType;
    static const unsigned int GhostType            = GridBaseType::GhostType;

    static const unsigned int BOUNDARY_SEGMENT_PARTITION_TYPE = GridStorageType::BOUNDARY_SEGMENT_PARTITION_TYPE;
    static const int NO_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::None;
    static const int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;


    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;




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

    typedef Dune::CurvGrid::CurvIndexSet<const Grid>    IndexSetImpl;
    typedef Dune::CurvGrid::CurvIdSet<const Grid>       IdSetImpl;

    /** \} */

    /** \name Miscellaneous Types
     * \{ */

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
      : gridbase_(&gridbase),
        mpihelper_(mpihelper),
        globalIdSet_(),
        leafIndexSet_(gridbase)
    {

    }


    /** \brief destructor */
    ~CurvilinearGrid ()  { delete gridbase_; }

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
      return gridbase_->nEntity(codim);
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
    	return type.isSimplex() ?  size ( dimension - type.dim() ) : 0;
    }

    /** \brief returns the number of boundary segments within the macro grid
     *
     *  \returns number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const
    {
    	return gridbase_->nEntity(FACE_CODIM, PartitionType::InteriorEntity, DOMAIN_BOUNDARY_TYPE);
    }
    /** \} */

    size_t numProcessBoundaries () const
    {
    	return gridbase_->nEntity(FACE_CODIM, PartitionType::BorderEntity);
    }

    size_t numInternal(int codim) const
    {
    	return gridbase_->nEntity(codim, PartitionType::InteriorEntity);
    }

    const GlobalIdSet &globalIdSet () const  { return globalIdSet_; }

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

    const LeafIndexSet &leafIndexSet () const  { return leafIndexSet_; }

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










    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
    	assert(level == 0);
    	const Dune::CurvGrid::CurvLevelIterator<cd, pitype, const Grid> iterImpl(gridbase_->entityDuneIndexBegin(cd, pitype), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator(iterImpl);
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const
    {
    	assert(level == 0);
    	const Dune::CurvGrid::CurvLevelIterator<cd, pitype, const Grid> iterImpl(gridbase_->entityDuneIndexEnd(cd, pitype), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator(iterImpl);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
    	assert(level == 0);
    	const Dune::CurvGrid::CurvLevelIterator<cd, All_Partition, const Grid> iterImpl(gridbase_->entityDuneIndexBegin(cd, All_Partition), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator(iterImpl);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
    	assert(level == 0);
    	const Dune::CurvGrid::CurvLevelIterator<cd, All_Partition, const Grid> iterImpl(gridbase_->entityDuneIndexEnd(cd, All_Partition), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator(iterImpl);
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const
    {
    	const Dune::CurvGrid::CurvLevelIterator<cd, pitype, const Grid> iterImpl(gridbase_->entityDuneIndexBegin(cd, pitype), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator(iterImpl);
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafend () const
    {
    	const Dune::CurvGrid::CurvLevelIterator<cd, pitype, const Grid> iterImpl(gridbase_->entityDuneIndexEnd(cd, pitype), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator(iterImpl);
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
    	const Dune::CurvGrid::CurvLevelIterator<cd, All_Partition, const Grid> iterImpl(gridbase_->entityDuneIndexBegin(cd, All_Partition), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator(iterImpl);
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
    	const Dune::CurvGrid::CurvLevelIterator<cd, All_Partition, const Grid> iterImpl(gridbase_->entityDuneIndexEnd(cd, All_Partition), *gridbase_);
    	return typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator(iterImpl);
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
    int ghostSize( int codim ) const  { return gridbase_->nEntity(codim, Dune::PartitionType::GhostEntity); }

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
    	std::cout << "called grid.communicate()" << std::endl;
    	Dune::CurvGrid::Communication<Grid> communicator(*gridbase_, mpihelper_);
    	if (dataHandle.contains(dimension, ELEMENT_CODIM)) { communicator.template communicate<DataHandle, Data, ELEMENT_CODIM>(dataHandle, interface, direction); }
    	if (dataHandle.contains(dimension, FACE_CODIM))    { communicator.template communicate<DataHandle, Data, FACE_CODIM>(dataHandle, interface, direction); }
    	if (dataHandle.contains(dimension, EDGE_CODIM))    { communicator.template communicate<DataHandle, Data, EDGE_CODIM>(dataHandle, interface, direction); }
    	if (dataHandle.contains(dimension, VERTEX_CODIM))  { communicator.template communicate<DataHandle, Data, VERTEX_CODIM>(dataHandle, interface, direction); }
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
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityPointer        EntityPointer;
      typedef typename Dune::CurvGrid::CurvEntityPointer<EntitySeed::codimension, const Grid>  EntityPointerImpl;

      return EntityPointer(EntityPointerImpl( seed, *gridbase_ ));
    }


    /** \brief obtain Entity from EntitySeed. */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::Entity
    entity ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::Entity    Entity;
      typedef Dune::CurvGrid::CurvEntity<EntitySeed::codimension, dimension, const Grid>  EntityImpl;
      typedef Dune::CurvGrid::CurvEntitySeed<EntitySeed::codimension, const Grid>   SeedImpl;

      SeedImpl seedImpl = Grid::getRealImplementation(seed);

      return Entity(EntityImpl( seedImpl, *gridbase_ ));
    }


    /** \} */


    const GridBaseType & gridbase() const  { return *gridbase_; }





    // Auxiliary Methods only existent in CurvilinearGeometry
    // *******************************************************************

    // User sets the relative error tolerance for computing volumes of curvilinear entities
    void geometryRelativeTolerance(double tolerance)  { gridbase_->geometryRelativeTolerance(tolerance); }

    // User gets the relative error tolerance for computing volumes of curvilinear entities
    double geometryRelativeTolerance()                { return gridbase_->geometryRelativeTolerance(); }


    // Obtains global index of an entity
    template<int codim>
    GlobalIndexType entityGlobalIndex(const typename Traits::template Codim< codim >::Entity &entity)
    {
    	GlobalIndexType globalIndex;
    	bool exist = gridbase_->findEntityGlobalIndex(codim, leafIndexSet().index(entity), globalIndex);
    	assert(exist);  // Check if the index exists for self-consistency
    	return globalIndex;
    }


    // Obtains physical tag of an entity
    template<int codim>
    PhysicalTagType entityPhysicalTag(const typename Traits::template Codim< codim >::Entity &entity)
    {
    	return gridbase_->physicalTag(codim, leafIndexSet().index(entity));
    }

    // Obtain interpolation order of an entity
    template<int codim>
    InterpolatoryOrderType entityInterpolationOrder (const typename Traits::template Codim< codim >::Entity &entity)
    {
    	return gridbase_->entityInterpolationOrder(codim, leafIndexSet().index(entity));
    }

    template<int codim>
    typename Codim<codim>::EntityGeometryMappingImpl entityBaseGeometry (const typename Traits::template Codim< codim >::Entity &entity)
    {
    	return gridbase_->template entityGeometry<codim>(leafIndexSet().index(entity));
    }




  private:
    GridBaseType * gridbase_;
    MPIHelper & mpihelper_;

    //mutable std::vector< LevelIndexSet *, typename Allocator::template rebind< LevelIndexSet * >::other > levelIndexSets_;
    IndexSetImpl leafIndexSet_;
    IdSetImpl globalIdSet_;
  };



  // CurvilinearGrid::Codim
  // -------------------

  template<class ct, int cdim, bool isCached, class LogMsg>
  template< int codim >
  struct CurvilinearGrid< ct, cdim, isCached, LogMsg>::Codim
    : public Base::template Codim< codim >
  {

	typedef typename GridBaseType::template Codim<codim>::EntityGeometry  EntityGeometryMappingImpl;

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
