// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_GRIDFAMILY_HH
#define DUNE_CURVGRID_GRIDFAMILY_HH

#include <dune/grid/common/grid.hh>
#include <dune/curvilineargrid/curvilineargrid/capabilities.hh>
#include <dune/curvilineargrid/curvilineargrid/declaration.hh>
#include <dune/curvilineargrid/curvilineargrid/entity.hh>
#include <dune/curvilineargrid/curvilineargrid/entityseed.hh>
#include <dune/curvilineargrid/curvilineargrid/entitypointer.hh>
#include <dune/curvilineargrid/curvilineargrid/geometry.hh>
#include <dune/curvilineargrid/curvilineargrid/gridview.hh>
#include <dune/curvilineargrid/curvilineargrid/intersection.hh>
#include <dune/curvilineargrid/curvilineargrid/intersectioniterator.hh>
#include <dune/curvilineargrid/curvilineargrid/iterator.hh>
#include <dune/curvilineargrid/curvilineargrid/idset.hh>
#include <dune/curvilineargrid/curvilineargrid/indexsets.hh>

namespace Dune
{

  /** \brief namespace containing the implementations of CurvilinearGrid
   *  \ingroup CurvGrid
   */
  namespace CurvGrid
  {


    // GridFamily
    // ----------

    template <int dim, int dimworld, class ct>
    struct GridFamily
    {
    	struct Traits
    	{
    		// GridBase typedefs
    		// ****************************************************

    		typedef Dune::CurvilinearGridBase<ct, dimworld>             GridBaseType;
    		typedef Dune::CurvilinearGridStorage<ct, dimworld>          GridStorageType;

    	    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    	    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
    	    typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
    	    typedef typename GridStorageType::StructuralType            StructuralType;
    	    typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
    	    typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;
    	    typedef typename GridStorageType::Vertex                    Vertex;

  	      	typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
  	      	typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;


  	      	template <int codim>
  	      	struct BaseCodim
  	      	{
  	      		typedef typename GridStorageType::template Codim<codim>::EntityGeometry	EntityGeometry;
  	      	};

      	    // Logging Message Typedefs
      	    static const unsigned int LOG_PHASE_DEV      = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
      	    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;



    		// Grid typedefs
    		// ****************************************************

      	    typedef CurvilinearGrid< dim, dimworld, ct> Grid;

            typedef typename ct ctype;

            static const int dimension       = dim;
            static const int dimensionworld  = dimworld;


            typedef Dune::CurvGrid::HierarchicIterator<const Grid>       HierarchicIteratorImpl;

            //!   Since refinement not implemented, there is only 1 intersection type
            typedef Dune::CurvGrid::Intersection<const Grid>  LeafIntersectionImpl;
            typedef Dune::CurvGrid::Intersection<const Grid>  LevelIntersectionImpl;
            typedef Dune::CurvGrid::IntersectionIterator<const Grid>  LeafIntersectionIteratorImpl;
            typedef Dune::CurvGrid::IntersectionIterator<const Grid>  LevelIntersectionIteratorImpl;


            typedef Dune::Intersection< const Grid, LeafIntersectionImpl > LeafIntersection;
            typedef Dune::Intersection< const Grid, LevelIntersectionImpl > LevelIntersection;


            typedef Dune::IntersectionIterator< const Grid, LeafIntersectionIteratorImpl, LeafIntersectionImpl >    LeafIntersectionIterator;
            typedef Dune::IntersectionIterator< const Grid, LevelIntersectionIteratorImpl, LevelIntersectionImpl >  LevelIntersectionIterator;
            typedef Dune::EntityIterator< 0, const Grid, HierarchicIteratorImpl>                                    HierarchicIterator;

            template< int codim >
            struct Codim
            {
                typedef Dune::CurvGrid::Geometry< dimension-codim, dimensionworld, const Grid > GeometryImpl;
                typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, GeometryImpl > Geometry;

                // TODO: Maybe implement a simple linear geometry wrapper for this one to save memory
                typedef typename GeometryImpl                                                   LocalGeometry;

                //
                typedef Dune::CurvGrid::Entity<codim, dim, const Grid>      EntityImpl;
                typedef Dune::CurvGrid::EntityPointer< const Grid, codim >  EntityPointerImpl;
                typedef Dune::CurvGrid::EntitySeed< codim, const Grid >     EntitySeedImpl;

                typedef Dune::Entity<codim, dim, const Grid, EntityImpl>       Entity;
                typedef Dune::EntityPointer< const Grid, EntityPointerImpl >   EntityPointer;
                typedef Dune::EntitySeed< const Grid, EntitySeedImpl >         EntitySeed;

                template< PartitionIteratorType pitype >
                struct Partition
                {
                	typedef Dune::CurvGrid::LeafIterator< codim, const Grid >  LeafIteratorImpl;
                	typedef Dune::CurvGrid::LevelIterator< codim, const Grid > LevelIteratorImpl;

                    typedef Dune::EntityIterator< codim, const Grid, LeafIteratorImpl > LeafIterator;
                    typedef Dune::EntityIterator< codim, const Grid, LevelIteratorImpl > LevelIterator;
                };

                typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
                typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
            };

            typedef Dune::CurvGrid::IndexSet<const Grid>   LeafIndexSet;
            typedef Dune::CurvGrid::IndexSet<const Grid >  LevelIndexSet;

            typedef Dune::CurvGrid::IdSet<const Grid>  GlobalIdSet;
            typedef Dune::CurvGrid::IdSet<const Grid>  LocalIdSet;

            typedef typename GridBaseType::IdType  IdType;

            typedef typename Dune::CollectiveCommunication<MPI_Comm> CollectiveCommunication;

            template< PartitionIteratorType pitype >
            struct Partition
            {
            	typedef Dune::CurvGrid::LeafGridView<const Grid, pitype>   LeafGridViewImpl;
            	typedef Dune::CurvGrid::LevelGridView<const Grid, pitype>  LevelGridViewImpl;

            	typedef Dune::DefaultLeafGridViewTraits<const Grid, pitype>   LeafGridViewTraitsImpl;
            	typedef Dune::DefaultLevelGridViewTraits<const Grid, pitype>  LevelGridViewTraitsImpl;

                typedef Dune::GridView<LeafGridViewTraitsImpl>    LeafGridView;
                typedef Dune::GridView<LevelGridViewTraitsImpl>   LevelGridView;
            };
      };
    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_GRIDFAMILY_HH
