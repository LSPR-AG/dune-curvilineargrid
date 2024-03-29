/***************************************************************************
                          curvilineargridstorage.hh
                             -------------------
    begin                : Tue Nov 25 2014
    copyright            : (C) 2014 by Aleksejs Fomins, LSPR AG
    description          : Upgraded the mesh to curvilinear grid
***************************************************************************/

#ifndef DUNE_CURVILINEARGRIDSTORAGE_HH
#define DUNE_CURVILINEARGRIDSTORAGE_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/impl/curvilinearlooseoctree.hh>


namespace Dune {

namespace CurvGrid {


template <class ct, int cdim, bool isCached>
class CurvilinearGridStorage
{

public:

    // Grid Variable Types
    // ******************************************************************
	typedef  ct       ctype;
	typedef  int      GlobalIndexType;
	typedef  int      LocalIndexType;
	typedef  CurvilinearGeometryHelper::InternalIndexType      InternalIndexType;

	typedef  int      StructuralType;
	typedef  int      PhysicalTagType;
	typedef  int      InterpolatoryOrderType;

	static const int dimension = cdim;

	// Grid Coordinate types
	// ******************************************************************
	typedef Dune::FieldVector<ctype, cdim>         GlobalCoordinate;
	typedef Dune::FieldVector<ctype, cdim>         LocalCoordinate3D;
	typedef Dune::FieldVector<ctype, cdim-1>      LocalCoordinate2D;
	typedef Dune::FieldVector<ctype, cdim-2>      LocalCoordinate1D;


    // Entity Definition Structures
    // ******************************************************************
	struct FaceBoundaryType
	{
		enum {
			None               = 0,   // Faces that do not belong to any boundaries
			DomainBoundary     = 1,   // Faces on the boundary of computational domain
			InteriorBoundary   = 2,   // Artificial user-defined boundary that may include interior faces
			PeriodicBoundary   = 3   // Periodic boundary
		};
	};


    // Entity Key Structures
    // ******************************************************************

	// Unique GlobalId for entities when they still do not possess GlobalIndex
	// Relies on vertices having GlobalIndex from the start, provided by mesh generator
	struct CurvilinearEntityMapKey
	{
	    // This is a minimal info necessary to recognize an edge among all processes
	    // The sorted globalId's of vertices of an edge
	    struct EdgeKey
	    {
	    	GlobalIndexType node0;
	    	GlobalIndexType node1;

	        void sort() {
	            if (node0 > node1)  { std::swap(node0, node1); }
	        }

	        // Allows comparing keys. Not automatically defined by C++
	        bool operator==(const EdgeKey& A) const
	        {
	            return (A.node0 == node0) && (A.node1 == node1);
	        }

	        // Necessary to construct a set/map of keys
	        // Compare with priority on the lower registers
	        bool operator<(const EdgeKey& A) const
	        {
	            if (A.node0 == node0)  { return A.node1 < node1; }
	            else                   { return A.node0 < node0; }
	        }
	    };

	    // This is a minimal info necessary to recognize a face among all processes
	    // The sorted globalId's of vertices of an edge
	    struct FaceKey
	    {
	    	GlobalIndexType node0;
	    	GlobalIndexType node1;
	    	GlobalIndexType node2;

	        void sort()
	        {
	            if (node0 > node1)  { std::swap(node0, node1); }
	            if (node1 > node2)  { std::swap(node1, node2); }
	            if (node0 > node1)  { std::swap(node0, node1); }
	        }

	        // Allows comparing keys. Not automatically defined by C++
	        bool operator==(const FaceKey& A) const
	        {
	            return (A.node0 == node0) && (A.node1 == node1) && (A.node2 == node2);
	        }

	        // Necessary to construct a set/map of keys
	        // Compare with priority on the lower registers
	        bool operator<(const FaceKey& A) const
	        {
	            if (A.node0 == node0)
	            {
	                if (A.node1 == node1) { return A.node2 < node2; }
	                else { return A.node1 < node1; }
	            }
	            else { return A.node0 < node0; }
	        }
	    };
	};

	// Unique GlobalId for entities that possess GlobalIndex
	struct IdType
	{
	    GlobalIndexType  globalindex_;

	    IdType()  { }
	    IdType(GlobalIndexType globalindex) : globalindex_(globalindex) { }

		bool operator== ( const IdType & other) const { return globalindex_ == other.globalindex_; }
		bool operator!= ( const IdType & other) const { return !(*this == other);  }
		bool operator<  ( const IdType & other) const { return (globalindex_ < other.globalindex_); }
		bool operator>  ( const IdType & other) const { return (globalindex_ > other.globalindex_); }

		//template< class C, class T >
		//std::basic_ostream< C, T > &operator<< ( std::basic_ostream< C, T > &, const Id & );


		// Printing of the Id
		friend inline std::ostream& operator<< (std::ostream& s, const IdType & x)
		{
		    s << "(" << x.globalindex_ << ")";
		    return s;
		}
	};


    // Entity Storage Structures
    // ******************************************************************
	struct VertexStorage
	{
		GlobalIndexType                 globalIndex;
		Dune::PartitionType             ptype;
    	GlobalCoordinate  coord;
	};

	struct EdgeStorage
	{
		GlobalIndexType      globalIndex;
		Dune::PartitionType  ptype;
		LocalIndexType       elementIndex;
		InternalIndexType    subentityIndex;
	};

	// Face stores indices to 2 intersecting elements, and subentity index for the first one
	// Note: element1Index is always an internal element index
	// Note: element2Index is an internal element index only for internal faces, it is a ghost element index for process boundaries, and it is -1 for domain boundaries.
	struct FaceStorage
	{
		Dune::GeometryType   geometryType;
		GlobalIndexType      globalIndex;
		Dune::PartitionType  ptype;
		StructuralType       boundaryType;
		LocalIndexType       element1Index;
		LocalIndexType       element2Index;
		InternalIndexType    element1SubentityIndex;
		InternalIndexType    element2SubentityIndex;
		PhysicalTagType      physicalTag;
	};

	struct EntityStorage
	{
		Dune::GeometryType geometryType;
		GlobalIndexType globalIndex;
		Dune::PartitionType ptype;
		std::vector<LocalIndexType> vertexIndexSet;
		InterpolatoryOrderType interpOrder;
		PhysicalTagType physicalTag;
	};


    // Public Type Definitions
    // ******************************************************************
	//typedef GridBase           GridBaseType;

    typedef typename CurvilinearEntityMapKey::EdgeKey   EdgeKey;
    typedef typename CurvilinearEntityMapKey::FaceKey   FaceKey;

    typedef std::map<LocalIndexType, LocalIndexType>    Local2LocalMap;
    typedef std::map<GlobalIndexType, LocalIndexType>   Global2LocalMap;
    typedef typename Local2LocalMap::iterator           Local2LocalIterator;
    typedef typename Global2LocalMap::iterator          Global2LocalIterator;
    typedef typename Local2LocalMap::const_iterator           Local2LocalConstIterator;
    typedef typename Global2LocalMap::const_iterator          Global2LocalConstIterator;

    typedef std::set<LocalIndexType>                    LocalIndexSet;
    typedef typename LocalIndexSet::iterator            IndexSetIterator;

    typedef std::vector<std::vector <int> >             EntityNeighborRankVector;

    typedef Dune::CurvGrid::CurvilinearOctreeNode<ctype, cdim, isCached>   NodeType;
    typedef Dune::CurvGrid::CurvilinearLooseOctree<ctype, cdim, NodeType>   CurvilinearLooseOctree;


    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = 3;
    static const int   EDGE_CODIM     = 2;
    static const int   FACE_CODIM     = 1;
    static const int   ELEMENT_CODIM  = 0;

    static const int BOUNDARY_SEGMENT_PARTITION_TYPE = 500;


    // Curvilinear Grid User Constants
    // ******************************************************************
     double GEOMETRY_TOLERANCE;



    // Curvilinear Grid Storage Variables
    // ******************************************************************

    // Properties of the grid
    bool withGhostElements_;        // Whether ghost elements will be constructed or not
    bool withElementGlobalIndex_;   // Whether the grid will be supplied with element global indices using insertElement


    // Storage of process bounding box, since its computation is expensive
    GlobalCoordinate boundingBoxCenter_;
    GlobalCoordinate boundingBoxExtent_;

    // Storage necessary for user access and computation of globalIndex
    int nEntityTotal_[4];

    // Stores vertex coordinates and globalIds
    std::vector<VertexStorage> point_;

    // Stores all necessary data about edges, faces and elements
    // Minimalism - edges and faces do not store interpolatory vertex ids, but refer to an element subentity
    std::vector<EdgeStorage> edge_;
    std::vector<FaceStorage> face_;
    std::vector<EntityStorage> element_;

    // Storage of local subentity indices for element
    std::vector< std::vector<LocalIndexType> > elementSubentityCodim2_;   // (element_ index -> vector{edge_ index} )
    std::vector< std::vector<LocalIndexType> > elementSubentityCodim1_;   // (element_ index -> vector{face_ index} )


    // Maps from global to local indices - all entities of given codim, regardless of structural type
    Global2LocalMap entityIndexMap_[4];

    // Define unique local index for corners, to satisfy dune
    Local2LocalMap  cornerIndexMap_;  // vertex index -> corner unique index
    Local2LocalMap  cornerIndexMapRev_;  // corner unique index -> vertex index

    // Entity local index -> local structural entity index
    Local2LocalMap  boundarySegmentIndexMap_;				// Domain Boundary Face Index -> Boundary Segment Index
    Local2LocalMap  boundarySegment2LocalIndexMap_;  // Boundary Segment Index -> Domain Boundary Face Index
    Local2LocalMap  periodicBoundaryIndexMap_;	 			// Unique local index for periodic boundaries

    Local2LocalMap  processBoundaryIndexMap_[4];			// Unique local index for process boundaries
    Local2LocalMap  boundaryInternalEntityIndexMap_[4];	// Unique local index for interior entities next to communicating boundary (process and periodic)
    Local2LocalMap  ghostIndexMap_[4];								// Unique local index for ghost entities (process and periodic)


    // Index sets for entities of a specific structural type
    // Used to iterate over the grid entities
    LocalIndexSet  entityAllIndexSet_[4];
    LocalIndexSet  entityInternalIndexSet_[4];
    LocalIndexSet  entityProcessBoundaryIndexSet_[4];
    LocalIndexSet  entityGhostIndexSet_[4];
    LocalIndexSet  faceDomainBoundaryIndexSet_;
    LocalIndexSet  faceInteriorBoundaryIndexSet_;
    LocalIndexSet  facePeriodicBoundaryIndexSet_;

    // Two additional composite sets to represent Dune-specific composite partition types
    LocalIndexSet  entityDuneInteriorIndexSet_[4];         // In Dune interior entities are (internal + domain boundaries)
    LocalIndexSet  entityDuneInteriorBorderIndexSet_[4];   // In Dune interior border entities are (internal + domain + process boundaries)

    // List of all the processes sharing an entity with this process, noting the type that entity has on each process
    // BI - Boundary Internal Entity. Subentity of internal element neighbouring a PB face
    // PB - Process Boundary Entity (corners, edges, faces. Not elements)
    // G -  Ghost
    EntityNeighborRankVector BI2GNeighborRank_[4];    // boundary internal entity index -> vector{neighbor ranks}
    EntityNeighborRankVector PB2PBNeighborRank_[4];   // (entityPBIndex<codim> -> vector{neighbour ranks})
    EntityNeighborRankVector PB2GNeighborRank_[4];    // process boundary entity index -> vector{neighbor ranks}
    EntityNeighborRankVector G2BIPBNeighborRank_[4];  // ghost entity index -> vector{neighbor ranks}
    EntityNeighborRankVector G2GNeighborRank_[4];     // ghost entity index -> vector{neighbor ranks}

    // [TODO] Currently only available for FACE_CODIM. If necessary, implement for edges and faces
    EntityNeighborRankVector PERB2PERBNeighborRank_[4];  // Periodic Boundary Neighbor ranks


    // Periodic boundary storage
    std::vector<bool> periodicCuboidDimensions_;		// Defines which of the X,Y,Z directions of the cuboid domain boundary are periodic
    GlobalCoordinate periodicCuboidLength_;					// Defines the length of each dimension of the periodic cuboid
    std::vector<unsigned int> periodicFaceMatchPermutationIndexInner_;  // PeriodicFaceIndex -> Permutation index for interior intersection. When both periodic neighbor faces are permuted, then they match
    std::vector<unsigned int> periodicFaceMatchPermutationIndexOuter_;  // PeriodicFaceIndex -> Permutation index for exterior intersection.

    // Octree used to efficiently locate elements in which the points are located
    CurvilinearLooseOctree * octree_;

    // MPIHelper used for parallel communication
    MPIHelper & mpihelper_;


    // Constructor and Destructor
    // ******************************************************************
    CurvilinearGridStorage (MPIHelper & mpihelper, bool withGhostElements, bool withElementGlobalIndex, std::vector<bool> periodicCuboidDimensions) :
    	mpihelper_(mpihelper),
    	withGhostElements_(withGhostElements),
    	withElementGlobalIndex_(withElementGlobalIndex),
    	nEntityTotal_ {0, 0, 0, 0},
    	octree_(0)
    {
    		// [TODO] Hardcode constant
    		GEOMETRY_TOLERANCE = 1.0e-5;    // Default value for tolerance. Can be adjusted using gridbase.setGeometryTolerance(tolerance)

    		// Check periodicity sanity
            if (periodicCuboidDimensions.size() != 0) {
            	periodicCuboidDimensions_ = periodicCuboidDimensions;
            	assert(periodicCuboidDimensions_.size() == cdim);  // It is a coordinate periodicity vector
            	assert(periodicCuboidDimensions_[0] || periodicCuboidDimensions_[1] || periodicCuboidDimensions_[2]); // At least one of 3D must be periodic
            }
    }

    ~CurvilinearGridStorage()
    {
    	if (octree_) { delete octree_; }
    }

};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDSTORAGE_HH
