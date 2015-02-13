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
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>


namespace Dune {




template <class ct, int cdim>
class CurvilinearGridStorage
{

public:

	typedef Dune::CurvilinearGridStorage<ct, cdim>   This;

    // Grid Variable Types
    // ******************************************************************
	typedef  int      GlobalIndexType;
	typedef  int      LocalIndexType;
	typedef  Dune::CurvilinearGeometryHelper::InternalIndexType      InternalIndexType;

	typedef  int      StructuralType;
	typedef  int      PhysicalTagType;
	typedef  int      InterpolatoryOrderType;


    // Entity Definition Structures
    // ******************************************************************
	struct PartitionType
	{
		enum {
			Internal           = 0,   // Entity that is not on the process boundary
			ProcessBoundary    = 1,   // Boundary entity shared by more than one process
			DomainBoundary     = 2,   // Boundary entity that is not a process boundary [codim > 0]
			InternalBoundary   = 3,   // Artificial user-defined boundary that is not the process boundary [Not Implemented]
			FrontBoundary      = 4,   // Faces of Overlap partition [Not Implemented]
			Ghost              = 5,   // Entities stored on this process but not owned by it [all codim]
			Overlap            = 6,   // Dune-Magic [Not Implemented]
			Periodic           = 7,   // Periodic boundary [Not Implemented]
		};
	};

	const std::string PartitonTypeName[7];


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
		std::pair<StructuralType, GlobalIndexType>  id_;

		 bool operator== ( const IdType & other) { return ((id_.first == other.id_.first) && (id_.second == other.id_.second)); }
		 bool operator!= ( const IdType & other) { return !(*this == other);  }
		 bool operator<  ( const IdType & other) { return (id_.first == other.id_.first) ? (id_.second < other.id_.second) : (id_.first < other.id_.first); }

		 //template< class C, class T >
		 //std::basic_ostream< C, T > &operator<< ( std::basic_ostream< C, T > &, const Id & );


		 // Printing of the Id
		 friend inline std::ostream& operator<< (std::ostream& s, const IdType & x)
		 {
		 	s << "(" << x.id_.first << "," << x.id_.second << ")";
		 	return s;
		 }
	};


    // Entity Storage Structures
    // ******************************************************************
	struct VertexStorage
	{
		GlobalIndexType               globalIndex;
		StructuralType                structuralType;
    	Dune::FieldVector<ct, cdim>   coord;
	};

	struct EdgeStorage
	{
		GlobalIndexType    globalIndex;
		StructuralType     structuralType;
		LocalIndexType     elementIndex;
		InternalIndexType  subentityIndex;
	};

	// Face stores indices to 2 intersecting elements, and subentity index for the first one
	// Note: element1Index is always an internal element index
	// Note: element2Index is an internal element index only for internal faces, it is a ghost element index for process boundaries, and it is -1 for domain boundaries.
	struct FaceStorage
	{
		Dune::GeometryType geometryType;
		GlobalIndexType   globalIndex;
		StructuralType    structuralType;
		LocalIndexType    element1Index;
		LocalIndexType    element2Index;
		InternalIndexType element1SubentityIndex;
		PhysicalTagType   physicalTag;
	};

	struct EntityStorage
	{
		Dune::GeometryType geometryType;
		GlobalIndexType globalIndex;
		StructuralType structuralType;
		std::vector<LocalIndexType> vertexIndexSet;
		InterpolatoryOrderType interpOrder;
		PhysicalTagType physicalTag;
	};


    // Public Type Definitions
    // ******************************************************************

    typedef Dune::FieldVector<ct, cdim>                 Vertex;

    typedef typename CurvilinearEntityMapKey::EdgeKey   EdgeKey;
    typedef typename CurvilinearEntityMapKey::FaceKey   FaceKey;

    typedef std::map<LocalIndexType, LocalIndexType>    Local2LocalMap;
    typedef std::map<GlobalIndexType, LocalIndexType>   Global2LocalMap;
    typedef typename Local2LocalMap::iterator           Local2LocalIterator;
    typedef typename Global2LocalMap::iterator          Global2LocalIterator;

    typedef std::set<LocalIndexType>                    LocalIndexSet;
    typedef typename LocalIndexSet::iterator            IndexSetIterator;

    typedef std::vector<std::vector <int> >             EntityNeighborRankVector;

    typedef Dune::CurvilinearOctreeNode<ct, cdim>                 NodeType;
    typedef Dune::CurvilinearLooseOctree<ct, cdim, NodeType>      CurvilinearLooseOctree;

    template <int codim>
    struct Codim
    {
    	typedef Dune::CachedCurvilinearGeometry<ct, cdim - codim, cdim>    EntityGeometry;
    };


    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = 3;
    static const int   EDGE_CODIM     = 2;
    static const int   FACE_CODIM     = 1;
    static const int   ELEMENT_CODIM  = 0;



    // Curvilinear Grid Storage Variables
    // ******************************************************************

    // Properties of the grid
    bool withGhostElements_;   // Whether ghost elements will be constructed or not


    // Storage of process bounding box, since its computation is expensive
    Vertex boundingBoxCenter_;
    Vertex boundingBoxExtent_;

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

    // Entity local index -> local structural entity index
    Local2LocalMap  boundarySegmentIndexMap_;  // This one only for Domain Boundary Faces

    Local2LocalMap  processBoundaryIndexMap_[4];
    Local2LocalMap  boundaryInternalEntityIndexMap_[4];
    Local2LocalMap  ghostIndexMap_[4];

    // Index sets for entities of a specific structural type
    // Used to iterate over the grid entities
    LocalIndexSet  entityAllIndexSet_[4];
    LocalIndexSet  entityInternalIndexSet_[4];
    LocalIndexSet  entityProcessBoundaryIndexSet_[4];
    LocalIndexSet  entityDomainBoundaryIndexSet_[4];
    LocalIndexSet  entityGhostIndexSet_[4];

    // Two additional composite sets to represent Dune-specific composite partition types
    LocalIndexSet  entityDuneInteriorIndexSet_[4];         // In Dune interior entities are (internal + domain boundaries)
    LocalIndexSet  entityDuneInteriorBorderIndexSet_[4];   // In Dune interior border entities are (internal + domain + process boundaries)

    // List of all the processes sharing an entity with this process, noting the type that entity has on each process
    // BI - Boundary Internal Entity. Subentity of internal element neighbouring a PB face
    // PB - Process Boundary Entity (corners, edges, faces. Not elements)
    // G -  Ghost
    EntityNeighborRankVector BI2GNeighborRank_[4];    // boundary internal entity index -> vector{neighbor ranks}
    EntityNeighborRankVector PB2PBNeighborRank_[4];  // (entityPBIndex<codim> -> vector{neighbour ranks})
    EntityNeighborRankVector PB2GNeighborRank_[4];    // process boundary entity index -> vector{neighbor ranks}
    EntityNeighborRankVector G2BIPBNeighborRank_[4];  // ghost entity index -> vector{neighbor ranks}
    EntityNeighborRankVector G2GNeighborRank_[4];     // ghost entity index -> vector{neighbor ranks}



    // Octree used to efficiently locate elements in which the points are located
    CurvilinearLooseOctree * octree_;


    // Constructor and Destructor
    // ******************************************************************
    CurvilinearGridStorage (bool withGhostElements) :
    	withGhostElements_(withGhostElements),
    	nEntityTotal_ {0, 0, 0, 0},
    	PartitonTypeName { "Internal", "ProcessBoundary", "DomainBoundary", "InternalBoundary", "FrontBoundary", "Ghost", "Overlap" },
    	octree_(0)
    {

    }

    ~CurvilinearGridStorage()
    {
    	delete octree_;
    }

};



} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDSTORAGE_HH
