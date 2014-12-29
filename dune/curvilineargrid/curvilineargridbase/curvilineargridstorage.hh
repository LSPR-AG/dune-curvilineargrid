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

    // Grid Variable Types
    // ******************************************************************
	typedef  int      GlobalIndexType;
	typedef  int      LocalIndexType;
	typedef  int      InternalIndexType;
	typedef  int      StructuralType;
	typedef  int      PhysicalTagType;
	typedef  int      InterpolatoryOrderType;


    // Entity Definition Structures
    // ******************************************************************
	struct EntityStructuralType
	{
		// Enumerates the structural type of faces of the mesh
		enum {
			Vertex = 0,                   // Any vertex

			InternalEdge = 1,             // An edge not on the process boundary
			ProcessBoundaryEdge = 2,      // Edge on the process boundary that has exactly 1 neighbor processes.
			DomainBoundaryEdge = 3,       // [ONLY 2D] Edge on the process boundary that has 0 neighbor processes.
			InternalBoundaryEdge = 4,     // [ONLY 2D] Internal edge marked by user as part of internal surface.
			ComplexBoundaryEdge = 5,      // An edge shared by more than 2 processes

			InternalFace = 6,             // Face not on the process boundary
			GhostFace = 7,                // [ONLY 2D] Face owned by other process, but connected through ProcessBoundaryEdge
			DomainBoundaryFace = 8,       // Process boundary face that has exactly 0 neighbors
			ProcessBoundaryFace = 9,      // Process boundary face that has exactly 1 neighbor
			InternalBoundaryFace = 10,    // Internal face marked by user as part of internal surface.

			InternalElement = 11,         // Any element owned by this process
			GhostElement = 12,            // Element owned by other process, but connected through ProcessBoundaryEdge

			Overlap = 13,                 // Internal Dune Magic (Just in case)
			Front = 14                    // Internal Dune Magic (Just in case)
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
		std::pair<StructuralType, GlobalIndexType>  id_;

		 bool operator== ( const IdType & other) { return ((id_.first == other.id_.first) && (id_.second == other.id_.second)); }
		 bool operator!= ( const IdType & other) { return !(*this == other);  }
		 bool operator<  ( const IdType & other) { return (id_.first == other.id_.first) ? (id_.second < other.id_.second) : (id_.first < other.id_.first); }

		 //template< class C, class T >
		 //std::basic_ostream< C, T > &operator<< ( std::basic_ostream< C, T > &, const Id & );
	};


    // Entity Storage Structures
    // ******************************************************************
	struct VertexStorage
	{
		GlobalIndexType globalIndex;
    	Dune::FieldVector<ct, cdim> coord;
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

    typedef std::map<GlobalIndexType, LocalIndexType>   Index2IndexMap;
    typedef typename Index2IndexMap::iterator           IndexMapIterator;

    typedef Dune::CurvilinearOctreeNode<ct, cdim>                 NodeType;
    typedef Dune::CurvilinearLooseOctree<ct, cdim, NodeType>      CurvilinearLooseOctree;

    template <int codim>
    struct Codim
    {
    	typedef Dune::CurvilinearGeometry<ct, cdim - codim, cdim>    EntityGeometry;
    };


    // Curvilinear Grid Storage Variables
    // ******************************************************************

    // Storage of process bounding box, since its computation is expensive
    Vertex boundingBoxCenter_;
    Vertex boundingBoxExtent_;

    // Storage necessary for user access and computation of globalIndex
    int nVertexTotal_;
    int nEdgeTotal_;
    int nFaceTotal_;
    int nElementTotal_;

    // Stores vertex coordinates and globalIds
    std::vector<VertexStorage> point_;

    // Stores all necessary data about edges, faces and elements
    // Minimalism - edges and faces do not store interpolatory vertex ids, but refer to an element subentity
    std::vector<EdgeStorage> edge_;
    std::vector<FaceStorage> face_;
    std::vector<EntityStorage> element_;

    // Storage of local subentity indices for element
    std::vector< std::vector<LocalIndexType> > elementSubentityCodim1_;   // (element_ index -> vector<edge_ index> )
    std::vector< std::vector<LocalIndexType> > elementSubentityCodim2_;   // (element_ index -> vector<face_ index> )

    // Maps from global to local indices - all entities of given codim, regardless of structural type
    Index2IndexMap vertexGlobal2LocalMap_;
    Index2IndexMap edgeGlobal2LocalMap_;
    Index2IndexMap faceGlobal2LocalMap_;
    Index2IndexMap elementGlobal2LocalMap_;

    // Structural type maps for elements
    Index2IndexMap internalElementGlobal2LocalMap_;
    Index2IndexMap ghostElementGlobal2LocalMap_;

    // Structural type maps for faces
    Index2IndexMap faceInternalGlobal2LocalMap_;
    Index2IndexMap faceDomainBoundaryGlobal2LocalMap_;
    Index2IndexMap faceProcessBoundaryGlobal2LocalMap_;

    // List of all ranks of processors neighboring processorBoundaries.
    std::map<LocalIndexType, int> processBoundaryNeighborProcess_;  // (face_ index -> neighbor rank)

    // Octree used to efficiently locate elements in which the points are located
    CurvilinearLooseOctree * octree_;


    // Constructor and Destructor
    // ******************************************************************
    CurvilinearGridStorage () :
    	nVertexTotal_(0),
    	nEdgeTotal_(0),
    	nFaceTotal_(0),
    	nElementTotal_(0),
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
