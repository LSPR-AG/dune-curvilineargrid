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









// Wraps the map keys necessary to construct Id for entities with minimal memory usage
struct CurvilinearEntityMapKey
{
    // This is a minimal info necessary to recognize an edge among all processes
    // The sorted globalId's of vertices of an edge
    struct EdgeKey
    {
        int node0;
        int node1;

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
        int node0;
        int node1;
        int node2;

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


// Enumerates the structural type of faces of the mesh
struct CurvilinearGridFaceType
{
  enum {
    Internal = 0,
    DomainBoundary = 1,
    ProcessBoundary = 2
  };
};


template <class ct>
class CurvilinearGridStorage
{

public:

    // Entity Storage Structures
    // ******************************************************************
	struct VertexStorage
	{
    	int globalIndex;
    	Dune::FieldVector<ct, 3> coord;
	};

	struct EdgeStorage
	{
		int globalIndex;
		int elementIndex;
		int subentityIndex;
	};

	// Face stores indices to 2 intersecting elements, and subentity index for the first one
	// Note: element1Index is always an internal element index
	// Note: element2Index is an internal element index only for internal faces, it is a ghost element index for process boundaries, and it is -1 for domain boundaries.
	struct FaceStorage
	{
		int globalIndex;
		int structuralType;
		int element1Index;
		int element2Index;
		int element1SubentityIndex;
		int physicalTag;
	};

	struct EntityStorage
	{
		Dune::GeometryType geometryType;
		int globalIndex;
		std::vector<int> vertexIndexSet;
		int interpOrder;
		int physicalTag;
	};


    // Public Type Definitions
    // ******************************************************************

    typedef Dune::FieldVector<ct, 3>               Vertex;

    typedef CurvilinearEntityMapKey::EdgeKey       EdgeKey;
    typedef CurvilinearEntityMapKey::FaceKey       FaceKey;

    typedef std::map<EdgeKey, int>                 EdgeKey2EdgeIdMap;
    typedef std::map<FaceKey, int>                 FaceKey2FaceIdMap;

    typedef Dune::CurvilinearGeometry<ct, 1, 3>    EdgeGeometry;
    typedef Dune::CurvilinearGeometry<ct, 2, 3>    FaceGeometry;
    typedef Dune::CurvilinearGeometry<ct, 3, 3>    ElementGeometry;

    typedef std::map<int, int> Index2IndexMap;

    typedef Dune::CurvilinearOctreeNode<ct>                       NodeType;
    typedef Dune::CurvilinearLooseOctree<ct, 3, NodeType>         CurvilinearLooseOctree;


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
    std::vector<EntityStorage> element_;
    std::vector<EntityStorage> ghostElement_;
    std::vector<EdgeStorage> edge_;
    std::vector<FaceStorage> face_;

    // Storage of local subentity indices for element
    std::vector< std::vector<int> > elementSubentityCodim1_;   // (element_ index -> vector<edge_ index> )
    std::vector< std::vector<int> > elementSubentityCodim2_;   // (element_ index -> vector<face_ index> )

    // Maps from global to local indices
    Index2IndexMap vertexGlobal2LocalMap_;
    Index2IndexMap edgeGlobal2LocalMap_;
    Index2IndexMap faceGlobal2LocalMap_;
    Index2IndexMap elementGlobal2LocalMap_;
    Index2IndexMap ghostGlobal2LocalMap_;

    Index2IndexMap faceInternalGlobal2LocalMap_;
    Index2IndexMap faceDomainBoundaryGlobal2LocalMap_;
    Index2IndexMap faceProcessBoundaryGlobal2LocalMap_;

    // List of all ranks of processors neighboring processorBoundaries. Index the same as processBoundaryFaceIndex_.
    std::vector<int> processBoundaryNeighborProcess_;  // (processBoundaryFaceIndex_ -> neighbor rank)

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
