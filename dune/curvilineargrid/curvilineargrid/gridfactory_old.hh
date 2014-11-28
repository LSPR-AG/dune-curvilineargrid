// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_CURVGRID_GRIDFACTORY_HH
#define DUNE_CURVGRID_GRIDFACTORY_HH

/** \file
 *  \author Aleksejs Fomins
 *  \brief  Implementation of Curvilinear Grid Factory
 */

#include <config.h>

#include <map>
#include <limits>
#include <vector>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundaryprojection.hh>

//#include <dune/grid/utility/globalindex.hh>

#include <dune/alugrid/common/transformation.hh>
#include <dune/alugrid/3d/alugrid.hh>
#include <dune/alugrid/3d/gridfactory.hh>
#include <dune/alugrid/3d/gridfactory.cc>



#include <parmetis.h>






namespace Dune
{


template< class ALUGrid >
class CurvilinearGridFactory
{
	/** \brief Factory class for 3d ALUGrids */

  public:
    typedef ALUGrid HostGrid;

    typedef typename HostGrid::ctype ctype;

    typedef unsigned int VertexGlobalId;
    typedef unsigned int VertexLocalIndex;
    typedef unsigned int CornerLocalIndex;

    typedef unsigned int ElementGlobalId;
    typedef unsigned int ElementLocalIndex;
    typedef unsigned int BoundaryLocalIndex;

    static const unsigned int dimension = HostGrid::dimension;
    static const unsigned int dimensionworld = HostGrid::dimensionworld;

    template< int codim >
    struct Codim
    {
      typedef typename HostGrid::template Codim< codim >::Entity Entity;
      typedef typename HostGrid::template Codim< codim >::EntityPointer EntityPointer;
    };


  private:
    bool verbose = true;

    struct ElementData
    {
    	GeometryType geometryType;
    	ElementGlobalId globalId;
    	int interpOrder;
    	std::vector< VertexLocalIndex > vertexIndex;
    };

    typedef FieldVector< ctype, dimensionworld > VertexCoordinate;
    typedef std::vector<VertexLocalIndex > ElementVertexIndex;
    typedef std::vector< std::pair< VertexCoordinate, VertexGlobalId > > VertexVector;
    typedef std::vector< ElementData > ElementVector;
    typedef std::vector<BoundaryLocalIndex > ElementBoundaryIndex;

    typedef std::map< VertexGlobalId, VertexLocalIndex > VertexGlobal2LocalMap;
    typedef std::map< VertexGlobalId, CornerLocalIndex > VertexGlobal2CornerLocalMap;
    typedef std::map< ElementLocalIndex, std::vector<BoundaryLocalIndex> > Element2BoundaryLinker;

    typedef Dune::ReferenceElement< ctype, dimension > ReferenceElement;
    typedef Dune::ReferenceElements< ctype, dimension > ReferenceElements;
    typedef Dune::ReferenceElements< ctype, dimension-1 > SubReferenceElements;

    typedef std::vector<VertexLocalIndex >              BoundaryKey;
    typedef std::vector< BoundaryKey > 					BoundaryKeyVector;
    typedef std::map< BoundaryKey, BoundaryLocalIndex > ProcessBorderMap;

    typedef std::pair<unsigned, unsigned>   IndexToProcess;




    // Auxiliary structures
    VertexGlobal2LocalMap       vertexIndexMap_;
    VertexGlobal2CornerLocalMap cornerIndexMap_;
    int cornercount_;

    // Data
    VertexVector vertices_;
    ElementVector elements_;
    ElementVector boundarySegments_;
    Element2BoundaryLinker elementBoundaries_;

    // Parallel Implementation
    MPIHelper &mpihelper_;
    static const int MASTER_RANK = 0;






  public:


    // TODO: Figure out the meaning of the first argument "const bool realGrid" in alugrid constructor
    CurvilinearGridFactory( MPIHelper &mpihelper ) :
      mpihelper_(mpihelper),
      cornercount_(0)
    {}

    ~CurvilinearGridFactory ()  {}


    // Insert vertex with globalID specification
    // NOTE: Vertices will be inserted into HostGrid at the InsertElement stage
    // NOTE: Has protection from double-insertion, because may be added multiple times during partition phase
    void insertVertex ( const VertexCoordinate &pos, const VertexGlobalId globalId )
    {
    	print_debug("Factory insertVertex localId = " + std::to_string(vertices_.size()) + ", globalID = " + std::to_string(globalId) + " with coordinates = " + vector2string(pos));

		if (vertexIndexMap_.find(globalId) == vertexIndexMap_.end())  {
			vertexIndexMap_[globalId] = vertices_.size();
			vertices_.push_back( std::make_pair( pos, globalId ) );
		}
    }

    // Store curvilinear element, insert linearised element to HostGridFactory
    // NOTE: Vertices given in globalID, converted intrinsically into
    //   - vertexInsertionIndex for curvilinear elements
    //   - cornerInsertionIndex for linearised elements
    void insertElement(
      GeometryType &geometry,
      const int globalId,
      const std::vector< VertexLocalIndex > &elementVertexId,
      const int elemOrder)
    {
        ElementData thisElement;

        print_debug("Factory insertElement localID = " + std::to_string(elements_.size()) + ", globalID = " + std::to_string(globalId) + " with vertices = " + vector2string(elementVertexId));

        thisElement.geometryType = geometry;
        thisElement.globalId = globalId;
        thisElement.interpOrder = elemOrder;
        thisElement.vertexIndex = elementVertexId;

        // Insert element into CurvGrid and HostGrid
        elements_.push_back(thisElement);
    }


    // Insert boundary as face of some element
    // void insertBoundary ( const int element, const int face, const int id )  { hostgridfactory_.insertBoundary(element, face, id );  }


    // Store curvilinear boundarySegment, insert linearised segment to HostGridFactory
    // NOTE: Vertices given in globalID, converted intrinsically into
    //   - vertexInsertionIndex for curvilinear elements
    //   - cornerInsertionIndex for linearised elements
    // NOTE: Requires the element this boundary is linked to
    void insertBoundarySegment(
        GeometryType &geometry,
        const int globalId,
        const std::vector< VertexGlobalId > &boundaryVertexId,
        const int elemOrder,
        const ElementLocalIndex linkedElement)
    {
        ElementData thisBoundarySegment;

        print_debug("Factory insertBoundarySegment localID = " + std::to_string(boundarySegments_.size()) + ", globalID = " + std::to_string(globalId) + " with vertices = " + vector2string(boundaryVertexId));

        thisBoundarySegment.geometryType = geometry;
        thisBoundarySegment.globalId = globalId;
        thisBoundarySegment.interpOrder = elemOrder;
        thisBoundarySegment.vertexIndex = boundaryVertexId;

        // Link a corresponding internal element to this boundary
        Element2BoundaryLinker::iterator linkedIter = elementBoundaries_.find(linkedElement);

        // If this element has no boundaries specified yet, add a new boundary vector
        if (linkedIter != elementBoundaries_.end()) { elementBoundaries_[linkedElement] = ElementBoundaryIndex(1, boundarySegments_.size());    print_debug("eeee-a"); }
        else
        {
        	// Otherwise extract the existing boundaries, add a new one, and push back to the map
			//ElementBoundaryIndex thisElementLinked = linkedIter->second;
			ElementBoundaryIndex thisElementLinked = elementBoundaries_[linkedElement];
			thisElementLinked.push_back(boundarySegments_.size());
			elementBoundaries_[linkedElement] = thisElementLinked;
        }

        // Store boundary segments
        boundarySegments_.push_back(thisBoundarySegment);
    }

    //insertBoundarySegment ( const std::vector< unsigned int >& vertices, const shared_ptr<BoundarySegment<3,3> >& boundarySegment )  {  }


    //CurvilinearGrid* createGrid ()
    void createGrid()
    {
    	// 1) Compute the partitioning of the elements among processes
      	// ***********************************************************
    	unsigned numElements = elements_.size();
    	unsigned numBSegments = boundarySegments_.size();

    	std::vector<unsigned> part(numElements);
    	//partition_compute(part);
    	//print_debug("Partition computation result: " + vector2string(part));

    	// 2) Communicate corresponding elements, create maps
    	// ***********************************************************
    	//partition_communicate(part);
    	//print_debug("Partitioning communicated");

    	// 3) Construct HostGrid
    	// ***********************************************************
    	HostGrid* Grid;
    	{
    		Dune::GridFactory<HostGrid> hostgridfactory;

    		// 3.1) Insert elements and associated vertices
    		for (int i = 0; i < numElements; i++)   { insertHostElement(hostgridfactory, elements_[i], false); }
    		print_debug("Inserted elements to HostGridFactory");

    		// 3.2) Insert boundarySegments and associated vertices
    		for (int i = 0; i < numBSegments; i++)  { insertHostElement(hostgridfactory, boundarySegments_[i], true); }
    		print_debug("Inserted boundary segments to HostGridFactory");

    		// 3.3) Compute and insert process boundaries
    		computeAndInsertProcessBorders(hostgridfactory);
    		print_debug("Inserted processBorders to HostGridFactory");

    		// 3.4) Create HostGrid
    		// According to Andreas Dedner, ALUGrid* createGrid ( const bool addMissingBoundaries, const std::string dgfName = "" );
    		// We do not really want to add missing boundaries as we hope to have inserted all of them
    		HostGrid* Grid = hostgridfactory.createGrid(false);
    		print_debug("Created HostGrid");
    	}

    	// Create CurvilinearGrid

    }







    // ????????????????????????????????????
    virtual unsigned int
    insertionIndex ( const typename HostGrid::LeafIntersection &intersection ) const
    {

    }

    // ????????????????????????????????????
    virtual bool
    wasInserted ( const typename HostGrid::LeafIntersection &intersection ) const
    {

    }

  private:

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::string tmp_str;
        for (int i = 0; i < V.size(); i++) { tmp_str += std::to_string(V[i]) + " "; }
        return tmp_str;
    }

    // Writes debug info to the command line
    // TODO: Use IFDEF to manipulate between no output, all output, or only master process output
    void print_debug(std::string s)
    {
        if (verbose) { std::cout << "Process_" << mpihelper_.rank() << ": " << s << std::endl; }
    }


    // Get GlobalID of a vertex
    size_t globalId ( const VertexLocalIndex &id ) const  {  return vertices_[ id ].second; }

    // Get Coordinate of a vertex
    const VertexCoordinate &position ( const VertexLocalIndex &id ) const  {  return vertices_[ id ].first;  }




    // Inserts corners into HostGridFactory. Uses global-to-local map to check if already inserted
    // Then inserts element into HostGridFactory using cornerInsertionIndex
    void insertHostElement(Dune::GridFactory<HostGrid> & hostgridfactory, ElementData elem, bool is_bsegment) {
    	// Compute corners, add to renumber, insert to HostGridFactory

    	print_debug("Adding element globalId=" + std::to_string(elem.globalId) + " to HostGridFactory as " + (is_bsegment ? "boundarySegment" : "internal element"));


    	std::vector< VertexLocalIndex > cornerVertexIndex = getCornerInfo(elem.geometryType, elem.interpOrder, elem.vertexIndex);
        std::vector< CornerLocalIndex > cornerIndex;

        if ((is_bsegment && (elem.geometryType.dim() != dimensionworld - 1)) || (!is_bsegment && (elem.geometryType.dim() != dimensionworld)))
        {
        	DUNE_THROW(Dune::IOError, "Curvilinear Grid Factory: Unexpected geometry type dimension");
        }

        for (int i = 0; i < cornerVertexIndex.size(); i++)
        {
        	// If this corner has already been inserted into HostGrid, reuse inserted vertex
        	// Otherwise insert new vertex to hostgrid and update map

        	VertexGlobal2CornerLocalMap::iterator cornerpos = cornerIndexMap_.find(cornerVertexIndex[i]);

            if (cornerpos == cornerIndexMap_.end())
        	  {
        		  // Insert corner into HostGridFactory
            	  hostgridfactory.insertVertex(
        				  vertices_[cornerVertexIndex[i]].first,
        				  vertices_[cornerVertexIndex[i]].second);

        		  cornerIndex.push_back(cornercount_);
        		  cornerIndexMap_[cornerVertexIndex[i]] = cornercount_++;
        	  } else {
        		  cornerIndex.push_back(cornerpos->second);
        	  }
        }

        // Insert Element or a boundarySegment
        if (!is_bsegment)  { hostgridfactory.insertElement(elem.geometryType, cornerIndex); }
        else               { hostgridfactory.insertBoundarySegment(cornerIndex);  }
    }

    // Gets corner information from Dune-style curvilinear vertex array
    // Templated because information can mean {id's, indices, vertex coordinates, etc}
    template <class T>
    std::vector<T> getCornerInfo(GeometryType & geomType, const int order, const std::vector<T> & vertexInfo)
    {
    	std::vector<T> rez;

    	// Get corner number of this element
    	int cornerNumber;
    	if (geomType.dim() == 2)  { cornerNumber = SubReferenceElements::general(geomType).size( geomType.dim() ); }
    	else                      { cornerNumber = ReferenceElements::general(geomType).size( geomType.dim() ); }

    	for (int i = 0; i < cornerNumber; i++) {
    		int cornerCurvIndex = Dune::CurvilinearElementInterpolator<ctype, dimension, dimensionworld>::cornerID(geomType, order, i);
    		rez.push_back(vertexInfo[cornerCurvIndex]);
    	}
    	return rez;
    }


    // Returns d-1 subentity corner local indices, sorted
    std::vector<int> GEOMETRY_ElementSubentityCornerInternalIndices(GeometryType & gt, const int ind)
    {
        std::vector<int> rez;
        if (gt.isTriangle())
        {
            switch (ind)
            {
            case 0 :  rez = std::vector<int> {0, 1};  break;
            case 1 :  rez = std::vector<int> {1, 2};  break;
            case 2 :  rez = std::vector<int> {0, 2};  break;
            default : DUNE_THROW(Dune::IOError, "GMSH Reader: Wrong input arguments for SubentityCorners " );
            }
        }
        else if (gt.isTetrahedron())
        {
            switch (ind)
            {
            case 0 :  rez = std::vector<int> {0, 1, 2};  break;
            case 1 :  rez = std::vector<int> {0, 1, 3};  break;
            case 2 :  rez = std::vector<int> {0, 2, 3};  break;
            case 3 :  rez = std::vector<int> {1, 2, 3};  break;
            default : DUNE_THROW(Dune::IOError, "GMSH Reader: Wrong input arguments for SubentityCorners " );
            }
        } else  {  DUNE_THROW(Dune::IOError, "GMSH Reader: Not implemented element subentityNo for this element type " ); }
        return rez;
    }

    // Takes all interpolatory vertex indices, extracts corner indices, sorts them in ascending order
    BoundaryKey GEOMETRY_MakeElementCornerKey(
    		GeometryType & gt,
    		const int interpOrder,
    		const ElementVertexIndex & vertexIndex)
	{
    	// 1) Get corners
    	BoundaryKey key = getCornerInfo(gt, interpOrder, vertexIndex);

    	// 2) sort by increasing index
        std::sort(key.begin(), key.end());

        return key;
	}


    BoundaryKeyVector GEOMETRY_MakeElementSubentityCornerKeys(
    		GeometryType & gt,
    		const int interpOrder,
    		const ElementVertexIndex & vertexIndex)
	{
    	BoundaryKeyVector keys(4);

    	// 1) Get corners
    	std::vector<VertexLocalIndex > corners = getCornerInfo(gt, interpOrder, vertexIndex);

    	// 2) sort by increasing index
        std::sort(corners.begin(), corners.end());

        // 3) Get this element subentity number
		int thisElmSubentities = ReferenceElements::general(gt).size(1);

		// 4) Create key for each codim-1 subentity of the element
        for (int iSub = 0; iSub < thisElmSubentities; iSub++)
        {
        	// 4.1) get all subentities associated with this element type
        	std::vector<int> this_subentity_ind = GEOMETRY_ElementSubentityCornerInternalIndices(gt, iSub);
        	for (int iCoord = 0; iCoord < this_subentity_ind.size(); iCoord++) { keys[iSub].push_back(corners[this_subentity_ind[iCoord]]); }
        }

        return keys;
	}

    // Calculates process borders by finding all subentities of all elements,
    // Subtracting those that count twice, and those which match with boundarySegments
    void computeAndInsertProcessBorders(Dune::GridFactory<HostGrid> & hostgridfactory)
    {
    	// Make an std::set of faces (face = indices of all facial vertices sorted in ascending order)
    	ProcessBorderMap processFaces;

    	// Go over all elements on this process
    	for (int i = 0; i < elements_.size(); i++)
    	{
    		print_debug("Making subentity keys for element " + std::to_string(i) + " given by: " + vector2string(elements_[i].vertexIndex));

        	// 1) Get index sets which correspond to all faces of this element
    		BoundaryKeyVector subentityKeys = GEOMETRY_MakeElementSubentityCornerKeys(elements_[i].geometryType, elements_[i].interpOrder, elements_[i].vertexIndex);

    		// 2) Add all faces on the boundary of this process to a map
            for (int iSub = 0; iSub < subentityKeys.size(); iSub++)
            {
            	print_debug("* Adding boundary key " + vector2string(subentityKeys[iSub]));

            	// Check if this face has been seen already
            	auto faceit = processFaces.find(subentityKeys[iSub]);

                // If this boundary not added yet, add it
                if (faceit == processFaces.end())  { processFaces[subentityKeys[iSub]] = i; }
                // If this boundary added twice, it is not a processBoundary so erase it
                else                               { processFaces.erase(faceit); }
            }
    	}

    	// 3) Subtract from the map all faces which are on the domain boundary
        for (int i = 0; i < boundarySegments_.size(); i++)
        {
        	print_debug("Making boundary key for boundary given by: " + vector2string(boundarySegments_[i].vertexIndex));

    		// 3.1) Make key of corners of this boundary
    		BoundaryKey key = GEOMETRY_MakeElementCornerKey(boundarySegments_[i].geometryType, boundarySegments_[i].interpOrder, boundarySegments_[i].vertexIndex);

    		// 3.2) Now delete this face
    		auto faceit = processFaces.find(key);

    		if (faceit == processFaces.end()) { print_debug("Error: Have not found boundary face " + vector2string(key)); }

    		processFaces.erase(faceit);
        }

        print_debug("Writing process boundaries");


        // 4) Loop over all process boundaries and insert them
        ProcessBorderMap::iterator procMapEnd = processFaces.end();
        for (ProcessBorderMap::iterator it=processFaces.begin(); it!=procMapEnd; ++it) {
        	hostgridfactory.insertProcessBorder(it->first);
        }
    }

    // Compute which partner to cling with in this round
    int ClinkAlgorithmFriend(int rank, int s, int round)
    {
    	if (rank == 0) { return s - round; }
    	else           { return rank + s - 3 - 2 * ((rank + round - 1) % (s - 1)); }
    }

    // Order of sorting Index-Process pairs, sort ascending wrt index
    static bool indexToProcessOrder(IndexToProcess a, IndexToProcess b)  { return a.first < b.first; }

    // Binary search for finding beginning of elements which correspond to requested process
    unsigned indexToProcessBinaryFindFirst(std::vector<IndexToProcess> & part_sorted, unsigned p)
    {
    	int l = 0;
    	int r = part_sorted.size() - 1;

    	while (r-l > 1) {
    		int mid = (l + r) / 2;
    		if (part_sorted[mid].first < p) { l = mid; } else { r = mid; }
    	}

    	return r;
    }

    // TODO: At the moment the boundarySegment must share vertices with some element
    void partition_communicate(std::vector<unsigned> & part)
    {
    	// 1) Sort communicated element indices by process
    	// ****************************************************
    	std::vector<IndexToProcess> part_sorted;
    	for (int i = 0; i < part.size(); i++) { part_sorted.push_back(IndexToProcess(part[i], i)); }
    	std::sort(part_sorted.begin(), part_sorted.end(), indexToProcessOrder);

    	// 2) Create new arrays for all the important information
    	// ****************************************************
        VertexVector vertices_tmp;
        ElementVector elements_tmp;
        ElementVector boundarySegments_tmp;
        Element2BoundaryLinker elementBoundaries_tmp;

        vertices_.swap(vertices_tmp);
        elements_.swap(elements_tmp);
        boundarySegments_.swap(boundarySegments_tmp);
        std::swap(elementBoundaries_, elementBoundaries_tmp);


        // 3) Fill in data to stay on this process, adjust new local vertex index
        // ****************************************************
        int rank = mpihelper_.rank();
        int size = mpihelper_.size();
        int thisDataInd = indexToProcessBinaryFindFirst(part_sorted, rank);

        // Nullify the vertex map, because we are inserting all vertices from scratch again
        vertexIndexMap_ = VertexGlobal2LocalMap();

        while (part_sorted[thisDataInd].first == rank)
        {
        	// 3.1) Get element
        	int thisElmLocalIndexOld = part_sorted[thisDataInd].second;
        	ElementData thisElem = elements_tmp[thisElmLocalIndexOld];

        	// 3.2) Get its vertex coordinates, insert them to the updated factory, update local vertex index
        	for (int i = 0; i < thisElem.vertexIndex.size(); i++)
        	{
        		std::pair< VertexCoordinate, VertexGlobalId > thisVertexPair = vertices_tmp[thisElem.vertexIndex[i]];
        		insertVertex(thisVertexPair.first, thisVertexPair.second);

        		thisElem.vertexIndex[i] = vertexIndexMap_[thisVertexPair.second];
        	}

        	// 3.3) Insert element to the updated factory
        	insertElement(thisElem.geometryType, thisElem.globalId, thisElem.vertexIndex, thisElem.interpOrder);


        	// 3.4) Insert linked boundary segments to updated factory
        	std::vector<BoundaryLocalIndex> thisElmBoundaryLocalIndicesOld = elementBoundaries_tmp[thisElmLocalIndexOld];

        	for (int i = 0; i < thisElmBoundaryLocalIndicesOld.size(); i++)
        	{
        		// 3.4.1) Get boundary segment
        		ElementData thisBoundary = boundarySegments_tmp[thisElmBoundaryLocalIndicesOld[i]];

            	// 3.4.2) Get its vertex coordinates, update their local vertex index
            	for (int j = 0; j < thisBoundary.vertexIndex.size(); j++)
            	{
            		VertexLocalIndex vertexIndexOld = thisBoundary.vertexIndex[j];
            		VertexGlobalId vertexGlobalId = vertices_tmp[vertexIndexOld].second;
            		VertexLocalIndex vertexIndexNew = vertexIndexMap_[vertexGlobalId];

            		thisBoundary.vertexIndex[j] = vertexIndexNew;
            	}

            	// 3.4.3) Insert boundary segment to the updated factory
                insertBoundarySegment(
                		thisBoundary.geometryType,
                		thisBoundary.globalId,
                		thisBoundary.vertexIndex,
                		thisBoundary.interpOrder,
                		elements_.size() - 1);
        	}
        }


        // 4) Communicate-Fill data from other processes
        // ****************************************************
        int comm_stage_max = size % 2 ? size - 1 : size;
        for (int iStage = 0; iStage < comm_stage_max; iStage++ )
        {
        	// 4.1) Find out which process will this be communicating with at this stage
        	// -------------------------------------------------------------------------
        	int rank_friend = ClinkAlgorithmFriend(rank, comm_stage_max, iStage);

        	// 4.2) Prepare data for communication: elements, boundary segments, vertices, boundarylinkers
        	// TODO: For communicated vertices use GlobalVertexID vectors
        	// -------------------------------------------------------------------------
        	int pDataInd = indexToProcessBinaryFindFirst(part_sorted, rank_friend);

        	VertexGlobal2LocalMap vertexCommMap;

        	// TODO: comm_send and comm_recv arrays come here
            VertexVector vertices_comm;
            ElementVector elements_comm;
            ElementVector boundarySegments_comm;
            std::vector<unsigned> boundaryLinker_comm;

            VertexVector vertices_recv;
            ElementVector elements_recv;
            ElementVector boundarySegments_recv;
            std::vector<unsigned> boundaryLinker_recv;

            while (part_sorted[pDataInd].first == rank_friend)
            {
            	// 4.2.1) Get element
            	int thisElmIndex = part_sorted[pDataInd++].second;
            	ElementData thisElem = elements_tmp[thisElmIndex];

            	for (int iVert = 0; iVert < thisElem.vertexIndex.size(); iVert++)
            	{
            		// 4.2.2) Add vertices to communication renumber, then to communication array
            		int globalId = vertices_tmp[thisElem.vertexIndex[iVert]].second;
            		if (vertexCommMap.find(globalId) == vertexCommMap.end())
            		{
            			vertices_comm.push_back(vertices_tmp[thisElem.vertexIndex[iVert]]);
            			vertexCommMap[globalId] = iVert;
            		}

            		// 4.2.3) Change all vertexIndex by globalID
            		thisElem.vertexIndex[iVert] = globalId;
            	}

            	// 4.2.4) Add element to communication array
            	elements_comm.push_back(thisElem);

            	// 4.2.5) Get Boundary Segments, change vertexIndex too, add to comm array
            	std::vector<BoundaryLocalIndex> thisElmBoundaryLocalIndicesOld = elementBoundaries_tmp[thisElmIndex];
            	for (int iBound = 0; iBound < thisElmBoundaryLocalIndicesOld.size(); iBound++)
            	{
            		// 4.2.5.1) Get boundary segment
            		ElementData thisBoundary = boundarySegments_tmp[thisElmBoundaryLocalIndicesOld[iBound]];


                	// 4.2.5.2) Get its vertex coordinates, update local vertex index
                	for (int iVert = 0; iVert < thisBoundary.vertexIndex.size(); iVert++)
                	{
                		int globalId = vertices_tmp[thisBoundary.vertexIndex[iVert]].second;
                		thisBoundary.vertexIndex[iVert] = globalId;
                	}

                	// 4.2.5.3) add boundary segment to the comm array
                	boundarySegments_comm.push_back(thisBoundary);

                	// 4.2.5.4) Link from this boundary to the element it belongs to, localProcCommIndex
                	boundaryLinker_comm.push_back(elements_comm.size() - 1);
            	}
            }


        	// 4.3) Communicate arrays
            // -------------------------------------------------------------------------
        	mpiExchange(
        			rank_friend,
        			vertices_comm, elements_comm, boundarySegments_comm, boundaryLinker_comm,
        			vertices_recv, elements_recv, boundarySegments_recv, boundaryLinker_recv);


        	// 4.4) Insert received data
        	// -------------------------------------------------------------------------
        	std::vector<unsigned> localCom2localProcElem;

        	// 4.4.1) Insert vertices
        	for (int iVert = 0; iVert < vertices_recv.size(); iVert++) { insertVertex(vertices_recv[iVert].first, vertices_recv[iVert].second);  }

        	// 4.4.2) Insert elements
        	for (int iElem = 0; iElem < elements_recv.size(); iElem++) {
        		insertElement(elements_recv[iElem].geometryType, elements_recv[iElem].globalId, elements_recv[iElem].vertexIndex, elements_recv[iElem].interpOrder);
        		localCom2localProcElem[iElem] = elements_.size() - 1;
        	}

        	// 4.4.3) Inseert boundary segments
        	for (int iBound = 0; iBound < boundarySegments_recv.size(); iBound++) {
        		insertBoundarySegment(
        				boundarySegments_recv[iBound].geometryType,
        				boundarySegments_recv[iBound].globalId,
        				boundarySegments_recv[iBound].vertexIndex,
        				boundarySegments_recv[iBound].interpOrder,
        				localCom2localProcElem[boundaryLinker_recv[iBound]]
        				);
        	}

        }
    }

    // Sends all comm arrays to friend, receives all recv arrays
    void mpiExchange(
    		int rank_friend,
    		VertexVector & vertices_comm, ElementVector & elements_comm, ElementVector & boundarySegments_comm, std::vector<unsigned> & boundaryLinker_comm,
    		VertexVector & vertices_recv, ElementVector & elements_recv, ElementVector & boundarySegments_recv, std::vector<unsigned> & boundaryLinker_recv)
    {
    	// MPI_RECV/SEND
    	// Or, even better MPI_Sendrecv for exchange
    }


    /** \brief Create an initial partitioning of a Dune grid, i.e., not taking into account communication cost
     *
     * This pipes a Dune grid into the method ParMETIS_V3_PartMeshKway (see the ParMetis documentation for details)
     *
     * \param gv The grid view to be partitioned
     * \param mpihelper The MPIHelper object, needed to get the MPI communicator
     *
     * \return std::vector with one uint per All_Partition element.  For each element, the entry is the
     *    number of the partition the element is assigned to.
     */
     // Partition mesh using ParMETIS
     void partition_compute(std::vector<unsigned> & part) {

        // ****************************************************
        // Preliminaries
        // ****************************************************
#if PARMETIS_MAJOR_VERSION < 4
      typedef idxtype idx_t;
      typedef float real_t;
#endif

      int elementNumber = elements_.size();
      int elementDim = elements_[0].geometryType.dim();
      int elementFaceCorners = ReferenceElements::general( elements_[0].geometryType ).size(0, 1, elementDim);

      // ****************************************************
      // Setup parameters for ParMETIS
      // ****************************************************
      idx_t wgtflag = 2;                                  // We use different weights for each element
      idx_t numflag = 0;                                  // we are using C-style arrays
      idx_t ncon = 1;                                     // number of balance constraints
      idx_t ncommonnodes = elementFaceCorners;            // number of nodes elements must have in common in order to be adjacent to each other
      idx_t nparts = mpihelper_.size();                   // number of parts equals number of processes
      std::vector<real_t> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
      std::vector<real_t> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)
      idx_t options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
      idx_t edgecut;                                      // will store number of edges cut by partition

      // ****************************************************
      // Communicate the number of elements on each process
      // ****************************************************
      std::vector<idx_t> elmdist;

      int* elmdist_tmp = new int[mpihelper_.rank()];

      // The index of elmdist_tmp should be the process number, the value the number of elements on each process
      MPI_Comm comm = Dune::MPIHelper::getCommunicator();
      Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

      collective_comm.allgather(&elementNumber, 1, elmdist_tmp);

      // elmdist should be an incremental array whose entries are the sum of all element numbers on previous processes
      elmdist.push_back(0);
      for (int i = 0; i < elementNumber; i++)  { elmdist.push_back(elmdist[i] + elmdist_tmp[i]); }
      delete[] elmdist_tmp;

      // ****************************************************
      // Construct element weights
      // The amount of computation associated with a curvilinear element is approx. calculated:
      //  1) The number of Lagrange Polynomials interpolating the element is equal to the number of interpolation points
      //  2) The number of basis functions to interpolate the field inside should be approximately that number too (why???)
      //  3) The number of new non-zero matrix elements is approx. number of basis functions squared
      // ****************************************************
      std::vector<idx_t> elmwgt;
      for (size_t i = 0; i < elementNumber; i++) { elmwgt.push_back(pow(elements_[i].interpOrder, 2)); }

      // ****************************************************
      // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
      // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
      // ****************************************************
      std::vector<idx_t> eptr, eind;
      int numVertices = 0;
      eptr.push_back(numVertices);

      for (size_t i = 0; i < elementNumber; i++)
      {
    	  std::vector< CornerLocalIndex > cornerIndex = getCornerInfo(elements_[i].geometryType, elements_[i].interpOrder, elements_[i].vertexIndex);

    	  int curNumCorners = cornerIndex.size();
    	  numVertices += curNumCorners;
    	  eptr.push_back(numVertices);

        for (size_t k = 0; k < curNumCorners; ++k)  { eind.push_back(cornerIndex[k]); }
      }

#if PARMETIS_MAJOR_VERSION >= 4
        const int OK =
#endif
        ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), elmwgt.data(), &wgtflag, &numflag,
                                 &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                                 options, &edgecut, reinterpret_cast<idx_t*>(part.data()), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
        if (OK != METIS_OK)
          DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif

    }


  };

} // End of Namespace Dune


#endif // #ifndef DUNE_CURVGRID_GRIDFACTORY_HH
