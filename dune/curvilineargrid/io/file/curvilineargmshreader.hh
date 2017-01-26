/*******************************************************************
 * Curvilinear GMSH Reader
 * 
 * author: Aleksejs Fomins
 * date: 01.08.2014 - created
 * 
 * description:
 * Reads curvilinear GMSH files in parallel. Partitiones the mesh in parallel using ParMetis.
 * - Writes mesh to a factory
 * - [optional] Writes the output to VTK file using Curvilinear VTK Writer
 * 
 *
 * Notifications:
 *   - GMSH provides GlobalIndex for all vertices, however, this index starts at 1. It is therefore lowered to 0 for all vertices and dependent entities
 *
 *
 *******************************************************************/


#ifndef DUNE_CURVILINEARGMSHREADER_HH
#define DUNE_CURVILINEARGMSHREADER_HH

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/loggingtimer.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/io/file/curvilinearvtkwriter.hh>
#include <dune/curvilineargrid/io/file/gmsh2dunemapper.hh>

#ifdef __cplusplus
extern "C" {
#endif
#include <parmetis.h>
#ifdef __cplusplus
}
#endif




namespace Dune
{

namespace CurvGrid {


  enum CurvilinearGmshReaderLoadBalanceStrategy {
	  LoadBalanceDefault,			// All elements of the grid have equal weight. ParMetis splits elements equally, then minimizes process boundary count
	  LoadBalanceBoundary		// Elements on domain boundary have higher weight.
  };


  //! dimension independent parts for CurvilinearGmshReaderParser
  // [TODO] Create new parameter - partition strategy
  //				- VolumeMethodUniformPartition - calculate weight only based on elements
  //				- HybridMethodUniformPartition - calculate reasonably large weight for boundary segments, and small for volume elements
  template<typename GridType, typename FactoryType>
  class CurvilinearGmshReaderParser
  {
  protected:

	typedef typename GridType::ctype  ctype;

    // static data
    static const bool isCached = GridType::is_cached;
	static const int dim_      = GridType::dimension;
    static const int dimWorld_ = GridType::dimensionworld;
    static_assert( (dimWorld_ <= 3), "GmshReader requires dimWorld <= 3." );

    // [TODO] Extend this notation when enabling 2D geometries
    static const int ELEMENT_CODIM = 0;
    static const int FACE_CODIM = 1;
    static const int EDGE_CODIM = 2;
    static const int VERTEX_CODIM = 3;


    // Stores all info associated with an element, except explicit vertex coordinates
    struct GmshEntityData
    {
        int entityIndex_;		// The global index of the entity assigned by GMSH
        int gmshType_;		// The GMSH type denoting the GMSH geometry type of the entity
        int physicalTag_;		// This tag is the same for multiple entities, denoting that they belong to the same large object (e.g. material)
        int elementTag_;		// Another tag always provided by GMSH. Not sure what it means. It is read from the file, but never used
        bool isOnDB_;			// If the entity has a face-contact with the domain boundary
        std::vector<int> processTagSet_;		// GMSH can partition the mesh itself. Then, there are process tags associated with each entity. Currently not used
        std::vector<int> vertexIndexSet_;		// Indices of all vertices that construct this entity. Order of vertices very important for orientation and curvilinear interpolation
    };

    // Logging Message Typedefs
    typedef Dune::CurvGrid::LoggingTimer<LoggingMessage>  LoggingTimer;

    // typedefs
    typedef Dune::FieldVector< ctype, dimWorld_ > GlobalCoordinate;
    typedef Dune::ReferenceElement< ctype, dim_ > ReferenceElement;
    typedef Dune::ReferenceElements< ctype, dim_ > ReferenceElements;

    typedef Dune::ReferenceElement< ctype, dim_-1 > SubReferenceElement;
    typedef Dune::ReferenceElements< ctype, dim_-1 > SubReferenceElements;

    typedef int LocalIndex;
    typedef int GlobalIndex;
    typedef std::set<GlobalIndex>  GlobalIndexSet;
    typedef std::vector<GlobalIndex>  GlobalIndexVector;
    typedef std::vector<LocalIndex>  LocalIndexVector;
    typedef std::vector<GlobalCoordinate>  GlobalCoordinateVector;

    typedef std::map<GlobalIndex, LocalIndex> Global2LocalIndexMap;
    typedef std::map<int, LocalIndex> Tag2LocalIndexMap;
    typedef std::map<GlobalIndexVector, LocalIndexVector> Global2LocalKeyMap;



  public:

    CurvilinearGmshReaderParser(
    	FactoryType & _factory,
    	MPIHelper &mpihelper,
    	bool insertBoundarySegment,
    	bool useGmshElementIndex,
    	bool writeVtkFile,
    	bool partitionMesh,
		CurvilinearGmshReaderLoadBalanceStrategy partStrat
    	) :
    		factory(_factory),
    		mpihelper_ (mpihelper),
    		writeVtkFile_(writeVtkFile),
    		insertBoundarySegment_(insertBoundarySegment),
    		useGmshElementIndex_(useGmshElementIndex),
    		partitionMesh_(partitionMesh),
			partStrat_(partStrat),
    		vtkCurvWriter_(mpihelper),
    		gmsh2dunemapper_()
    {
        // Initialize process parameters
        rank_=mpihelper.rank();
        size_=mpihelper.size();
   }

    int totalVertex()   { return nVertexTotal_; }

    int totalEntity()   { return nEntityTotal_; }

    int totalInteriorElement()   { return nInteriorElementTotal_; }

    int totalBoundarySegment()   { return nBoundarySegmentTotal_; }

    int totalDomanBoundarySegment()   { return nDomainBoundarySegmentTotal_; }

    int totalInteriorBoundarySegment()   { return nInteriorBoundarySegmentTotal_; }


    // This reads the GMSH format to parse the node and element structure
    void read (const std::string& f)
    {
    	fileName = f;

    	std::string log_string;
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: using file " + fileName);
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: reading " + std::to_string(dim_) + "d curvilinear gmsh grid...");

        // open file name, we use C I/O
        // ***********************************************
        FILE* file = fopen(fileName.c_str(),"r");
        if (file==0)  { DUNE_THROW(Dune::IOError, "Could not open " + fileName); }

        // Reading MeshFormat Header
        // ***********************************************
        double version_number;        // process header
        int file_type, data_size;

        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$MeshFormat")!=0)   { DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line"); }
        fscanf(file, "%lg %d %d ", &version_number, &file_type, &data_size);
        if( (version_number < 2.0) || (version_number > 2.3) )  { DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files"); }
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, ":: version " + std::to_string(version_number) + " Gmsh file detected");
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndMeshFormat")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndMeshFormat"); }

        // Reading Node data
        // ***********************************************
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Reading vertex header ------------------------");
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$Nodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $Nodes"); }
        fscanf(file, "%d\n", &nVertexTotal_);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, ":: file contains " + std::to_string(nVertexTotal_) + " vertices");


        //==========================================================
        // VERTEX PASS 1: Put file pointer and skip all vertices
        //==========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Vertex-Pass 1: skip all vertices, since we need element data first ---");
        LoggingTimer::time("CurvilinearGMSHReader: Vertex Pass 1");

        long section_vertex_offset = ftell(file);
        for (int iVertex = 0; iVertex < nVertexTotal_; iVertex++ )
        {
        	fgets(buf_, 512, file );
        	//LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, std::string(buf_));
        }
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndNodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndNodes"); }

        LoggingTimer::time("CurvilinearGMSHReader: Vertex Pass 1");
        //==========================================================
        // VERTEX PASS 1: Finished
        //==========================================================


        // Reading Element Data
        // *************************************************
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Reading elements-header -----------------------");
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$Elements")!=0)  { DUNE_THROW(Dune::IOError, "expected $Elements"); }
        fscanf(file, "%d\n", &nEntityTotal_);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, ":: file contains " + std::to_string(nEntityTotal_) + " elements");



        //=========================================================
        // ELEMENT PASS 1: Count the number of boundary segments
        //=========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Elements-Pass 1: counting elements on the boundary---");
        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 1 - Counting Elements");

        long fileOffsetElementSection = ftell(file);
        for (int i = 0; i < nEntityTotal_; i++)
        {
            int gmshIndex, gmshEntityType;
            fscanf(file, "%d %d ", &gmshIndex, &gmshEntityType);
            int entityIndex = gmsh2dunemapper_.gmsh2DuneIndex(gmshIndex);

            // Ignore the rest of data on this line
            fgets(buf_, 512, file );

            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, std::to_string(entityIndex) + " " + std::to_string(gmshEntityType) );

            // A boundary segment is defined here as any element with dimension less than world dimension
            GeometryType entityType = gmsh2dunemapper_.geometryType(gmshEntityType);

            // Ignore all entities of the types we can not currently handle
            if (gmsh2dunemapper_.canHandleEntityType(gmshEntityType)) {
                if (entityType.dim() < dimWorld_ )     { nBoundarySegmentTotal_++; }
                else                          { nInteriorElementTotal_++; }
            }
        }

        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 1 - Counting Elements");


        //=========================================================
        // ELEMENT PASS 2: Read all linear internal elements
        // Immediately partition them and distribute among processes
        //=========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Elements-Pass 2: reading internal elements for this process ---");
        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 2 - Reading and Partitioning Linear Elements");

        std::set<int> domainBoundaryTagIndexSet;
        GlobalIndexSet thisProcessInteriorElementIndexSet;
        Tag2LocalIndexMap boundaryGlobalTag2IndexMap;

        {
        	std::vector<int> boundaryTagVector;
        	std::vector<GmshEntityData> baseElementVector;

        	fseek(file, fileOffsetElementSection, SEEK_SET);
        	readBaseElements(file, nEntityTotal_, nInteriorElementTotal_, boundaryTagVector, baseElementVector);

        	fseek(file, fileOffsetElementSection, SEEK_SET);
        	markElementsNextToBoundary(file, nEntityTotal_, baseElementVector, boundaryTagVector, boundaryGlobalTag2IndexMap, domainBoundaryTagIndexSet);
        	partitionBaseElements(baseElementVector, thisProcessInteriorElementIndexSet);
        }

        int nBoundaryTag = boundaryGlobalTag2IndexMap.size();

        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 2 - Reading and Partitioning Linear Elements");

        //==========================================================
        // ELEMENT PASS 3: Read all data associated with element.
        // Note all necessary vertex global indices into a map
        // Map all d-1 subentities of all elements to the element localID's.
        //    - Needed to find boundaries corresponding to each element.
        //==========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Elements-Pass 2: reading internal elements for this process ---");
        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 3 - Reading Element Data");

        fseek(file, fileOffsetElementSection, SEEK_SET);
        GlobalIndexSet thisProcessVertexIndexSet;
        std::vector<GmshEntityData> internalElementVector;
        Global2LocalKeyMap boundaryKey2LinkedElementLocalIndexMap;
        readInteriorElements(file, nEntityTotal_, nInteriorElementTotal_,
        		thisProcessInteriorElementIndexSet,
				internalElementVector,
				thisProcessVertexIndexSet,
				boundaryKey2LinkedElementLocalIndexMap);
        int nInternalElement = internalElementVector.size();

        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 3 - Reading Element Data");


        //==========================================================
        // ELEMENT PASS 4: Read all data associated with boundary elements.
        // Only read boundary element if it is associated with this process
        //    - Note: Can not read internal and boundary elements at the same time
        //            because need d-1 subentity map from all internal elements first
        //==========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Elements-Pass 3: reading boundary elements for this process ---");
        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 4 - Reading Boundary Segment Data");

        fseek(file, fileOffsetElementSection, SEEK_SET);
        std::vector<std::vector<GmshEntityData> > tagIndex2boundarySegmentVector(nBoundaryTag);
        std::vector<std::vector<GlobalIndexVector> > tagIndex2linkedElementLocalIndexVectorVector(nBoundaryTag);
        readBoundaryElements(file, nEntityTotal_,
        		boundaryGlobalTag2IndexMap,
				boundaryKey2LinkedElementLocalIndexMap,
				tagIndex2boundarySegmentVector,
				tagIndex2linkedElementLocalIndexVectorVector);

        // If there is no domain boundary on this process, we might not know its index
        int nDomainBoundaryElement = 0.0;
        for (auto && tagIndex : domainBoundaryTagIndexSet) {
        	nDomainBoundaryElement += tagIndex2boundarySegmentVector[tagIndex].size();
        }

        int nBoundaryElement = 0;
        for (int i = 0; i < nBoundaryTag; i++)  { nBoundaryElement += tagIndex2boundarySegmentVector[i].size(); }
        int nInteriorBoundaryElement = nBoundaryElement - nDomainBoundaryElement;

        LoggingTimer::time("CurvilinearGMSHReader: Element Pass 4 - Reading Boundary Segment Data");


        //==========================================================
        // VERTEX PASS 2: Read the vertices
        // But only the ones that correspond to elements on this process
        //==========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Vertex-Pass 2: reading all vertices necessary for this process ---");
        LoggingTimer::time("CurvilinearGMSHReader: Vertex Pass 4 - Reading Associated Vertices");

        fseek(file, section_vertex_offset, SEEK_SET);
        std::map<GlobalIndex, GlobalCoordinate> vertexIndex2CoordinateMap;                // Only for testing purposes
        Global2LocalIndexMap vertexGlobal2LocalIndexMap;
        int nVertex = readVertices(file, nVertexTotal_, vertexIndex2CoordinateMap, thisProcessVertexIndexSet, vertexGlobal2LocalIndexMap);

        LoggingTimer::time("CurvilinearGMSHReader: Vertex Pass 4 - Reading Associated Vertices");

        //==========================================================
        // Final Step: Insert boundary segments and elements
        //==========================================================
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "----- Adding internal boundary elements to factory ---");
        LoggingTimer::time("CurvilinearGMSHReader: Inserting Entities into the factory");

        addInternalElements(vertexGlobal2LocalIndexMap, vertexIndex2CoordinateMap, internalElementVector, nBoundarySegmentTotal_);
        for (int iTag = 0; iTag < nBoundaryTag; iTag++) {
        	if (tagIndex2boundarySegmentVector[iTag].size() > 0) {
        		bool isDomainBoundary =  (domainBoundaryTagIndexSet.find(iTag) != domainBoundaryTagIndexSet.end());
        		addBoundaryElements(
        				vertexGlobal2LocalIndexMap,
						vertexIndex2CoordinateMap,
						tagIndex2boundarySegmentVector[iTag],
						tagIndex2linkedElementLocalIndexVectorVector[iTag],
						isDomainBoundary
        		);
        	}
        }

        Dune::CollectiveCommunication<MPI_Comm> comm = mpihelper_.getCollectiveCommunication();
        int nElementParallelSum  = comm.sum(nInternalElement);
        nDomainBoundarySegmentTotal_ = comm.sum(nDomainBoundaryElement);
        nInteriorBoundarySegmentTotal_ = comm.sum(nInteriorBoundaryElement);

        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "Finished Reading Mesh");
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: total vertices          = " + std::to_string(nVertexTotal_)          + " of which on this process " + std::to_string(nVertex) );
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: total internal elements = " + std::to_string(nInteriorElementTotal_) + " of which on this process " + std::to_string(nInternalElement) );
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: total domain boundary segments = " + std::to_string(nDomainBoundarySegmentTotal_) + " of which on this process " + std::to_string(nDomainBoundaryElement) );
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: total interior boundary segments = " + std::to_string(nInteriorBoundarySegmentTotal_) + " of which on this process " + std::to_string(nInteriorBoundaryElement) );
        LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, ":: total boundary segments = " + std::to_string(nBoundarySegmentTotal_) + " of which on this process " + std::to_string(nBoundaryElement) );


        assert(nBoundaryElement == nDomainBoundaryElement + nInteriorBoundaryElement);
        // Note: In case of Interior boundaries the below equality does not hold - a signle interior boundary can be shared by 1 or 2 processes
        //assert(nBoundarySegmentTotal_ == nDomainBoundarySegmentTotal_ + nInteriorBoundarySegmentTotal_);


        if ((nElementParallelSum  != nInteriorElementTotal_) && (rank_ == 0))
        {
        	std::cout << "The initially-calculated nInternal=" << nInteriorElementTotal_ << " does not match the inserted one = " << nElementParallelSum << std::endl;
        	DUNE_THROW(IOError, "Wrong number of internal elements");
        }
        /*
        if ((nBoundaryParallelSum  != nDomainBoundarySegmentTotal_) && (comm.rank() == 0))
        {
        	std::cout << "The initially-calculated nBoundary=" << nDomainBoundarySegmentTotal_ << " does not match the inserted one = " << nBoundaryParallelSum << std::endl;
        	DUNE_THROW(IOError, "Wrong number of boundary elements");
        }*/


        // TESTING SECTION - WRITES TEST ELEMENTS TO .VTK FILE
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (writeVtkFile_)
        {
        	LoggingTimer::time("CurvilinearGMSHReader: Writing VTK output");
        	//vtkCurvWriter_.writeVTK("./curvreader_output_process_" + std::to_string(rank_) + ".vtk");
        	vtkCurvWriter_.writeParallelVTU("./", "curvreader_output", VTU_DATA_FORMAT_ASCII);
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__,  "Curvilinear VTK Writer finished writing" );
            LoggingTimer::time("CurvilinearGMSHReader: Writing VTK output");
        }

        LoggingTimer::time("CurvilinearGMSHReader: Inserting Entities into the factory");

        // Close file
        fclose(file);
      }





  protected:

    // ***********************************************************************
    // Auxiliary Methods
    // ***********************************************************************

    // Testing if the current element type can be handled by DUNE
    // [TODO] Move to Gmsh2DuneMapper.hh
    bool checkEntityAllowed(int gmshIndex)
    {
    	GeometryType gt = gmsh2dunemapper_.geometryType(gmshIndex);
        bool isAllowedEntity = true;

        // Only allow simplex geometries for now
        isAllowedEntity &= gt.isSimplex();

        // Check if element is polynomial-complete (ask GMSH what that even means I dont know :) )
        isAllowedEntity &= !gmsh2dunemapper_.hasIncompleteOrder(gmshIndex);

        // test whether we support the element type at the moment
        if (!isAllowedEntity) { DUNE_THROW(Dune::IOError, "GMSH Reader: Have read an element of unexpected type "); }

        return isAllowedEntity;
    }


    // Checks whether this element (or boundary element) belongs on this parallel process
    bool elementOnProcess(int eIndex, int eTotal) {
        int eFirst = (eTotal * rank_) / size_;
        int eLast = (eTotal * (rank_+1)) / size_;

        std::string log_string = " == checkprocess if " + std::to_string(eIndex) + " in [" + std::to_string(eFirst) + "," + std::to_string(eLast) + "]";
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

        return ((eIndex >= eFirst)&&(eIndex < eLast));
    }


    // Reads continuous set of vertices from a file. Must be at the correct position in file already
    void readEntityVertices(FILE* file, GmshEntityData & entityData, int nVertex) {
        for (int iVertex = 0; iVertex < nVertex; iVertex++) {
            int gmshIndex;
            fscanf(file, "%d", &gmshIndex);
            int vertexGlobalIndex = gmsh2dunemapper_.gmsh2DuneIndex(gmshIndex);
            entityData.vertexIndexSet_.push_back(vertexGlobalIndex);
        }
        std::string log_string = "  --- read entity vertices = " + VectorHelper::vector2string(entityData.vertexIndexSet_);
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
    }




    // ***********************************************************************
    // Reader Sub-procedures
    // ***********************************************************************


    /** \brief Reads vertices into a map, given a set of indices of vertices that should be read
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  nVertexTotal_                 the total number of vertices specified in the file
     *  \param[in]  vertices                       the map from vertex globalID to vertex coordinate. Here the vertices will be stored
     *  \param[in]  vertexIndexSet               the set of globalID's of all vertices that belong to this process
     *  \param[in]  vertexGlobal2LocalIndexMap  the map from vertex globalID to localID
     *
     *  \return The total number of vertices on this process
     *  \note Assumes it is in the correct position in the file
     *
     */
    int readVertices(
            FILE* file,
            int nVertexTotal_,
            std::map<GlobalIndex, GlobalCoordinate> & vertexIndex2CoordinateMap,
            std::set<GlobalIndex> & vertexIndexSet,
            Global2LocalIndexMap & vertexGlobal2LocalIndexMap
            )
    {
        int gmshIndex;
        GlobalCoordinate x;

        // Iterator starts from 1 because GMSH numbers vertices [1,n]
        for( int i = 0; i < nVertexTotal_; ++i )
        {
        	LoggingMessage::writePatience(" Reading vertices...", i, nVertexTotal_);

            // If this vertex does not belong to this process, just skip it
            if (vertexIndexSet.count(i) == 0)  { fgets(buf_, 512, file ); }
            else
            {
                fscanf(file, "%d ", &gmshIndex);
                int vertexIndex = gmsh2dunemapper_.gmsh2DuneIndex(gmshIndex);
                std::string tmp_out;

                if( vertexIndex != i )  { DUNE_THROW( Dune::IOError, "Expected id " << i << ", got id " << vertexIndex << "." ); }
                for (int d = 0; d < dimWorld_; d++) {
                    double tmp_coord;
                    fscanf(file, "%lg", &tmp_coord);
                    x[d] = tmp_coord;

                    tmp_out += std::to_string(tmp_coord) + " ";
                }

                // The local id of this vertex is equalt to the number of vertices added so far
                vertexGlobal2LocalIndexMap[i] = vertexIndex2CoordinateMap.size();

                LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "  * Have read vertex " + tmp_out);

                // Maps global id to global coordinate
                vertexIndex2CoordinateMap[i] = x;

                // Insert vertex into a factory, noting its globalId.
                // Its localId in the factory should be given by the order the vertices are added to the factory
                factory.insertVertex(x, i);

                //fscanf(file, "\n");
                fgets(buf_, 512, file );
                // [TODO] Assert that buffer is empty - there should be nothing left on this line
            }
        }

        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndNodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndNodes"); }

        return vertexIndex2CoordinateMap.size();
    }


    // Reads the data about this entity, everything except the interpolation vertex id's
    // Assumes it is in the correct position in the file
    GmshEntityData readEntitySpec(FILE* file)
    {
        int nTag;
        int gmshIndex;
        GmshEntityData thisEntity;

        fscanf(file, "%d %d %d ", &gmshIndex, &thisEntity.gmshType_, &nTag);
        thisEntity.entityIndex_ = gmsh2dunemapper_.gmsh2DuneIndex(gmshIndex);

        std::stringstream log_string;
        log_string << "    * element " << thisEntity.entityIndex_ << " has " << nTag << " tags";
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string.str());

        /** \brief Reading tags
         *
         * Possible tags
         *   k == 0: physical entity tag (for example, material properties of the element)
         *   k == 1: elementary entity (not used here)
         *   if version_number < 2.2:
         *     k == 2: mesh partition 0
         *   else
         *     k == 3: number of mesh partitions
         *     k => 3: mesh partition k-3
         *
         * **/
        fscanf(file, "%d ", &thisEntity.physicalTag_);
        fscanf(file, "%d ", &thisEntity.elementTag_);

        // Possible functionality for mesh partitioning using GMSH
        // TODO: Currently this functionality not available
        for (int k = 2; k < nTag; k++)
        {
            int tmp_tag;
            fscanf(file, "%d ", &tmp_tag);
            thisEntity.processTagSet_.push_back(tmp_tag);
        }

        // Note: Initially mark all entities as non-boundary
        thisEntity.isOnDB_ = false;

        return thisEntity;
    }


    /** \brief Communicate the mesh as partitioned by partitionCompute
     * First, communicates how many elements are sent to each process from this process.
     * Then, communicates all elements (only globalId's)
     *
     * \param[in]  file                           file pointer to read from
     * \param[in]  nEntityTotal_                 the total number of elements specified in the file
     * \param[in]  nInteriorElementTotal_        the number of internal elements specified in the file
     * \param[in]  thisProcessInteriorElementIndexSet           set of globalId's of all elements present on this process
     *
     */
    void readBaseElements(
            FILE* file,
            int nEntityTotal_,
            int nInteriorElementTotal_,
			std::vector<int> & boundaryTagVector,
			std::vector<GmshEntityData> & baseElementVector
			)
    {
    	// Count the number of elements currently stored on this process
        int iSelectElem = 0;

        // Record the different boundary segment tags
        std::set<int> boundaryTagSet;

        // Reading element info - tag information and vertex global indices
        // *************************************************************
        for (int i = 0; i < nEntityTotal_; i++)
        {
        	LoggingMessage::writePatience(" Reading linear elements for partitioning...", i, nEntityTotal_);

            // Read the first part of the element info
            GmshEntityData thisElement = readEntitySpec(file);


            if (gmsh2dunemapper_.canHandleEntityType(thisElement.gmshType_)) {
            	// Find if this is a boundary element
				GeometryType elemType = gmsh2dunemapper_.geometryType(thisElement.gmshType_);

				// If this is a boundary segment, we should note its tag for future use
				if (elemType.dim() < dimWorld_) {
					boundaryTagSet.insert(thisElement.physicalTag_);

				// If it is an interior element that belongs to this process, it should be read and used for partitioning
				} else if (elementOnProcess(iSelectElem++, nInteriorElementTotal_)) {
					// Testing if the current element type can be handled by DUNE
					// *****************************************************************
					checkEntityAllowed(thisElement.gmshType_);
					std::string log_string = "    * element " + std::to_string(thisElement.entityIndex_) + " can be treated by Dune grid ";
					LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

					// Obtain all necessary info not to use gmshElementIndex in the following steps
					// *****************************************************
					int thisElmOrder               = gmsh2dunemapper_.elementOrder(thisElement.gmshType_);
					int thisElmDofNo               = gmsh2dunemapper_.dofNumber(thisElement.gmshType_);
					int thisElmCorners             = ReferenceElements::general(elemType).size(elemType.dim());
					int thisElmSubentities         = ReferenceElements::general(elemType).size(1);

					// Reading Corners of the Element
					// Note: In GMSH notation corners go first
					// *************************************************
					readEntityVertices(file, thisElement, thisElmCorners);

					// Store elements just read
					baseElementVector.push_back(thisElement);
				}
            }

            // Read until the end of the line
            fgets(buf_, 512, file );
        }
        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Finished reading base elements");


        // Convert the boundary tag set to a sorted tag array for future use
    	// Very Important: The tags need to be sorted, so that their index is the same on all processes
        // *************************************************************
    	boundaryTagVector = std::vector<int>(boundaryTagSet.begin(), boundaryTagSet.end());
    	std::sort(boundaryTagVector.begin(), boundaryTagVector.end());
    	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "Found boundary segments with tags: " + VectorHelper::vector2string(boundaryTagVector));


        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    void markElementsNextToBoundary(
            FILE* file,
            int nEntityTotal_,
			std::vector<GmshEntityData> & baseElementVector,
			std::vector<int> & boundaryTagVector,
			Tag2LocalIndexMap & boundaryGlobalTag2IndexMap,
			std::set<int> & domainBoundaryTagIndexSet)
    {
    	typedef std::map<std::vector<int>, int> Key2IndexMap;
    	typedef typename Key2IndexMap::iterator  Key2IndexIter;

    	Dune::GeometryType gtElem(Dune::GeometryType::BasicType::simplex, dim_);
    	int nFacePerElem = ReferenceElements::general(gtElem).size(FACE_CODIM);

    	// To save work, determine which tags are certainly IB to save on global communication later
    	int nPhysicalTag = boundaryTagVector.size();
    	std::vector<int> isTagInteriorBoundary(nPhysicalTag, 0);

    	// Step 1:  Make map facekey->faceindex->baseElementLocalInd
    	// ********************************************************************
    	Key2IndexMap faceKey2IndexMap;
    	std::vector<std::vector<int>> faceIndex2assocElemIndex;
    	for (int iElem = 0; iElem < baseElementVector.size(); iElem++) {
    		for (int iFace = 0; iFace < nFacePerElem; iFace++) {

    			// Extract the corner set of all faces of this element, transform into a key
    			std::vector<int> faceKey;
    			std::vector<int> faceCornerInterIndex = CurvilinearGeometryHelper::template subentityInternalCoordinateSet<ctype, dim_>(gtElem, 1, FACE_CODIM, iFace);
    			for (int iCorner = 0; iCorner < faceCornerInterIndex.size(); iCorner++) {
    				faceKey.push_back(baseElementVector[iElem].vertexIndexSet_[faceCornerInterIndex[iCorner]]);
    			}
    			std::sort(faceKey.begin(), faceKey.end());

    			// Insert all keys into the key map
    			Key2IndexIter it = faceKey2IndexMap.find(faceKey);
    			if (it == faceKey2IndexMap.end()) {
    				int pos = faceIndex2assocElemIndex.size();
    				faceKey2IndexMap[faceKey] = pos;
    				faceIndex2assocElemIndex.push_back(std::vector<int> {iElem});
    			} else {
    				faceIndex2assocElemIndex[it->second].push_back(iElem);
    			}
    		}
    	}

    	// Step 2: Process the boundary tags - make a map from sorted boundary tags to a new index
    	// ********************************************************************
    	for (unsigned int iTag = 0; iTag < boundaryTagVector.size(); iTag++) {
    		boundaryGlobalTag2IndexMap[boundaryTagVector[iTag]] = iTag;
    	}


        // Step 3: Read all boundary segments, subselect only those that correspond to element faces currently on this process
    	// ********************************************************************
        std::vector<Global2LocalIndexMap> tagInd2FaceGlobal2LocalIndexVector(nPhysicalTag);
        for (int i = 0; i < nEntityTotal_; i++)
        {
        	LoggingMessage::writePatience(" Reading boundary segments 1st time...", i, nEntityTotal_);

            // Read the first part of the boundary face info
            GmshEntityData thisElement = readEntitySpec(file);

            if (gmsh2dunemapper_.canHandleEntityType(thisElement.gmshType_)) {
                // Find if this is a boundary element
                GeometryType elemType = gmsh2dunemapper_.geometryType(thisElement.gmshType_);

                // If this element is not on the boundary just skip the rest of info on this line
                if (elemType.dim() == dimWorld_) { fgets(buf_, 512, file ); }
                else
                {
                	/************************************************
                	 *   Read the boundary segment vertices
                	 ************************************************/

                    // Testing if the current element type can be handled by DUNE
                    checkEntityAllowed(thisElement.gmshType_);
                    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__,  "    * element " + std::to_string(i) + " can be treated by Dune grid ");

                    // Obtain all necessary info not to use gmshElementIndex in the following steps
                    int thisElmCorners = SubReferenceElements::general(elemType).size(elemType.dim());

                    // Reading DoF's. Extract only corners. Note that corners come first in GMSH notation, so easy to read
                    readEntityVertices(file, thisElement, thisElmCorners);

                    // Finish reading line
                    fgets(buf_, 512, file );
                    // fscanf(file, "\n");

                	/************************************************
                	 *   Process the boundary - seek it among the faces of the elements currently on this process
                	 ************************************************/
                    std::vector<int> faceKey = thisElement.vertexIndexSet_;
                    std::sort(faceKey.begin(), faceKey.end());

                    // Subselect only the BS that are faces of the elements present on this process
                    // If it is already obvious that a BS is shared by 2 elements, mark it directly as IB
                    Key2IndexIter faceIter = faceKey2IndexMap.find(faceKey);
                    if (faceIter != faceKey2IndexMap.end()) {
                    	typename Tag2LocalIndexMap::iterator tagIter = boundaryGlobalTag2IndexMap.find(thisElement.physicalTag_);
                    	assert(tagIter != boundaryGlobalTag2IndexMap.end());

                    	tagInd2FaceGlobal2LocalIndexVector[tagIter->second][thisElement.entityIndex_] = faceIter->second;
                    	if (faceIndex2assocElemIndex[faceIter->second].size() == 2)  { isTagInteriorBoundary[tagIter->second] = 1; }
                    }
                }
            } else {
                // Finish reading line
                fgets(buf_, 512, file );
            }
        }

        // Step 4. Determine which of the tags is the Domain Boundary tag
        // Note: It is essential that all processes know this tag, even if they do not have any entities of this tag
        // ********************************************************************

        if (size_ > 1) {  // If the code is serial, this information should be obvious without communication
        	for (auto && tagIter : boundaryGlobalTag2IndexMap) {
        		int thisTag = tagIter.first;
        		int thisTagInd = tagIter.second;
        		int thisTagBSCount = tagInd2FaceGlobal2LocalIndexVector[thisTagInd].size();

        		int tagDBIndicator;
        		if (isTagInteriorBoundary[thisTagInd] == 1)	{ tagDBIndicator = -2; } // This process knows for sure that this field is IB
        		else if (thisTagBSCount == 0)						{ tagDBIndicator = -1; } // This process has no boundary segments of this tag, so it knows nothing
        		else																	{ tagDBIndicator = tagInd2FaceGlobal2LocalIndexVector[thisTagInd].begin()->first; }  // Check by comparing any face (e.g. first one) with the rest of processes

   	        	int tagAssoc[size_];

#if HAVE_MPI
    MPI_Comm comm = Dune::MPIHelper::getCommunicator();
    MPI_Allgather (&tagDBIndicator, 1, MPI_INT, tagAssoc, 1, MPI_INT, comm);
#endif

    			/*
					std::cout << "brrrrrrr (" << rank_ << "," << i << "): " << tagDBIndicator;
					for (int j = 0; j < size_; j++) { std::cout << " " << tagAssoc[j]; }
					std::cout << std::endl;
    			*/

    			// Only do sth if we don't know if the face is IB yet
    			// However, if there are no BS of this tag at all on the processor, it is irrelevant
    			if ((isTagInteriorBoundary[thisTagInd] == 0) && (thisTagBSCount > 0)) {
    				// Check if any other process already believes it is IB
    				for (int iProc = 0; iProc < size_; iProc++) { if (tagAssoc[iProc] == -2)  { isTagInteriorBoundary[thisTagInd] = 1; } }

    				// If nobody yet believes it is IB, it could still be that the IB is entirely contained in the interprocessor boundary by a lucky chance
    				// Need to check if the faces of this boundary appear only once over all processes
    				if (isTagInteriorBoundary[thisTagInd] == 0) {
        				int iProc = 0;
        				int tagParallelCount = 0;
    					while ((tagParallelCount < 2) && (iProc < size_)) {
    						if (tagAssoc[iProc] >= 0)  {
    							typename Global2LocalIndexMap::iterator indexIter = tagInd2FaceGlobal2LocalIndexVector[thisTagInd].find(tagAssoc[iProc]);
    							if (indexIter != tagInd2FaceGlobal2LocalIndexVector[thisTagInd].end())  { tagParallelCount++; }
    						}
    						iProc++;
    					}
    					assert(tagParallelCount > 0);  // Should at least find the contribution it had sent itself,if it had actually checked for it
    					if (tagParallelCount == 2) { isTagInteriorBoundary[thisTagInd] = 1; }
    				}
    			}
            }
        }

        // After above iteration, for each IB there must exist at least 2 processes that are certain about it
        // However, no process is guaranteed to know about all IB.
        // Thus it is essential to combine the knowledge of all processes to obtain the full picture
        Dune::CollectiveCommunication<MPI_Comm> comm = mpihelper_.getCollectiveCommunication();
        for (int iTagInd = 0; iTagInd < nPhysicalTag; iTagInd++) {
        	int tagDBVote  = comm.sum(isTagInteriorBoundary[iTagInd]);
        	assert(tagDBVote >= 0);

        	if (tagDBVote == 0) {
        		domainBoundaryTagIndexSet.insert(iTagInd);
        		//domainBoundaryTagVector.push_back(boundaryTagVector[iTagInd]);
        	}
        }

        // Check that there is exactly 1 domain boundary
        // NOTE: IN FACT, WE MUST ALLOW TO HAVE SEVERAL USER-DEFINED DB TAGS
        // assert(domainBoundaryTagIndexSet.size() == 1);


        // Step 5. Mark all elements neighbouring a Domain Boundary that they indeed do
        // ********************************************************************
        for (auto && tagIndex : domainBoundaryTagIndexSet) {
			for (auto && iter : tagInd2FaceGlobal2LocalIndexVector[tagIndex]) {
				std::vector<int> & assocElemInd = faceIndex2assocElemIndex[iter.second];
				assert(assocElemInd.size() == 1); // Domain Boundaries must have only one neighbouring element
				baseElementVector[assocElemInd[0]].isOnDB_ = true;
			}
        }
    }



    void partitionBaseElements(
    		std::vector<GmshEntityData> & baseElementVector,
        	GlobalIndexSet & thisProcessInteriorElementIndexSet)
    {
            if ((size_ > 1) && partitionMesh_)
            {
    #if HAVE_MPI
            	// Partition the elements
            	// Communicate the partitioning
            	// Fill the element set
            	// *************************************************************
            	std::vector<unsigned> part(baseElementVector.size(), 0);
        		partitionCompute(part, baseElementVector);
        		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Finished computing partition");

        		partitionCommunicate(part, baseElementVector, thisProcessInteriorElementIndexSet);
        		LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Finished communicating partition");

            	std::vector<int> test2 (thisProcessInteriorElementIndexSet.begin(), thisProcessInteriorElementIndexSet.end());
            	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " elements after partition: " + VectorHelper::vector2string(test2));
    #endif
            } else
            {
            	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, "No MPI found! Running sequential case without partitioning");
            	for (unsigned int i = 0; i < baseElementVector.size(); i++) { thisProcessInteriorElementIndexSet.insert(baseElementVector[i].entityIndex_); }
            }
    }


    /** \brief Reads all internal element data into a vector. Only reads elements which should be on this process
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  nEntityTotal_                 the total number of elements specified in the file
     *  \param[in]  nInteriorElementTotal_        the number of internal elements specified in the file
     *  \param[in]  thisProcessInteriorElementIndexSet          set of globalId's of all elements present on this process
     *  \param[in]  internalElementVector              A vector in which the globalID's of internal elements will be stored
     *  \param[in]  vertexIndexSet               the set of globalID's of all vertices that belong to this process
     *  \param[in]  boundaryKey2LinkedElementLocalIndexMap        A map from a sorted array of globalID's of vertices that make up a boundary to an array of localID's of internal elements to whom this boundary belongs.
     *
     *  \note Assumes it is in the correct position in the file
     *
     *  TODO: Make factory.insertVertex() work with globalId
     */
    void readInteriorElements(
            FILE* file,
            int nEntityTotal_,
            int nInteriorElementTotal_,
            GlobalIndexSet & thisProcessInteriorElementIndexSet,
            std::vector<GmshEntityData> & internalElementVector,
			GlobalIndexSet & vertexIndexSet,
            Global2LocalKeyMap & boundaryKey2LinkedElementLocalIndexMap)
    {
    	std::string log_string;

        // Reading element info - tag information and vertex global indices
        for (int i = 0; i < nEntityTotal_; i++)
        {
        	LoggingMessage::writePatience(" Reading internal elements...", i, nEntityTotal_);

            // Read the first part of the element info
            GmshEntityData thisElement = readEntitySpec(file);

            if (gmsh2dunemapper_.canHandleEntityType(thisElement.gmshType_)) {
            	// Find if this is a boundary element
				GeometryType elemType          = gmsh2dunemapper_.geometryType(thisElement.gmshType_);
				//bool onBoundary = (elemType.dim() < dimWorld_);

				// Check if this element belongs on this process
				if ((elemType.dim() < dimWorld_)||(thisProcessInteriorElementIndexSet.count(thisElement.entityIndex_) == 0)) { fgets(buf_, 512, file ); }
				else
				{
					// Testing if the current element type can be handled by DUNE
					// *****************************************************************
					checkEntityAllowed(thisElement.gmshType_);
					log_string = "    * element " + std::to_string(thisElement.entityIndex_) + " can be treated by Dune grid ";
					LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

					// Obtain all necessary info not to use gmshElementIndex in the following steps
					// *****************************************************
					int thisElmOrder               = gmsh2dunemapper_.elementOrder(thisElement.gmshType_);
					int thisElmDofNo               = gmsh2dunemapper_.dofNumber(thisElement.gmshType_);
					int thisElmCorners             = ReferenceElements::general(elemType).size(elemType.dim());
					int thisElmSubentities         = ReferenceElements::general(elemType).size(1);
					int thisElmSubCorners          = ReferenceElements::general(elemType).size(0, 1, elemType.dim());


					// Reading DoF's
					// *************************************************
					readEntityVertices(file, thisElement, thisElmDofNo);

					// Store all vertices of this element on this process
					// Note: set ignores request to add an already existing vertex
					for (int iDof = 0; iDof < thisElmDofNo; iDof++) { vertexIndexSet.insert(thisElement.vertexIndexSet_[iDof]); }

					// Finish reading line
					fgets(buf_, 512, file );
					// [TODO] assert that the buffer is empty - there should be nothing left in this string
					//fscanf(file, "\n");


					// Add all subentities of this element to a map, such that it is easy afterwards to find
					// boundary elements associated with this process
					// *************************************************************************************

					// 1) get all corners of this process
					// Note that in GMSH convention the corners come first in the Dof array
					GlobalIndexVector cornerVector;
					for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) { cornerVector.push_back(thisElement.vertexIndexSet_[iCorner]); }

					// 2) sort by increasing globalID - Sorting really unnecessary since it is only used to produce the key,
					//    and the key has to be sorted regardless
					//std::sort(cornerVector.begin(), cornerVector.end());

					// 3) Calculate the localID of this element
					int localID = internalElementVector.size();

					for (int iSub = 0; iSub < thisElmSubentities; iSub++)
					{
						log_string = " for iSub = " + std::to_string(iSub) + " out of " + std::to_string(thisElmSubentities) + " have sub_size = " + std::to_string(thisElmSubCorners);
						LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

						// 4) get all subsets associated with this element type
						//std::vector<int> thisSubentityCornerIndexSet = CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet;
						GlobalIndexVector key;

						for (int iCoord = 0; iCoord < thisElmSubCorners; iCoord++) {
							int thisCornerInternalSubIndex = ReferenceElements::general(elemType).subEntity(iSub, 1, iCoord, elemType.dim());
							key.push_back(cornerVector[thisCornerInternalSubIndex]);
						}

						// Sort the key - It must have the same orientation for both neighboring elements
						std::sort(key.begin(), key.end());


						// 5) Check if map empty for this entry, then add to the map
						if (boundaryKey2LinkedElementLocalIndexMap.find(key) == boundaryKey2LinkedElementLocalIndexMap.end()) {
							log_string = " -- element " + std::to_string(localID) + " has added boundary " + VectorHelper::vector2string(key);
							LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

							boundaryKey2LinkedElementLocalIndexMap[key] = LocalIndexVector (1, localID);
						} else
						{
							// 6) In this case this is the 2nd element sharing this boundary
							LocalIndexVector tmp = boundaryKey2LinkedElementLocalIndexMap[key];
							tmp.push_back(localID);
							boundaryKey2LinkedElementLocalIndexMap[key] = tmp;   // Should overwrite previous value

							log_string = " -- element " + std::to_string(localID) + " shares boundary " + VectorHelper::vector2string(key) + " with element " + std::to_string(tmp[0]);
							LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
						}
					}

					// correct differences between gmsh and Dune in the local vertex numbering
					// *************************************************
					gmsh2dunemapper_.gmsh2DuneElementDofNumbering(elemType, thisElmOrder, thisElement.vertexIndexSet_);

					internalElementVector.push_back(thisElement);
				}
            } else {
                // Finish reading line
                fgets(buf_, 512, file );
            }
        }

        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    /** \brief Reads all internal element data into a vector. Only reads elements which should be on this process
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  nEntityTotal_                 the total number of elements specified in the file
     *  \param[in]  tagIndex2boundarySegmentVector              A vector that stores global indices of boundary segments as function of a boundary subdomain physical tag
     *  \param[in]  boundaryKey2LinkedElementLocalIndexMap        A map from a sorted array of globalID's of vertices that make up a boundary to an array of localID's of internal elements to whom this boundary belongs.
     *  \param[in]  linkedElementLocalIndexVector   A vector that stores a vector of localID's of all elements linked this boundary, for each boundary localID
     *  \param[in]  domainBoundaryTagIndexSet   The tag index associated with domain boundary, to distinguish from other interior boundaries
     *
     *  \note Assumes it is in the correct position in the file
     */
    void readBoundaryElements(
            FILE* file,
            int nEntityTotal_,
			Tag2LocalIndexMap & boundaryGlobalTag2IndexMap,
			Global2LocalKeyMap & boundaryKey2LinkedElementLocalIndexMap,
			std::vector<std::vector<GmshEntityData> > & tagIndex2boundarySegmentVector,
            std::vector<std::vector<LocalIndexVector> > & tagIndex2linkedElementLocalIndexVectorVector)
    {
    	std::string log_string;

        int iSelectElem = 0;

        int nPhysicalTag = boundaryGlobalTag2IndexMap.size();
        assert(tagIndex2boundarySegmentVector.size() == nPhysicalTag);
        assert(tagIndex2linkedElementLocalIndexVectorVector.size() == nPhysicalTag);

        std::vector<bool> isTagInteriorBoundary(nPhysicalTag, false);

        // Reading element info - tag information and vertex global indices
        // ********************************************************************
        for (int i = 0; i < nEntityTotal_; i++)
        {
        	LoggingMessage::writePatience(" Reading boundary segments...", i, nEntityTotal_);

            // Read the first part of the element info
            GmshEntityData thisElement = readEntitySpec(file);

            if (gmsh2dunemapper_.canHandleEntityType(thisElement.gmshType_)) {
                // Find if this is a boundary element
                GeometryType elemType          = gmsh2dunemapper_.geometryType(thisElement.gmshType_);

                // If this element is not on the boundary just skip the rest of info on this line
                if (elemType.dim() == dimWorld_) { fgets(buf_, 512, file ); }
                else
                {
                    // Testing if the current element type can be handled by DUNE
                    // *****************************************************************
                    checkEntityAllowed(thisElement.gmshType_);
                    log_string = "    * element " + std::to_string(i) + " can be treated by Dune grid ";
                    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

                    // Obtain all necessary info not to use gmshElementIndex in the following steps
                    // *****************************************************
                    int thisElmOrder               = gmsh2dunemapper_.elementOrder(thisElement.gmshType_);
                    int thisElmDofNo               = gmsh2dunemapper_.dofNumber(thisElement.gmshType_);
                    int thisElmCorners             = SubReferenceElements::general(elemType).size(elemType.dim());

                    // Reading DoF's
                    // *************************************************
                    readEntityVertices(file, thisElement, thisElmDofNo);

                    // Finish reading line
                    fgets(buf_, 512, file );
                    // [TODO] Assert that buf is empty - there should be nothing left on this line
                    //fscanf(file, "\n");



                    // Finding the associated index in the tag array
                    // *****************************************************
                    typename Tag2LocalIndexMap::iterator tagIter = boundaryGlobalTag2IndexMap.find(thisElement.physicalTag_);
                    assert(tagIter != boundaryGlobalTag2IndexMap.end());  // Should not find non-existing tags
                    int tagIndex = tagIter->second;


                    // Check if associated with any internal element
                    // ****************************************************************

                    // 1) Obtain corners
                    GlobalIndexVector cornerVector;
                    for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) { cornerVector.push_back(thisElement.vertexIndexSet_[iCorner]); }

                    // 2) sort by increasing globalID
                    std::sort(cornerVector.begin(), cornerVector.end());

                    // 3) Obtain this boundary element's localID
                    //int localID = tagIndex2boundarySegmentVector[tagIndex].size();

                    log_string = " -- boundary " + std::to_string(thisElement.entityIndex_) + " checking key " + VectorHelper::vector2string(cornerVector);
                    LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);



                    // If this boundary already exists in the map, then add this element
                    // If not, it is likely on another process
                    if (boundaryKey2LinkedElementLocalIndexMap.find(cornerVector) != boundaryKey2LinkedElementLocalIndexMap.end())
                    {
                    	log_string = " -found b.e localID = " + VectorHelper::vector2string(boundaryKey2LinkedElementLocalIndexMap[cornerVector]);
                        LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);

                        // correct differences between gmsh and Dune in the local vertex numbering
                        // *************************************************
                        gmsh2dunemapper_.gmsh2DuneElementDofNumbering(elemType, thisElmOrder, thisElement.vertexIndexSet_);

                        // Add boundary element
                        tagIndex2boundarySegmentVector[tagIndex].push_back(thisElement);

                        // If this boundary is linked to an element, write its localID to the internal_element_boundaries for corresponding element
                        // DUNE_THROW(Dune::IOError, "Have " + std::to_string(thisBoundaryLinkedElements.size()) + " elements associated with 1 boundary. Not implemented yet");
                        LocalIndexVector thisBoundaryLinkedElements = boundaryKey2LinkedElementLocalIndexMap[cornerVector];
                        if (thisBoundaryLinkedElements.size() > 1)  { isTagInteriorBoundary[tagIndex] = true; }

                        tagIndex2linkedElementLocalIndexVectorVector[tagIndex].push_back(thisBoundaryLinkedElements);
                    }
                }
            } else {
                // Finish reading line
                fgets(buf_, 512, file );
            }
        }

        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    /** \brief Adds all internal elements to factory, also writes them to .vtk file.
     *
     *  \param[in]  vertexGlobal2LocalIndexMap    A map from global vertex index to vertex index local to this process
     *  \param[in]  vertexIndex2CoordinateMap     the map from vertex globalID to vertex coordinate
     *  \param[in]  internalElementVector         A vector in which the globalID's of internal elements are stored
     *  \param[in]  nBoundarySegmentTotal         Number of boundary segments over all processes
     *
     *  [TODO] Currently, physical tag is inserted explicitly, which may be unsatisfactory for Dune community.
     *          A possible workaround is to create compiler directive -DHAVE_PHYSICAL_TAG, which would determine
     *          whether to use extra argument in the factory.insertElement() routine
     *
     *  [TODO] The vertex and element vectors are stored twice - once inside the read procedure and once in factory
     *  as they are being added. Maybe possible to save space
     */
    void addInternalElements(
            Global2LocalIndexMap & vertexGlobal2LocalIndexMap,
            std::map<GlobalIndex, GlobalCoordinate> & vertexIndex2CoordinateMap,
            std::vector<GmshEntityData> & internalElementVector,
            unsigned int nBoundarySegmentTotal
            )
    {
        // Write elements to factory
        for (unsigned int i = 0; i < internalElementVector.size(); i++)
        {
        	LoggingMessage::writePatience(" Inserting internal elements into factory...", i, internalElementVector.size());

            // Obtain all necessary info not to use gmshElementIndex in the following steps
            // *****************************************************
            GeometryType elemType = gmsh2dunemapper_.geometryType(internalElementVector[i].gmshType_);
            int elemOrder         = gmsh2dunemapper_.elementOrder(internalElementVector[i].gmshType_);
            int elemDofNo         = gmsh2dunemapper_.dofNumber(internalElementVector[i].gmshType_);

            int elemDim           = elemType.dim();
            int elemCornerNo      = ReferenceElements::general(elemType).size(elemDim);

            std::string log_string = "    * internal_element " + std::to_string(i) + " has dimension " + std::to_string(elemDim) + " and vertex number " + std::to_string(elemCornerNo) + " and physical entity number " + std::to_string(internalElementVector[i].physicalTag_);
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);



            LocalIndexVector localDofVector;
            GlobalCoordinateVector elementNodeVector;

            for (int iDof = 0; iDof < elemDofNo; iDof++) {
            	int vertexGlobalIndex = internalElementVector[i].vertexIndexSet_[iDof];

                localDofVector.push_back(vertexGlobal2LocalIndexMap[vertexGlobalIndex]);    // Compute local DoF vector
                elementNodeVector.push_back(vertexIndex2CoordinateMap[vertexGlobalIndex]);  // Compute vertices of this element
            }

            //std::cout << "process_" << rank_ << " element=" << i << " globalVertexIndices=(" << VectorHelper::vector2string(internalElementVector[i].vertexIndexSet_) << ")" << " localVertexIndices=(" << VectorHelper::vector2string(localDofVector) << ")" << std::endl;


            // TESTING SECTION FOR TETRAHEDRA
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (writeVtkFile_)
            {
            	addElementToVTK<dimWorld_>(elemType, elementNodeVector, elemOrder, internalElementVector[i].physicalTag_, false);

            	log_string = "    * internal_element " + std::to_string(i) + " has been added to the VTK triangles  ";
            	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
            }


            // Compute Element Global Index from gmsh entity index.
            // In GMSH, boundary segments come first as a continuous block, so they need to be subtracted
            // ******************************************************
            int globalIndex = internalElementVector[i].entityIndex_ - nBoundarySegmentTotal;


            //Insert internal element
            //****************************************************
            // Note: Global index available through internalElementVector[i].entityIndex_ is not necessary, since it is
            // shared with domain boundary triangles, and the grid requires separate global index for all entity codimensions
            factory.insertElement(elemType, localDofVector, globalIndex, elemOrder, internalElementVector[i].physicalTag_);


            if (!insertBoundarySegment_) {
                // This should not happen because CurvGridFactory demands insertion of all boundary segments
                DUNE_THROW(Dune::IOError, "You must insert boundary segments, you do not have a choice :D" );
            }

            log_string = "    * internal_element " + std::to_string(i) + " has been added to the Geometry Factory ";
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
        }
    }


    /** \brief Adds all boundary elements to factory, also writes them to .vtk file.
     *
     *  \param[in]  vtk_curv_writer              A class that writes debug output to .vtk file(s)
     *  \param[in]  vertexGlobal2LocalIndexMap   the map from vertex globalID to localID
     *  \param[in]  vertexIndex2CoordinateMap    the map from vertex globalID to vertex coordinate.
     *  \param[in]  boundaryElementVector        A vector in which the globalID's of boundary elements are stored
     *  \param[in]  linkedElementLocalIndexVector   A vector that stores a vector of localID's of all elements linked this boundary, for each boundary localID
     *
     *  [TODO] Currently factory.insertBoundarySegment() inserts index of element associated with this element, which Dune might not like:
     *  Possible solutions:
     *    * Extend Dune-interface with this function. Otherwise inserting boundary segment is pointless - having to find which element this boundary segment is associated with
     *    * Recompute in CurvGrid - not too expensive to redo, but annoying
     *    * Introduce compiler directive -DHAVE_CURVREADER_BOUNDARY_SEGMENT_ASSOCIATION
     *
     *  [TODO] The vertex and element vectors are stored twice - once inside the read procedure and once in factory
     *  as they are being added. Maybe possible to save space
     */
    void addBoundaryElements(
            Global2LocalIndexMap & vertexGlobal2LocalIndexMap,
            std::map<GlobalIndex, GlobalCoordinate> & vertexIndex2CoordinateMap,
            std::vector<GmshEntityData> & boundaryElementVector,
            std::vector<LocalIndexVector> & linkedElementLocalIndexVector,
			bool isDomainBoundary
            )
    {
        // Write elements to factory
        for (unsigned int i = 0; i < boundaryElementVector.size(); i++)
        {
        	LoggingMessage::writePatience(" Inserting boundary segments into factory...", i, boundaryElementVector.size());

            // Obtain all necessary info not to use gmshElementIndex in the following steps
            // *****************************************************
            GeometryType boundaryType = gmsh2dunemapper_.geometryType(boundaryElementVector[i].gmshType_);
            int boundaryOrder         = gmsh2dunemapper_.elementOrder(boundaryElementVector[i].gmshType_);
            int boundaryDofNo         = gmsh2dunemapper_.dofNumber(boundaryElementVector[i].gmshType_);
            int boundaryDim           = boundaryType.dim();
            int boundaryCornerNo      = SubReferenceElements::general(boundaryType).size(boundaryDim);

            std::string log_string = "    * boundary_element " + std::to_string(i) + " has dimension " + std::to_string(boundaryDim) + " and vertex number " + std::to_string(boundaryCornerNo) + " and physical entity number " + std::to_string(boundaryElementVector[i].physicalTag_);
            LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);



            LocalIndexVector localDofVector;
            GlobalCoordinateVector elementNodeVector;

            for (int iDof = 0; iDof < boundaryDofNo; iDof++) {
            	// Compute local DoF vector
                localDofVector.push_back(vertexGlobal2LocalIndexMap[boundaryElementVector[i].vertexIndexSet_[iDof]]);

                // Compute vertices of this element
                elementNodeVector.push_back(vertexIndex2CoordinateMap[boundaryElementVector[i].vertexIndexSet_[iDof]]);
            }


            // TESTING SECTION FOR TETRAHEDRA
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (writeVtkFile_)
            {
            	addElementToVTK<dimWorld_-1>(boundaryType, elementNodeVector, boundaryOrder, boundaryElementVector[i].physicalTag_, true);

            	log_string = "    * boundary_element " + std::to_string(i) + " has been added to the VTK triangles  ";
            	LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
            }
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            //Insert boundary segments and elements
            //****************************************************

            // Adding boundarySegment to factory
            // Note: Global index available through boundaryElementVector[i].entityIndex_ is not useful, since GMSH one is
            // shared between interior elements and domain boundary triangles, while the grid requires separate global index for all entity codimensions
            if (isDomainBoundary && (!insertBoundarySegment_)) {
                // This should not happen because CurvGridFactory demands insertion of all boundary segments at the moment
                DUNE_THROW(Dune::IOError, "You must insert boundary segments, you do not have a choice :D" );
            } else {
                // Note: Linked interior element is no longer necessary for curvilinear grid constructor.
            	// factory.insertBoundarySegment(boundaryType, localDofVector, boundaryOrder, linkedElementLocalIndexVector[i][0], boundaryElementVector[i].physicalTag_, isDomainBoundary);
            	factory.insertBoundarySegment(boundaryType, localDofVector, boundaryOrder, boundaryElementVector[i].physicalTag_, isDomainBoundary);

                std::string log_string = "    * boundary entity " + std::to_string(i) + " of type " +
                		((isDomainBoundary) ? "Domain Boundary" : "Interior Boundary") +
                		" has been added to the Geometry Factory ";
                LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, log_string);
            }
        }
    }


    // ***********************************************************************
    // VTK Writer Sub-procedures
    // ***********************************************************************


    /** \brief Adds an element to a vtk writer
     *
     *  \param[in]  isBoundary                Whether this is an internal or a boundary element
     *
     */
    template<int mydim>
    void addElementToVTK(const GeometryType & elemType, const GlobalCoordinateVector & elemNodeVector, const int elemOrder, const int physicalTag, const bool isBoundary)
    {
    	const unsigned int INTERIOR_TYPE = Dune::PartitionType::InteriorEntity;
    	const unsigned int BOUNDARY_TYPE = BOUNDARY_SEGMENT_PARTITION_TYPE;

    	int VTK_DISCRETIZATION_POINTS = 6;    // Sampling frequency over curved element. min=2 is linear sampling
    	bool VTK_INTERPOLATE = true;          // Whether to use lagrange interpolation or intrinsic interpolatory vertices
    	bool VTK_EXPLODE = true;              // Whether to make gaps between all elements by scaling them away from center

    	// Defines what structural purpose this element has in the grid.
    	// Different elements will have different structural tags
    	int VTK_ELEMENT_STRUCTURAL_TYPE = isBoundary ? BOUNDARY_TYPE : INTERIOR_TYPE;

    	std::vector<int> elemTags  { physicalTag, VTK_ELEMENT_STRUCTURAL_TYPE, rank_ };

    	vtkCurvWriter_.template addCurvilinearElement<mydim, mydim>(
    			elemType,
    			elemNodeVector,
    			elemTags,
    			elemOrder,
    			VTK_DISCRETIZATION_POINTS,
    			VTK_INTERPOLATE,
    			VTK_EXPLODE);
    }


    // ***********************************************************************
    // Partitioner Sub-procedures
    // ***********************************************************************

    /** \brief Partition a mesh based on corner indices.
     * Computes for each element of this process (given by globalId) to which process it should belong.
     *
     * \param[in] part            return vector. Contains process ranks indexed over elements
     * \param[in] baseElementVector    vector of elements with their data, in particular corner id's and element id.
     */
     void partitionCompute(
    		 std::vector<unsigned> & part,
    		 std::vector<GmshEntityData> & baseElementVector
     ) {

        // ****************************************************
        // Preliminaries
        // ****************************************************
#if PARMETIS_MAJOR_VERSION < 4
      typedef ::idxtype ParmetisIndexType;
      typedef float     ParmetisRealType;
#else
      typedef ::idx_t    ParmetisIndexType;
      typedef ::real_t   ParmetisRealType;
#endif

      typedef std::vector<ParmetisIndexType>  ParmetisIndexVector;
      typedef std::vector<ParmetisRealType>   ParmetisRealVector;

      GeometryType elementType = gmsh2dunemapper_.geometryType(baseElementVector[0].gmshType_);
      int elementNumber = baseElementVector.size();
      int elementDim = elementType.dim();
      int elementFaceCorners = ReferenceElements::general(elementType).size(0, 1, elementDim);

      int constrNumber;
      switch (partStrat_) {
      case LoadBalanceDefault :		constrNumber = 1;		break;
      case LoadBalanceBoundary :	constrNumber = 2;		break;
      default : DUNE_THROW(IOError, "Unexpected LoadBalance strategy");
      }

      // ****************************************************
      // Setup parameters for ParMETIS
      // ****************************************************
      ParmetisIndexType wgtflag = 2;                                  // We use different weights for each element
      ParmetisIndexType numflag = 0;                                 // we are using C-style arrays
      ParmetisIndexType ncon = constrNumber;                 // number of balance constraints
      ParmetisIndexType ncommonnodes = elementFaceCorners;            // number of nodes elements must have in common in order to be adjacent to each other
      ParmetisIndexType nparts = size_;                               // number of parts equals number of processes
      ParmetisRealVector tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
      ParmetisRealVector ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)
      ParmetisIndexType options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
      ParmetisIndexType edgecut;                                      // will store number of edges cut by partition

      // ****************************************************
      // Communicate the number of elements on each process
      // ****************************************************
      std::stringstream logstr;
      logstr << "Preparing for ParMETIS mesh partition with parameters:";
      logstr << " partitionStrat=" << partStrat_;
      logstr << " numConstraint=" << constrNumber;
      logstr << " numElement=" << baseElementVector.size();
      LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, logstr.str());

      ParmetisIndexVector elmdist;
      ParmetisIndexVector elmdist_tmp (size_, 0);

      // The index of elmdist_tmp should be the process number, the value the number of elements on each process
#if HAVE_MPI
      MPI_Comm comm = Dune::MPIHelper::getCommunicator();
      Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

      collective_comm.allgather(&elementNumber, 1, reinterpret_cast<ParmetisIndexType*>(elmdist_tmp.data()));
#endif

      // elmdist should be an incremental array whose entries are the sum of all element numbers on previous processes
      elmdist.push_back(0);
      //for (int i = 0; i < elementNumber; i++)  { elmdist.push_back(elmdist[i] + elmdist_tmp[i]); }
      for (int i = 0; i < size_; i++)  { elmdist.push_back(elmdist[i] + elmdist_tmp[i]); }

      // ****************************************************
      // Construct element weights
      // The amount of computation associated with a curvilinear element is approx. calculated:
      //  1) The number of Lagrange Polynomials interpolating the element is equal to the number of interpolation points
      //  2) The number of basis functions to interpolate the field inside should be approximately that number too (why???)
      //  3) The number of new non-zero matrix elements is approx. number of basis functions squared
      // ****************************************************
      ParmetisIndexVector elmwgt;
      for (size_t i = 0; i < elementNumber; i++) {
    	  int ord = gmsh2dunemapper_.elementOrder(baseElementVector[i].gmshType_);
    	  int complete3DBasisSize = (ord + 1) * (ord + 2) * (ord + 3) / 2;
    	  elmwgt.push_back(pow(complete3DBasisSize, 2));

    	  // In case of boundary-priority, we need additional constraint to prioritize boundary segments
    	  if (partStrat_ == LoadBalanceBoundary) {
    		  int boundaryWgt = baseElementVector[i].isOnDB_ ? 1 : 0;
    		  elmwgt.push_back(boundaryWgt);
    	  }
      }

      // ****************************************************
      // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
      // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
      // ****************************************************
      ParmetisIndexVector eptr, eind;
      int numVertices = 0;
      eptr.push_back(numVertices);

      for (size_t i = 0; i < elementNumber; i++)
      {
    	  int curNumCorners = baseElementVector[i].vertexIndexSet_.size();
    	  numVertices += curNumCorners;
    	  eptr.push_back(numVertices);

    	  for (size_t k = 0; k < curNumCorners; ++k)  { eind.push_back(baseElementVector[i].vertexIndexSet_[k] + 1); }
      }

      LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Run parmetis routine");

      std::stringstream logstr2;
      logstr2 << elmdist.size() << " " << eptr.size() << " " << eind.size() << " " << elmwgt.size();
      LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Run parmetis routine " + logstr2.str());


#if HAVE_MPI
#if PARMETIS_MAJOR_VERSION >= 4
        const int OK =
#endif
        ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(), elmwgt.data(), &wgtflag, &numflag,
                                 &ncon, &ncommonnodes, &nparts, tpwgts.data(), ubvec.data(),
                                 options, &edgecut, reinterpret_cast<ParmetisIndexType*>(part.data()), &comm);

#if PARMETIS_MAJOR_VERSION >= 4
        if (OK != METIS_OK)
          DUNE_THROW(Dune::Exception, "ParMETIS returned an error code.");
#endif
#endif

        //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Parmetis-suggested processes for elements: " + VectorHelper::vector2string(part));
    }


     /** \brief Communicate the mesh as partitioned by partitionCompute
      * First, communicates how many elements are sent to each process from this process.
      * Then, communicates all elements (only globalId's)
      *
      * \param[in] part                     Contains process ranks indexed over elements.
      * \param[in] baseElementVector             vector of elements with their data, in particular corner id's and element id.
      * \param[in] thisProcessInteriorElementIndexSet    set of globalId's of all elements present on this process
      *
      */
     void partitionCommunicate(
    		 std::vector<unsigned> & part,
    		 std::vector<GmshEntityData> & baseElementVector,
    		 GlobalIndexSet & thisProcessInteriorElementIndexSet
     )
     {
    	 typedef std::pair<unsigned, unsigned> ETP;
#if HAVE_MPI
    	 MPI_Comm comm = Dune::MPIHelper::getCommunicator();
#endif   // In sequential case this method will not be run at all, so no specialized impl


    	 // 1) Construct a vector of globalId's sorted by the corresponding process number
    	 // *****************************************************************************

    	 LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Sort golbalId array");

    	 // Combine element globalId's with processes to which these elements go
    	 std::vector<ETP> elementToProcess;
    	 for (unsigned int i = 0; i < baseElementVector.size(); i++) {
    		 elementToProcess.push_back(std::make_pair(baseElementVector[i].entityIndex_, part[i]));
    	 }
    	 // Sort according to increasing process order
    	 struct comparator {
    		 bool operator() (ETP A, ETP B) { return (A.second < B.second); }
    	 } mycomparator;
    	 std::sort(elementToProcess.begin(), elementToProcess.end(), mycomparator);


    	 // 2) Compute how many elements are send to each process. Communicate this to all processes
    	 // *****************************************************************************
    	 LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Communicate number of elements to sent do each process");

    	 std::vector<int> sendcounts (size_, 0);
    	 std::vector<int> recvcounts (size_, 0);
    	 for (unsigned int i = 0; i < part.size(); i++) {
    		 //std::cout << " step " << i << " requests process " << part[i] << " of total " << size_ << std::endl;
    		 sendcounts[part[i]] += 1;
    	 }

#if HAVE_MPI
    	 MPI_Alltoall (sendcounts.data(), 1, MPI_INT, reinterpret_cast<int*>(recvcounts.data()), 1, MPI_INT, comm);
#endif


    	 // 3) Construct send and receive displacements (sdispls)
    	 // *****************************************************************************

    	 LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Compute buffer sizes");

    	 std::vector<int> sdispls (size_, 0);
    	 std::vector<int> rdispls (size_, 0);

    	 for (int i = 1; i < size_; i++) {
    		 sdispls[i] = sdispls[i - 1] + sendcounts[i - 1];
    		 rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];
    	 }


    	 // 4) Construct sendbuffer and recvbuffer
    	 // *****************************************************************************

    	 int recvbuf_size = 0;
    	 std::vector<int> sendbuf;
    	 for (unsigned int i = 0; i < elementToProcess.size(); i++) { sendbuf.push_back(elementToProcess[i].first); }
    	 for (unsigned int i = 0; i < recvcounts.size(); i++)       { recvbuf_size += recvcounts[i]; }

    	 std::vector<int> recvbuf(recvbuf_size, 0);


    	 // 5) Communicate global indices, put place them into part vector as return value
    	 // *****************************************************************************
    	 LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Communicate globalId's");

#if HAVE_MPI
   	     MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );
#endif

   	     //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Sending buffer: " + VectorHelper::vector2string(sendbuf));
   	     //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Sending counts: " + VectorHelper::vector2string(sendcounts));
   	     //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Sending displs: " + VectorHelper::vector2string(sdispls));
   	     //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Receiving buffer: " + VectorHelper::vector2string(recvbuf));
   	     //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Receiving counts: " + VectorHelper::vector2string(recvcounts));
   	     //LoggingMessage::template write<LOG_MSG_DVERB>( __FILE__, __LINE__, " Receiving displs: " + VectorHelper::vector2string(rdispls));

   	     thisProcessInteriorElementIndexSet = std::set<int> (recvbuf.begin(), recvbuf.end());

     }








     // *******************************************************************************
    // Variables
    // *******************************************************************************

  private:
    // Parallel Implementation
    MPIHelper &mpihelper_;
    int rank_;
    int size_;

    // Grid Factory
    //GridFactory<GridType>& factory;
    //CurvilinearGridBaseFactory<GridType> & factory;
    FactoryType & factory;

    // Reading file
    std::string fileName;
    char buf_[512];

    // Constructor constants

    bool insertBoundarySegment_;   // If to insert boundary segments into factory
    bool useGmshElementIndex_;     // If to reuse the GMSH element index to construct the GMSH global index
    bool writeVtkFile_;            // If to save mesh to VTK format after reading
    bool partitionMesh_;           // If to partition mesh using parmetis
    CurvilinearGmshReaderLoadBalanceStrategy partStrat_;  // The strategy of mesh partitioning

    // Total data about the mesh
    int nVertexTotal_ = 0;
    int nEntityTotal_ = 0;
    int nInteriorElementTotal_ = 0;
    int nBoundarySegmentTotal_ = 0;
    int nDomainBoundarySegmentTotal_ = 0;
    int nInteriorBoundarySegmentTotal_ = 0;

    // Mapping from GMSH to Dune conventions
    Gmsh2DuneMapper gmsh2dunemapper_;

    // Testing capabilities for writing to VTK.
    CurvilinearVTKWriter<GridType> vtkCurvWriter_;

  };






  /**
     \ingroup Gmsh

     \brief Read Gmsh mesh file with curvilinear finite elements

     Read a .msh file generated using Gmsh and construct a grid using the grid factory interface.

     The file format used by gmsh can hold grids that are more general than the simplex grids that
     the gmsh grid generator is able to construct.  We try to read as many grids as possible, as
     long as they are valid files.  You can test this by checking whether gmsh will load the file
     and display its content.

     All grids in a gmsh file live in three-dimensional Euclidean space.  If the world dimension
     of the grid type that you are reading the file into is less than three, the remaining coordinates
     are simply ignored.
   */
  template<typename GridType>
  class CurvilinearGmshReader
  {
  public:
    typedef GridType Grid;

  public:
    // Note: This reader is intended to be used by multiple different grid managers
    // Therefore we avoid methods explicitly returning pointers to the grid


    /** \brief Reads .GMSH grid, factory provided as argument
     *  Also receives physical_tag vector for both internal and boundary elements
     * */
    template <typename FactoryType>
    static void read (FactoryType & factory,
                      const std::string& fileName,
                      MPIHelper &mpihelper,
                      bool useGmshElementIndex = true,  // If the reader will reuse the gmsh element index to create element global index. If false, the element global index will be calculated by appropriately shifting the element local index
                      bool partitionMesh = true,
					  CurvilinearGmshReaderLoadBalanceStrategy partStrat = LoadBalanceDefault,  // What weighting strategy to use for Grid
					  bool curvReaderWriteVTK = false  // The reader does not immediately output the mesh to VTK by default
    )
    {
        // [FIXME] It should not be necessary to know about boundary segments
    	// The reader should always insert them if they are available
    	//const bool DEFAULT_CURV_GMSH_READER_WRITE_VTK               = false;   // If the reader will write mesh to .vtk after reading it
    	const bool DEFAULT_CURV_GMSH_READER_INSERT_BOUNDARY_SEGMENT = true;    // If the reader will insert boundary segments

    	std::string log_string = "[[Started CurvilinearGmshReader. This rank " + std::to_string(mpihelper.rank()) + " with total processes " + std::to_string(mpihelper.size());
    	LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, log_string);

        // create parse object
        CurvilinearGmshReaderParser<Grid, FactoryType> parser(
        	factory,
        	mpihelper,
        	DEFAULT_CURV_GMSH_READER_INSERT_BOUNDARY_SEGMENT,
        	useGmshElementIndex,
			curvReaderWriteVTK,
        	partitionMesh,
			partStrat
        );

        parser.read(fileName);

        // Insert compulsory total number of vertices and elements into the curvilinear factory
        factory.insertNVertexTotal(parser.totalVertex());
        factory.insertNElementTotal(parser.totalInteriorElement());

    	LoggingMessage::template write<LOG_MSG_PRODUCTION>( __FILE__, __LINE__, "...Finished CurvilinearGmshReader]]");
    }
  };

} // namespace CurvGrid

} // namespace Dune

#endif /** DUNE_CURVILINEARGMSHREADER_HH **/
