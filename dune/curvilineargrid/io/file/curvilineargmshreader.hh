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

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


#include <dune/curvilineargrid/feedback/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargrid/curvilineargridfactory.hh>
#include <dune/curvilineargrid/io/file/curvilinearvtkwriter.hh>

#include <parmetis.h>




namespace Dune
{

  // Stores all info associated with an element, except explicit vertex coordinates
  struct GmshElementData
  {
      int elementId_;
      int gmshIndex_;
      int physicalEntityTag_;
      int elementTag_;
      std::vector<int> processTagSet_;
      std::vector<int> elementDofSet_;
  };



  //! dimension independent parts for CurvilinearGmshReaderParser
  template<typename GridType>
  class CurvilinearGmshReaderParser
  {
  protected:


    // static data
    static const int dim_ = GridType::dimension;
    static const int dimWorld_ = GridType::dimensionworld;
    static_assert( (dimWorld_ <= 3), "GmshReader requires dimWorld <= 3." );

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;

    // typedefs
    typedef Dune::FieldVector< double, dimWorld_ > GlobalVector;
    typedef Dune::ReferenceElement< double, dim_ > ReferenceElement;
    typedef Dune::ReferenceElements< double, dim_ > ReferenceElements;

    typedef Dune::ReferenceElement< double, dim_-1 > SubReferenceElement;
    typedef Dune::ReferenceElements< double, dim_-1 > SubReferenceElements;




    // Methods
    // ***********************************************************************

    // Constructs a DUNE geometry type based on GMSH element index
    GeometryType gmshGeometryType(int gmshIndex)
    {
        GeometryType rez;
        int gi = gmshIndex;

             if (gi == 15)                                                                                                     { rez.makeVertex(); }
        else if ((gi == 1) || (gi == 8) || (gi == 26) || (gi == 27) || (gi == 28))                                             { rez.makeLine(); }
        else if ((gi == 2) || (gi == 9) || (gi == 20) || (gi == 21) || (gi == 22) || (gi == 23) || (gi == 24) || (gi == 25))   { rez.makeTriangle(); }
        else if ((gi == 3) || (gi == 10) || (gi == 16))                                                                        { rez.makeQuadrilateral(); }
        else if ((gi == 4) || (gi == 11) || (gi == 29) || (gi == 30) || (gi == 31))                                            { rez.makeTetrahedron(); }
        else if ((gi == 5) || (gi == 12) || (gi == 17) || (gi == 92) || (gi == 93))                                            { rez.makeHexahedron(); }
        else if ((gi == 6) || (gi == 13) || (gi == 18))                                                                        { rez.makePrism(); }
        else if ((gi == 7) || (gi == 14) || (gi == 19))                                                                        { rez.makePyramid(); }
        else  { DUNE_THROW(Dune::IOError, "Unexpected geometry type");  }

/*             // Note that hexahedron GeometryType is missing
             switch (thisElmName)
             {
             case GMSH_EDGE           : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
             case GMSH_TRIANGLE       : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
             case GMSH_QUADRANGLE     : factory.insertElement(GeometryType(GeometryType::cube,dim_),    corners);   break;
             case GMSH_TETRAHEDRON    : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
             case GMSH_HEXAHEDRON     : factory.insertElement(GeometryType(GeometryType::simplex,dim_), corners);   break;
             case GMSH_PRISM          : factory.insertElement(GeometryType(GeometryType::prism,dim_),   corners);   break;
             case GMSH_PYRAMID        : factory.insertElement(GeometryType(GeometryType::pyramid,dim_), corners);   break;
             }*/

        return rez;
    }

    // Returns the type name of the element given its GMSH_index
    int gmshElementOrder(int gmshIndex)
    {
        // Hexahedra of high dimension have funny index
        if (gmshIndex == 92) { return 3; }
        if (gmshIndex == 93) { return 4; }

        // Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
        const int elemOrder[32]          = {-1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 3, 4, 5, 3, 4, 5};

        return elemOrder[gmshIndex];
    }

    bool gmshElementIsIncomplete(int gmshIndex)
    {
        // Both high-order hexahedrons are complete
        if ((gmshIndex == 92) || (gmshIndex == 93)) { return 0; }

        // Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
        const bool elemIncomplete[32]          = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};

        return elemIncomplete[gmshIndex];
    }

    // Returns the number of degrees of freedom of the element given its GMSH_index
    // note: This info can not simply be obtained from referenceElement, because some of the elements in GMSH have incomplete order, so less DoF than expected
    int gmshElementDofNumber(int gmshIndex)
    {
        // Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
        const int nDofs[32]              = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56};

        return nDofs[gmshIndex];
    }

    // Returns the total number of DoF associated with all subentities of a given dimension for this element, subtracting the ones that come from the corners
    int gmshElementSubentityExtraDim(int gmshIndex, int dim)
    {
        const int nDofsExtraEdge[32] = {-1, 0, 0, 0, 0, 0, 0, 0, 1, 3, 4, 6, 12, 9, 8, 0, 4, 12, 9, 8, 6, 6, 9, 9, 12, 12, 2, 3, 4, 12, 18, 24};
        const int nDofsExtraFace[32] = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 6, 3, 1, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 6, 0, 0, 0, 4, 12, 24};
        const int nDofsExtraElem[32] = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4};

        switch(dim)
        {
        case 1: return nDofsExtraEdge[gmshIndex];  break;
        case 2: return nDofsExtraFace[gmshIndex];  break;
        case 3: return nDofsExtraElem[gmshIndex];  break;
        }

        return -1;
    }

    // correct differences between gmsh and Dune in the local vertex numbering
    // NOTE: THIS METHOD DOES NOT WORK WITH INCOMPLETE ORDER GMSH ELEMENTS AT THE MOMENT
    void gmsh2DuneElementDofNumbering(GeometryType gt, int thisElmOrder, std::vector<int> &elementDofSet) {
        int thisElmDofNo = elementDofSet.size();
        std::vector<int> tmp;

        if (!gt.isSimplex() || gt.dim() < 2 || gt.dim() > 3)  {
        	DUNE_THROW(Dune::IOError, "CURVILINEAR_GMSH_READER: gmsh2duneRenumbering() only implemented for Simplex 2D and 3D geometries at the moment");
        }

        if (gt.isTriangle())
        {
            for (int i = 0; i < thisElmDofNo; i++) { tmp.push_back(elementDofSet[triangularInterpolatoryVertexGmsh2DuneMap[thisElmOrder - 1][i]]); }
            for (int i = 0; i < thisElmDofNo; i++) { elementDofSet[i] = tmp[i]; }
        }
        else if (gt.isTetrahedron())
        {
            for (int i = 0; i < thisElmDofNo; i++) { tmp.push_back(elementDofSet[tetrahedralInterpolatoryVertexGmsh2DuneMap[thisElmOrder - 1][i]]); }
            for (int i = 0; i < thisElmDofNo; i++) { elementDofSet[i] = tmp[i]; }
        }
    }



    // Testing if the current element type can be handled by DUNE
    bool checkElementAllowed(int gmshIndex)
    {
    	GeometryType gt = gmshGeometryType(gmshIndex);
        bool isAllowedElement = true;

        // Only allow simplex geometries for now
        isAllowedElement &= gt.isSimplex();

        // Check if element is polynomial-complete (ask GMSH what that even means I dont know :) )
        isAllowedElement &= !gmshElementIsIncomplete(gmshIndex);

        // test whether we support the element type at the moment
        if (!isAllowedElement) { DUNE_THROW(Dune::IOError, "GMSH Reader: Have read an element of unexpected type "); }

        return isAllowedElement;
    }

    // Checks whether this element (or boundary element) belongs on this parallel process
    bool elementOnProcess(int eIndex, int eTotal) {
        int eFirst = (eTotal * rank_) / size_;
        int eLast = (eTotal * (rank_+1)) / size_;

        std::string log_string = " == checkprocess if " + std::to_string(eIndex) + " in [" + std::to_string(eFirst) + "," + std::to_string(eLast) + "]";
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

        return ((eIndex >= eFirst)&&(eIndex < eLast));
    }


    /** \brief Partition a mesh based on corner indices.
     * Computes for each element of this process (given by globalId) to which process it should belong.
     *
     * \param[in] part            return vector. Contains process ranks indexed over elements
     * \param[in] baseElementVector    vector of elements with their data, in particular corner id's and element id.
     */
     void partitionCompute(
    		 std::vector<unsigned> & part,
    		 std::vector<GmshElementData> & baseElementVector
     ) {

        // ****************************************************
        // Preliminaries
        // ****************************************************
#if PARMETIS_MAJOR_VERSION < 4
      typedef idxtype idx_t;
      typedef float real_t;
#endif

      GeometryType elementType = gmshGeometryType(baseElementVector[0].gmshIndex_);
      int elementNumber = baseElementVector.size();
      int elementDim = elementType.dim();
      int elementFaceCorners = ReferenceElements::general(elementType).size(0, 1, elementDim);

      // ****************************************************
      // Setup parameters for ParMETIS
      // ****************************************************
      idx_t wgtflag = 2;                                  // We use different weights for each element
      idx_t numflag = 0;                                  // we are using C-style arrays
      idx_t ncon = 1;                                     // number of balance constraints
      idx_t ncommonnodes = elementFaceCorners;            // number of nodes elements must have in common in order to be adjacent to each other
      idx_t nparts = size_;                               // number of parts equals number of processes
      std::vector<real_t> tpwgts(ncon*nparts, 1./nparts); // load per subdomain and weight (same load on every process)
      std::vector<real_t> ubvec(ncon, 1.05);              // weight tolerance (same weight tolerance for every weight there is)
      idx_t options[4] = {0, 0, 0, 0};                    // use default values for random seed, output and coupling
      idx_t edgecut;                                      // will store number of edges cut by partition

      // ****************************************************
      // Communicate the number of elements on each process
      // ****************************************************
      Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Communicating element numbers on all processes to each process");

      std::vector<idx_t> elmdist;
      std::vector<idx_t> elmdist_tmp (size_, 0);

      // The index of elmdist_tmp should be the process number, the value the number of elements on each process
#if HAVE_MPI
      MPI_Comm comm = Dune::MPIHelper::getCommunicator();
      Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

      collective_comm.allgather(&elementNumber, 1, reinterpret_cast<idx_t*>(elmdist_tmp.data()));
#endif

      // elmdist should be an incremental array whose entries are the sum of all element numbers on previous processes
      elmdist.push_back(0);
      for (int i = 0; i < elementNumber; i++)  { elmdist.push_back(elmdist[i] + elmdist_tmp[i]); }

      // ****************************************************
      // Construct element weights
      // The amount of computation associated with a curvilinear element is approx. calculated:
      //  1) The number of Lagrange Polynomials interpolating the element is equal to the number of interpolation points
      //  2) The number of basis functions to interpolate the field inside should be approximately that number too (why???)
      //  3) The number of new non-zero matrix elements is approx. number of basis functions squared
      // ****************************************************
      std::vector<idx_t> elmwgt;
      for (size_t i = 0; i < elementNumber; i++) {
    	  int elementOrder = gmshElementOrder(baseElementVector[i].gmshIndex_);
    	  elmwgt.push_back(pow(elementOrder, 2));
      }

      // ****************************************************
      // Create and fill arrays "eptr", where eptr[i] is the number of vertices that belong to the i-th element, and
      // "eind" contains the vertex-numbers of the i-the element in eind[eptr[i]] to eind[eptr[i+1]-1]
      // ****************************************************
      std::vector<idx_t> eptr, eind;
      int numVertices = 0;
      eptr.push_back(numVertices);

      for (size_t i = 0; i < elementNumber; i++)
      {
    	  int curNumCorners = baseElementVector[i].elementDofSet_.size();
    	  numVertices += curNumCorners;
    	  eptr.push_back(numVertices);

    	  for (size_t k = 0; k < curNumCorners; ++k)  { eind.push_back(baseElementVector[i].elementDofSet_[k]); }
      }

      Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Run parmetis routine");


#if HAVE_MPI
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
#endif

        //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Parmetis-suggested processes for elements: " + vector2string(part));
    }


     /** \brief Communicate the mesh as partitioned by partitionCompute
      * First, communicates how many elements are sent to each process from this process.
      * Then, communicates all elements (only globalId's)
      *
      * \param[in] part                     Contains process ranks indexed over elements.
      * \param[in] baseElementVector             vector of elements with their data, in particular corner id's and element id.
      * \param[in] thisProcessElementIndexSet    set of globalId's of all elements present on this process
      *
      */
     void partitionCommunicate(
    		 std::vector<unsigned> & part,
    		 std::vector<GmshElementData> & baseElementVector,
    		 std::set<int> & thisProcessElementIndexSet
     )
     {
    	 typedef std::pair<unsigned, unsigned> ETP;
#if HAVE_MPI
    	 MPI_Comm comm = Dune::MPIHelper::getCommunicator();
#endif   // In sequential case this method will not be run at all, so no specialized impl


    	 // 1) Construct a vector of globalId's sorted by the corresponding process number
    	 // *****************************************************************************

    	 Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Sort golbalId array");

    	 // Combine element globalId's with processes to which these elements go
    	 std::vector<ETP> elementToProcess;
    	 for (int i = 0; i < baseElementVector.size(); i++) {
    		 elementToProcess.push_back(std::make_pair(baseElementVector[i].elementId_, part[i]));
    	 }
    	 // Sort according to increasing process order
    	 struct comparator {
    		 bool operator() (ETP A, ETP B) { return (A.second < B.second); }
    	 } mycomparator;
    	 std::sort(elementToProcess.begin(), elementToProcess.end(), mycomparator);


    	 // 2) Compute how many elements are send to each process. Communicate this to all processes
    	 // *****************************************************************************
    	 Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Communicate number of elements to sent do each process");

    	 std::vector<int> sendcounts (size_, 0);
    	 std::vector<int> recvcounts (size_, 0);
    	 for (int i = 0; i < part.size(); i++) {
    		 //std::cout << " step " << i << " requests process " << part[i] << " of total " << size_ << std::endl;
    		 sendcounts[part[i]] += 1;
    	 }

#if HAVE_MPI
    	 MPI_Alltoall (sendcounts.data(), 1, MPI_INT, reinterpret_cast<int*>(recvcounts.data()), 1, MPI_INT, comm);
#endif


    	 // 3) Construct send and receive displacements (sdispls)
    	 // *****************************************************************************

    	 Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Compute buffer sizes");

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
    	 for (int i = 0; i < elementToProcess.size(); i++) { sendbuf.push_back(elementToProcess[i].first); }
    	 for (int i = 0; i < recvcounts.size(); i++)       { recvbuf_size += recvcounts[i]; }

    	 std::vector<int> recvbuf(recvbuf_size, 0);


    	 // 5) Communicate global indices, put place them into part vector as return value
    	 // *****************************************************************************
    	 Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Communicate globalId's");

#if HAVE_MPI
   	     MPI_Alltoallv (sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvbuf.data()), recvcounts.data(), rdispls.data(), MPI_INT, comm );
#endif

   	     //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Sending buffer: " + vector2string(sendbuf));
   	     //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Sending counts: " + vector2string(sendcounts));
   	     //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Sending displs: " + vector2string(sdispls));
   	     //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Receiving buffer: " + vector2string(recvbuf));
   	     //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Receiving counts: " + vector2string(recvcounts));
   	     //Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Receiving displs: " + vector2string(rdispls));

   	     thisProcessElementIndexSet = std::set<int> (recvbuf.begin(), recvbuf.end());

     }


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
     *  TODO: Make factory.insertVertex() work with globalId
     *
     */
    int readVertices(
            FILE* file,
            int nVertexTotal_,
            std::map<int, GlobalVector> & vertexIndex2CoordinateMap,
            std::set<int> & vertexIndexSet,
            std::map<int, int> & vertexGlobal2LocalIndexMap
            )
    {
        int id;
        GlobalVector x;

        // Iterator starts from 1 because GMSH numbers vertices [1,n]
        for( int i = 1; i <= nVertexTotal_; ++i )
        {
            // If this vertex does not belong to this process, just skip it
            if (vertexIndexSet.count(i) == 0)  { fgets(buf_, 512, file ); }
            else
            {
                fscanf(file, "%d ", &id);
                std::string tmp_out;

                if( id != i )  { DUNE_THROW( Dune::IOError, "Expected id " << i << "(got id " << id << "." ); }
                for (int d = 0; d < dimWorld_; d++) {
                    double tmp_coord;
                    fscanf(file, "%lg", &tmp_coord);
                    x[d] = tmp_coord;

                    tmp_out += std::to_string(tmp_coord) + " ";
                }

                // The local id of this vertex is equalt to the number of vertices added so far
                vertexGlobal2LocalIndexMap[i] = vertexIndex2CoordinateMap.size();

                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "  * Have read vertex " + tmp_out);

                // Maps global id to global coordinate
                vertexIndex2CoordinateMap[i] = x;

                // Insert vertex into a factory, noting its globalId.
                // Its localId in the factory should be given by the order the vertices are added to the factory
                factory.insertVertex(x, i);

                fscanf(file, "\n");
            }
        }

        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndNodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndNodes"); }

        return vertexIndex2CoordinateMap.size();
    }


    // Reads the data about this element, everything except the interpolation vertex id's
    // Assumes it is in the correct position in the file
    GmshElementData readElementSpec(FILE* file)
    {
        int nTag;
        GmshElementData thisElement;

        fscanf(file, "%d %d %d ", &thisElement.elementId_, &thisElement.gmshIndex_, &nTag);

        std::string log_string =  "    * element " + std::to_string(thisElement.elementId_) + " has " + std::to_string(nTag) + " tags";
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

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
        fscanf(file, "%d ", &thisElement.physicalEntityTag_);
        fscanf(file, "%d ", &thisElement.elementTag_);

        // Possible functionality for mesh partitioning using GMSH
        // TODO: Currently this functionality not available
        for (int k = 2; k < nTag; k++)
        {
            int tmp_tag;
            fscanf(file, "%d ", &tmp_tag);
            thisElement.processTagSet_.push_back(tmp_tag);
        }

        return thisElement;

    }


    /** \brief Communicate the mesh as partitioned by partitionCompute
     * First, communicates how many elements are sent to each process from this process.
     * Then, communicates all elements (only globalId's)
     *
     * \param[in]  file                           file pointer to read from
     * \param[in]  nElementTotal_                 the total number of elements specified in the file
     * \param[in]  nInternalElementTotal_        the number of internal elements specified in the file
     * \param[in]  thisProcessElementIndexSet           set of globalId's of all elements present on this process
     *
     */
    void readAndPartitionBaseElements(
            FILE* file,
            int nElementTotal_,
            int nInternalElementTotal_,
    		std::set<int> & thisProcessElementIndexSet)
    {
        int iSelectElem = 0;
        std::vector<GmshElementData> baseElementVector;

        // Reading element info - tag information and vertex global indices
        // *************************************************************
        for (int i = 1; i <= nElementTotal_; i++)
        {
            // Read the first part of the element info
            GmshElementData thisElement = readElementSpec(file);

            // Find if this is a boundary element
            GeometryType elemType          = gmshGeometryType(thisElement.gmshIndex_);
            int elemDim = elemType.dim();
            bool onBoundary = (elemDim < dimWorld_);

            // Check if we want to read this element at this stage
            if (onBoundary || !elementOnProcess(iSelectElem++, nInternalElementTotal_)) { fgets(buf_, 512, file ); }
            else
            {
                // Testing if the current element type can be handled by DUNE
                // *****************************************************************
                checkElementAllowed(thisElement.gmshIndex_);
                std::string log_string = "    * element " + std::to_string(thisElement.elementId_) + " can be treated by Dune grid ";
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                // Obtain all necessary info not to use gmshElementIndex in the following steps
                // *****************************************************
                int thisElmOrder               = gmshElementOrder(thisElement.gmshIndex_);
                int thisElmDofNo               = gmshElementDofNumber(thisElement.gmshIndex_);
                int thisElmCorners             = ReferenceElements::general(elemType).size(elemType.dim());
                int thisElmSubentities         = ReferenceElements::general(elemType).size(1);

                // Reading Corners of the Element
                // Note: In GMSH notation corners go first
                // *************************************************
                for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) {
                    int tmpVertexGlobalId;
                    fscanf(file, "%d", &tmpVertexGlobalId);

                    log_string = "  --- have read corner " + std::to_string(tmpVertexGlobalId);
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                    thisElement.elementDofSet_.push_back(tmpVertexGlobalId);
                }
                fgets(buf_, 512, file );

                // Store elements just read
                baseElementVector.push_back(thisElement);
            }
        }

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Finished reading base elements");


        //std::cout << "Process " << rank_ << ": elements before partition: ";
        //for (int i = 0; i < baseElementVector.size(); i++) { std::cout << baseElementVector[i].elementId_ << " "; }
        //std::cout << std::endl;

#if HAVE_MPI
        // Partition the elements
        // Communicate the partitioning
        // Fill the element set
        // *************************************************************
    	std::vector<unsigned> part(baseElementVector.size());
    	partitionCompute(part, baseElementVector);
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Finished computing partition");

    	partitionCommunicate(part, baseElementVector, thisProcessElementIndexSet);
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " Finished communicating partition");

    	std::vector<int> test2 (thisProcessElementIndexSet.begin(), thisProcessElementIndexSet.end());
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " elements after partition: " + vector2string(test2));
#else
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "No MPI found! Running sequential case without partitioning");
    	for (int i = 0; i < baseElementVector.size(); i++) { thisProcessElementIndexSet.insert(baseElementVector[i].elementId_); }
#endif
        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    /** \brief Reads all internal element data into a vector. Only reads elements which should be on this process
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  nElementTotal_                 the total number of elements specified in the file
     *  \param[in]  nInternalElementTotal_        the number of internal elements specified in the file
     *  \param[in]  thisProcessElementIndexSet          set of globalId's of all elements present on this process
     *  \param[in]  internalElementVector              A vector in which the globalID's of internal elements will be stored
     *  \param[in]  vertexIndexSet               the set of globalID's of all vertices that belong to this process
     *  \param[in]  boundaryKey2LinkedElementSet        A map from a sorted array of globalID's of vertices that make up a boundary to an array of localID's of internal elements to whom this boundary belongs.
     *
     *  \note Assumes it is in the correct position in the file
     *
     *  TODO: Make factory.insertVertex() work with globalId
     */
    void readInternalElements(
            FILE* file,
            int nElementTotal_,
            int nInternalElementTotal_,
            std::set< int > & thisProcessElementIndexSet,
            std::vector<GmshElementData> & internalElementVector,
            std::set<int> & vertexIndexSet,
            std::map< std::vector<int>, std::vector<int> > & boundaryKey2LinkedElementSet
            )
    {
    	std::string log_string;

        // Reading element info - tag information and vertex global indices
        for (int i = 1; i <= nElementTotal_; i++)
        {
            // Read the first part of the element info
            GmshElementData thisElement = readElementSpec(file);

            // Find if this is a boundary element
            GeometryType elemType          = gmshGeometryType(thisElement.gmshIndex_);
            int elemDim = elemType.dim();
            bool onBoundary = (elemDim < dimWorld_);

            // Check if this element belongs on this process
            // Note: There is no need to check if the element is a boundary segment, because
            if (thisProcessElementIndexSet.count(thisElement.elementId_) == 0) { fgets(buf_, 512, file ); }
            else
            {
                // Testing if the current element type can be handled by DUNE
                // *****************************************************************
                checkElementAllowed(thisElement.gmshIndex_);
                log_string = "    * element " + std::to_string(thisElement.elementId_) + " can be treated by Dune grid ";
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                // Obtain all necessary info not to use gmshElementIndex in the following steps
                // *****************************************************
                int thisElmOrder               = gmshElementOrder(thisElement.gmshIndex_);
                int thisElmDofNo               = gmshElementDofNumber(thisElement.gmshIndex_);
                int thisElmCorners             = ReferenceElements::general(elemType).size(elemType.dim());
                int thisElmSubentities         = ReferenceElements::general(elemType).size(1);

                // Reading DoF's
                // *************************************************
                for (int iDof = 0; iDof < thisElmDofNo; iDof++) {
                    int tmpVertexGlobalId;
                    fscanf(file, "%d", &tmpVertexGlobalId);

                    log_string = "  --- have read DoF " + std::to_string(tmpVertexGlobalId);
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                    thisElement.elementDofSet_.push_back(tmpVertexGlobalId);

                    // Insert all used global vertex indexes into set
                    // Note: set ignores request to add an already existing element
                    vertexIndexSet.insert(tmpVertexGlobalId);
                }
                fscanf(file, "\n");


                // Add all subentities of this element to a map, such that it is easy afterwards to find
                // boundary elements associated with this process
                // *************************************************************************************

                // 1) get all corners of this process
                // Note that in GMSH convention the corners come first in the Dof array
                std::vector<int> cornerVector;
                for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) { cornerVector.push_back(thisElement.elementDofSet_[iCorner]); }

                // 2) sort by increasing globalID
                std::sort(cornerVector.begin(), cornerVector.end());

                // 3) Calculate the localID of this element
                int localID = internalElementVector.size();

                for (int iSub = 0; iSub < thisElmSubentities; iSub++)
                {
                    // 4) get all subsets associated with this element type
                    std::vector<int> thisSubentityCornerIndexSet = Dune::CurvilinearGeometryHelper::linearElementSubentityCornerInternalIndexSet(elemType, 1, iSub);
                    std::vector<int> key;

                    log_string = " for iSub = " + std::to_string(iSub) + " out of " + std::to_string(thisElmSubentities) + " have sub_size = " + std::to_string(thisSubentityCornerIndexSet.size());
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                    for (int iCoord = 0; iCoord < thisSubentityCornerIndexSet.size(); iCoord++) { key.push_back(cornerVector[thisSubentityCornerIndexSet[iCoord]]); }

                    // 5) Check if map empty for this entry, then add to the map
                    if (boundaryKey2LinkedElementSet.find(key) == boundaryKey2LinkedElementSet.end()) {
                    	log_string = " -- element " + std::to_string(localID) + " has added boundary " + vector2string(key);
                        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                        boundaryKey2LinkedElementSet[key] = std::vector<int> (1, localID);
                    } else
                    {
                        // 6) In this case this is the 2nd element sharing this boundary
                        std::vector<int>  tmp = boundaryKey2LinkedElementSet[key];
                        tmp.push_back(localID);
                        boundaryKey2LinkedElementSet[key] = tmp;   // Should overwrite previous value

                        log_string = " -- element " + std::to_string(localID) + " shares boundary " + vector2string(key) + " with element " + std::to_string(tmp[0]);
                        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
                    }
                }

                // correct differences between gmsh and Dune in the local vertex numbering
                // *************************************************
                gmsh2DuneElementDofNumbering(elemType, thisElmOrder, thisElement.elementDofSet_);

                internalElementVector.push_back(thisElement);
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
     *  \param[in]  nElementTotal_                 the total number of elements specified in the file
     *  \param[in]  boundaryElementVector              A vector in which the globalID's of boundary elements will be stored
     *  \param[in]  boundaryKey2LinkedElementSet        A map from a sorted array of globalID's of vertices that make up a boundary to an array of localID's of internal elements to whom this boundary belongs.
     *  \param[in]  linkedElementLocalIndexSet   A vector that stores a vector of localID's of all elements linked this boundary, for each boundary localID
     *
     *  \note Assumes it is in the correct position in the file
     */
    void readBoundaryElements(
            FILE* file,
            int nElementTotal_,
            std::vector<GmshElementData> & boundaryElementVector,
            std::map< std::vector<int>, std::vector<int> > & boundaryKey2LinkedElementSet,
            std::vector< std::vector<int> > & linkedElementLocalIndexSet
            )
    {
    	std::string log_string;

        int iSelectElem = 0;

        // Reading element info - tag information and vertex global indices
        for (int i = 1; i <= nElementTotal_; i++)
        {
            // Read the first part of the element info
            GmshElementData thisElement = readElementSpec(file);

            // Find if this is a boundary element
            GeometryType elemType          = gmshGeometryType(thisElement.gmshIndex_);
            int elemDim = elemType.dim();
            bool onBoundary = (elemDim < dimWorld_);

            // If this element is not on the boundary just skip the rest of info on this line
            if (!onBoundary) { fgets(buf_, 512, file ); }
            else
            {
                // Testing if the current element type can be handled by DUNE
                // *****************************************************************
                checkElementAllowed(thisElement.gmshIndex_);
                log_string = "    * element " + std::to_string(i) + " can be treated by Dune grid ";
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                // Obtain all necessary info not to use gmshElementIndex in the following steps
                // *****************************************************
                int thisElmOrder               = gmshElementOrder(thisElement.gmshIndex_);
                int thisElmDofNo               = gmshElementDofNumber(thisElement.gmshIndex_);
                int thisElmCorners             = SubReferenceElements::general(elemType).size(elemType.dim());

                // Reading DoF's
                // *************************************************
                for (int iDof = 0; iDof < thisElmDofNo; iDof++) {
                    int tmpVertexGlobalId;
                    fscanf(file, "%d", &tmpVertexGlobalId);

                    log_string = "  --- have read DoF " + std::to_string(tmpVertexGlobalId);
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                    thisElement.elementDofSet_.push_back(tmpVertexGlobalId);
                }
                fscanf(file, "\n");



                // Check if associated with any internal element
                // ****************************************************************

                // 1) Obtain corners
                std::vector<int> cornerVector;
                for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) { cornerVector.push_back(thisElement.elementDofSet_[iCorner]); }

                // 2) sort by increasing globalID
                std::sort(cornerVector.begin(), cornerVector.end());

                // 3) Obtain this boundary element's localID
                int localID = boundaryElementVector.size();

                log_string = " -- boundary " + std::to_string(localID) + " checking key " + vector2string(cornerVector);
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);



                // If this boundary already exists in the map, then add this element
                // If not, it is likely on another process
                if (boundaryKey2LinkedElementSet.find(cornerVector) != boundaryKey2LinkedElementSet.end())
                {
                	log_string = " -found b.e localID = " + vector2string(boundaryKey2LinkedElementSet[cornerVector]);
                    Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

                    // correct differences between gmsh and Dune in the local vertex numbering
                    // *************************************************
                    gmsh2DuneElementDofNumbering(elemType, thisElmOrder, thisElement.elementDofSet_);

                    // Add boundary element
                    boundaryElementVector.push_back(thisElement);

                    // If this boundary is linked to an element, write its localID to the internal_element_boundaries for corresponding element
                    std::vector<int> thisBoundaryLinkedElements = boundaryKey2LinkedElementSet[cornerVector];
                    if (thisBoundaryLinkedElements.size() > 1) { DUNE_THROW(Dune::IOError, "Have 2 boundaries associated with 1 element. Not implemented yet"); }

                    linkedElementLocalIndexSet.push_back(thisBoundaryLinkedElements);
                }
            }
        }

        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    /** \brief Adds all internal elements to factory, also writes them to .vtk file.
     *
     *  \param[in]  vtk_curv_writer                A class that writes debug output to .vtk file(s)
     *  \param[in]  vertexGlobal2LocalIndexMap  the map from vertex globalID to localID
     *  \param[in]  vertexIndex2CoordinateMap                       the map from vertex globalID to vertex coordinate.
     *  \param[in]  internalElementVector              A vector in which the globalID's of internal elements are stored
     *
     *  TODO: The vertex and element vectors are stored twice - once inside the read procedure and once in factory
     *  as they are being added. Maybe possible to save space
     */
    void addInternalElements(
            std::map<int, int> & vertexGlobal2LocalIndexMap,
            std::map<int, GlobalVector> & vertexIndex2CoordinateMap,
            std::vector< GmshElementData > & internalElementVector
            )
    {
        // Write elements to factory
        for (int i = 0; i < internalElementVector.size(); i++)
        {
              // Obtain all necessary info not to use gmshElementIndex in the following steps
            // *****************************************************
            GeometryType elemType = gmshGeometryType(internalElementVector[i].gmshIndex_);
            int elemOrder               = gmshElementOrder(internalElementVector[i].gmshIndex_);
            int elemDofNo               = gmshElementDofNumber(internalElementVector[i].gmshIndex_);

            int elemDim                 = elemType.dim();
            int elemCornerNo            = ReferenceElements::general(elemType).size(elemDim);

            std::string log_string = "    * internal_element " + std::to_string(i) + " has dimension " + std::to_string(elemDim) + " and vertex number " + std::to_string(elemCornerNo) + " and physical entity number " + std::to_string(internalElementVector[i].physicalEntityTag_);
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);



            std::vector<int> localDofVector;
            std::vector<GlobalVector> elementNodeVector;

            for (int iDof = 0; iDof < elemDofNo; iDof++) {
            	// Compute local DoF vector
                localDofVector.push_back(vertexGlobal2LocalIndexMap[internalElementVector[i].elementDofSet_[iDof]]);

                // Compute vertices of this element
                elementNodeVector.push_back(vertexIndex2CoordinateMap[internalElementVector[i].elementDofSet_[iDof]]);
            }

            //std::cout << "process_" << rank_ << " element=" << i << " globalVertexIndices=(" << vector2string(internalElementVector[i].elementDofSet_) << ")" << " localVertexIndices=(" << vector2string(localDofVector) << ")" << std::endl;


            // TESTING SECTION FOR TETRAHEDRA
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (writeVtkFile_)
            {
            	addElementToVTK(elemType, elementNodeVector, elemOrder, internalElementVector[i].physicalEntityTag_, false);

            	log_string = "    * internal_element " + std::to_string(i) + " has been added to the VTK triangles  ";
            	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
            }



            //Insert internal element
            //****************************************************

            internalElement2PhysicalEntityIndex.push_back(internalElementVector[i].physicalEntityTag_);
            // Insert element to factory
            factory.insertElement(elemType, internalElementVector[i].elementId_, localDofVector, elemOrder, internalElementVector[i].physicalEntityTag_);


            if (!insertBoundarySegment) {
                // This should not happen because CurvGridFactory demands insertion of all boundary segments
                DUNE_THROW(Dune::IOError, "You must insert boundary segments, you do not have a choice :D" );
            }

            log_string = "    * internal_element " + std::to_string(i) + " has been added to the Geometry Factory ";
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
        }
    }


    /** \brief Adds all boundary elements to factory, also writes them to .vtk file.
     *
     *  \param[in]  vtk_curv_writer                A class that writes debug output to .vtk file(s)
     *  \param[in]  vertexGlobal2LocalIndexMap  the map from vertex globalID to localID
     *  \param[in]  vertexIndex2CoordinateMap                       the map from vertex globalID to vertex coordinate.
     *  \param[in]  boundaryElementVector              A vector in which the globalID's of boundary elements are stored
     *  \param[in]  linkedElementLocalIndexSet   A vector that stores a vector of localID's of all elements linked this boundary, for each boundary localID
     *
     *  TODO: The vertex and element vectors are stored twice - once inside the read procedure and once in factory
     *  as they are being added. Maybe possible to save space
     */
    void addBoundaryElements(
            std::map<int, int> & vertexGlobal2LocalIndexMap,
            std::map<int, GlobalVector> & vertexIndex2CoordinateMap,
            std::vector< GmshElementData > & boundaryElementVector,
            std::vector< std::vector<int> > & linkedElementLocalIndexSet
            )
    {
        // Write elements to factory
        for (int i = 0; i < boundaryElementVector.size(); i++)
        {
            // Obtain all necessary info not to use gmshElementIndex in the following steps
            // *****************************************************
            GeometryType boundaryType = gmshGeometryType(boundaryElementVector[i].gmshIndex_);
            int boundaryOrder               = gmshElementOrder(boundaryElementVector[i].gmshIndex_);
            int boundaryDofNo               = gmshElementDofNumber(boundaryElementVector[i].gmshIndex_);
            int boundaryDim                 = boundaryType.dim();
            int boundaryCornerNo            = SubReferenceElements::general(boundaryType).size(boundaryDim);

            std::string log_string = "    * boundary_element " + std::to_string(i) + " has dimension " + std::to_string(boundaryDim) + " and vertex number " + std::to_string(boundaryCornerNo) + " and physical entity number " + std::to_string(boundaryElementVector[i].physicalEntityTag_);
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);



            std::vector<int> localDofVector;
            std::vector<GlobalVector> elementNodeVector;

            for (int iDof = 0; iDof < boundaryDofNo; iDof++) {
            	// Compute local DoF vector
                localDofVector.push_back(vertexGlobal2LocalIndexMap[boundaryElementVector[i].elementDofSet_[iDof]]);

                // Compute vertices of this element
                elementNodeVector.push_back(vertexIndex2CoordinateMap[boundaryElementVector[i].elementDofSet_[iDof]]);
            }


            // TESTING SECTION FOR TETRAHEDRA
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (writeVtkFile_)
            {
            	addElementToVTK(boundaryType, elementNodeVector, boundaryOrder, boundaryElementVector[i].physicalEntityTag_, true);

            	log_string = "    * boundary_element " + std::to_string(i) + " has been added to the VTK triangles  ";
            	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
            }
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            //Insert boundary segments and elements
            //****************************************************

            // Adding boundarySegment to factory
            if (insertBoundarySegment)
            {
                factory.insertBoundarySegment(boundaryType, boundaryElementVector[i].elementId_, localDofVector, boundaryOrder, linkedElementLocalIndexSet[i][0], boundaryElementVector[i].physicalEntityTag_);

                // Adding physical tag
                boundaryElement2PhysicalEntityIndex.push_back(boundaryElementVector[i].physicalEntityTag_);

                log_string = "    * boundary_element " + std::to_string(i) + " has been added to the Geometry Factory ";
                Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);
            }
            else
            {
                // This should not happen because CurvGridFactory demands insertion of all boundary segments
                DUNE_THROW(Dune::IOError, "You must insert boundary segments, you do not have a choice :D" );
            }
        }
    }

    /** \brief Adds an element to a vtk writer
     *
     *  \param[in]  isBoundary                Whether this is an internal or a boundary element
     *
     */
    void addElementToVTK(const GeometryType & elemType, const std::vector<GlobalVector> & elementNodeVector, const int elemOrder, const int physicalTag, const bool isBoundary)
    {
    	const int VTK_INTERNAL = Dune::VtkEntityStructuralType::Internal;
    	const int VTK_BOUNDARY = Dune::VtkEntityStructuralType::DomainBoundary;

    	int VTK_DISCRETIZATION_POINTS = 2;    // Sampling frequency over curved element. min=2 is linear sampling
    	bool VTK_INTERPOLATE = true;          // Whether to use lagrange interpolation or intrinsic interpolatory vertices
    	bool VTK_EXPLODE = true;              // Whether to make gaps between all elements by scaling them away from center
    	bool VTK_WRITE_EDGES = false;         // Whether to write volume discretization edges to file
    	bool VTK_WRITE_TRIANGLES = true;      // Whether to write surface discretization triangles to file

    	// Defines what structural purpose this element has in the grid.
    	// Different elements will have different structural tags
    	int VTK_ELEMENT_STRUCTURAL_TYPE = isBoundary ? VTK_BOUNDARY : VTK_INTERNAL;

    	vtkCurvWriter_.addCurvilinearElement(
    			elemType,
    			elementNodeVector,
    			elemOrder,
    			physicalTag,
    			VTK_DISCRETIZATION_POINTS,
    			VTK_INTERPOLATE,
    			VTK_EXPLODE,
    			VTK_WRITE_EDGES,
    			VTK_WRITE_TRIANGLES,
    			VTK_ELEMENT_STRUCTURAL_TYPE);
    }

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::string tmp_str;
        for (int i = 0; i < V.size(); i++) { tmp_str += std::to_string(V[i]) + " "; }
        return tmp_str;
    }

  public:

    CurvilinearGmshReaderParser(
    	Dune::CurvilinearGridFactory<GridType> & _factory,
    	bool verbose,
    	bool processVerbose,
    	bool insertBoundarySegment,
    	bool writeVtkFile,
    	MPIHelper &mpihelper) :
    		factory(_factory),
    		verbose_(verbose),
    		processVerbose_(processVerbose),
    		writeVtkFile_(writeVtkFile),
    		insertBoundarySegment(insertBoundarySegment),
    		mpihelper_ (mpihelper),
    		vtkCurvWriter_(verbose_, processVerbose_, mpihelper_)
    {
        // Initialize process parameters
        rank_=mpihelper.rank();
        size_=mpihelper.size();

        std::string log_string = "I am process " + std::to_string(rank_) + " with total processes " + std::to_string(size_);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, log_string);

        // Initialize triangular point renumberings for GMSH->DUNE convention
        triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 1, 2} );
        triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 1, 5, 4, 2} );
        triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 4, 1, 8, 9, 5, 7, 6, 2} );
        triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 4, 5, 1, 11, 12, 13, 6, 10, 14, 7, 9, 8, 2} );
        triangularInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 4, 5, 6, 1, 14, 15, 18, 16, 7, 13, 20, 19, 8, 12, 17, 9, 11, 10, 2} );

        // Initialize tetrahedral point renumberings for GMSH->DUNE convention
        tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 3, 1, 2} );
        tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 7, 3, 4, 9, 1, 6, 8, 5, 2} );
        tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 11, 10, 3, 4, 17, 14, 5, 15, 1, 9, 18, 12, 16, 19, 6, 8, 13, 7, 2} );
        tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 15, 14, 13, 3, 4, 25, 27, 19, 5, 26, 20, 6, 21, 1, 12, 28, 29, 16, 22, 34, 31, 24, 32, 7, 11, 30, 17, 23, 33, 8, 10, 18, 9, 2} );
        tetrahedralInterpolatoryVertexGmsh2DuneMap.push_back( std::vector<int> {0, 19, 18, 17, 16, 3, 4, 34, 39, 36, 24, 5, 37, 38, 25, 6, 35, 26, 7, 27, 1, 15, 40, 43, 41, 20, 28, 52, 55, 46, 33, 53, 49, 30, 47, 8, 14, 45, 44, 21, 31, 54, 51, 32, 50, 9, 13, 42, 22, 29, 48, 10, 12, 23, 11, 2} );
   }

    int & totalVertex()   { return nVertexTotal_; }

    int & totalElement()   { return nElementTotal_; }

    int & totalInternalElement()   { return nInternalElementTotal_; }

    std::vector<int> & boundaryIdMap()   { return boundaryElement2PhysicalEntityIndex; }

    std::vector<int> & elementIndexMap() { return internalElement2PhysicalEntityIndex; }

    // This reads the GMSH format to parse the node and element structure
    void read (const std::string& f)
    {
    	std::string log_string;

        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: using file " + fileName);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: reading" + std::to_string(dim_) + "d curvilinear gmsh grid...");

        // open file name, we use C I/O
        // ***********************************************
        fileName = f;
        FILE* file = fopen(fileName.c_str(),"r");
        if (file==0)  { DUNE_THROW(Dune::IOError, "Could not open " << fileName); }

        // Reading MeshFormat Header
        // ***********************************************
        double version_number;        // process header
        int file_type, data_size;

        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$MeshFormat")!=0)   { DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line"); }
        fscanf(file, "%lg %d %d ", &version_number, &file_type, &data_size);
        if( (version_number < 2.0) || (version_number > 2.3) )  { DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files"); }
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: version " + std::to_string(version_number) + " Gmsh file detected");
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndMeshFormat")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndMeshFormat"); }

        // Reading Node data
        // ***********************************************
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Reading vertex header ------------------------");
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$Nodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $Nodes"); }
        fscanf(file, "%d\n", &nVertexTotal_);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: file contains " + std::to_string(nVertexTotal_) + " vertices");


        //==========================================================
        // VERTEX PASS 1: Put file pointer and skip all vertices
        //==========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Vertex-Pass 1: skip all vertices, since we need element data first ---");
        long section_vertex_offset = ftell(file);
        for (int iVertex = 0; iVertex < nVertexTotal_; iVertex++ )
        {
        	fgets(buf_, 512, file );
        	//Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, std::string(buf_));
        }
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$EndNodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndNodes"); }
        //==========================================================
        // VERTEX PASS 1: Finished
        //==========================================================


        // Reading Element Data
        // *************************************************
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Reading elements-header -----------------------");
        fscanf(file, "%s\n", buf_);
        if (strcmp(buf_,"$Elements")!=0)  { DUNE_THROW(Dune::IOError, "expected $Elements"); }
        fscanf(file, "%d\n", &nElementTotal_);
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: file contains " + std::to_string(nElementTotal_) + " elements");



        //=========================================================
        // ELEMENT PASS 1: Count the number of boundary segments
        //=========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Elements-Pass 1: counting elements on the boundary---");
        long fileOffsetElementSection = ftell(file);

        for (int i = 0; i < nElementTotal_; i++)
        {
            int id, gmshElementIndex;
            fscanf(file, "%d %d ", &id, &gmshElementIndex);

            // Ignore the rest of data on this line
            fgets(buf_, 512, file );

            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, std::to_string(id) + " " + std::to_string(gmshElementIndex) );

            // A boundary segment is defined here as any element with dimension less than world dimension
            GeometryType elemType          = gmshGeometryType(gmshElementIndex);
            int elemDim = elemType.dim();

            if (elemDim < dimWorld_ )     { nBoundaryElementTotal_++; }
            else                          { nInternalElementTotal_++; }
        }


        //=========================================================
        // ELEMENT PASS 2: Read all linear internal elements
        // Immediately partition them and distribute among processes
        //=========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Elements-Pass 2: reading internal elements for this process ---");
        fseek(file, fileOffsetElementSection, SEEK_SET);
        std::set< int > thisProcessElementIndexSet;

        readAndPartitionBaseElements(file, nElementTotal_, nInternalElementTotal_, thisProcessElementIndexSet);

        //==========================================================
        // ELEMENT PASS 3: Read all data associated with element.
        // Note all necessary vertex global indices into a map
        // Map all d-1 subentities of all elements to the element localID's.
        //    - Needed to find boundaries corresponding to each element.
        //==========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Elements-Pass 2: reading internal elements for this process ---");
        fseek(file, fileOffsetElementSection, SEEK_SET);
        std::set<int> thisProcessVertexIndexSet;
        std::vector< GmshElementData > internalElementVector;
        std::map< std::vector<int>, std::vector<int> > boundaryKey2LinkedElementSet;
        readInternalElements(file, nElementTotal_, nInternalElementTotal_, thisProcessElementIndexSet, internalElementVector, thisProcessVertexIndexSet, boundaryKey2LinkedElementSet);
        int nInternalElement = internalElementVector.size();


        //==========================================================
        // ELEMENT PASS 4: Read all data associated with boundary elements.
        // Only read boundary element if it is associated with this process
        //    - Note: Can not read internal and boundary elements at the same time
        //            because need d-1 subentity map from all internal elements first
        //==========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Elements-Pass 3: reading boundary elements for this process ---");
        fseek(file, fileOffsetElementSection, SEEK_SET);
        std::vector< GmshElementData > boundaryElementVector;
        std::vector< std::vector<int> > linkedElementLocalIndexSet;
        readBoundaryElements(file, nElementTotal_, boundaryElementVector, boundaryKey2LinkedElementSet, linkedElementLocalIndexSet);
        int nBoundaryElement = boundaryElementVector.size();



        //==========================================================
        // VERTEX PASS 2: Read the vertices
        // But only the ones that correspond to elements on this process
        //==========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Vertex-Pass 2: reading all vertices necessary for this process ---");
        fseek(file, section_vertex_offset, SEEK_SET);
        std::map<int, GlobalVector> vertexIndex2CoordinateMap;                // Only for testing purposes
        std::map<int, int> vertexGlobal2LocalIndexMap;
        int nVertex = readVertices(file, nVertexTotal_, vertexIndex2CoordinateMap, thisProcessVertexIndexSet, vertexGlobal2LocalIndexMap);


        //==========================================================
        // Final Step: Insert boundary segments and elements
        //==========================================================
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "----- Adding internal boundary elements to factory ---");

        addInternalElements(vertexGlobal2LocalIndexMap, vertexIndex2CoordinateMap, internalElementVector);
        addBoundaryElements(vertexGlobal2LocalIndexMap, vertexIndex2CoordinateMap, boundaryElementVector, linkedElementLocalIndexSet);


        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: total vertices          = " + std::to_string(nVertexTotal_)          + " of which on this process " + std::to_string(nVertex) );
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: total internal elements = " + std::to_string(nInternalElementTotal_) + " of which on this process " + std::to_string(nInternalElement) );
        Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, ":: total boundary elements = " + std::to_string(nBoundaryElementTotal_) + " of which on this process " + std::to_string(nBoundaryElement) );


        Dune::CollectiveCommunication<MPI_Comm> comm = mpihelper_.getCollectiveCommunication();
        int nElementParallelSum  = comm.sum(nInternalElement);
        int nBoundaryParallelSum = comm.sum(nBoundaryElement);

        if (rank_ == MASTER_RANK) {
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, " I am Master Process. Sum of internal elements = " +  std::to_string(nElementParallelSum) + ", sum of boundary elements = " +  std::to_string(nBoundaryParallelSum));
        }


        // TESTING SECTION - WRITES TEST ELEMENTS TO .VTK FILE
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (writeVtkFile_)
        {
        	vtkCurvWriter_.writeVTK("./curvreader_output_process_" + std::to_string(rank_) + ".vtk");
        	vtkCurvWriter_.writeParallelVTU("./curvreader_output");
            Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__,  "Curvilinear VTK Writer finished writing" );
        }

        // Close file
        fclose(file);
      }



    // Variables
    // *******************************************************************************

  protected:
    bool verbose_;         // If to show debug output
    bool processVerbose_;  // If to log using all processes or only master process
    bool writeVtkFile_;    // If to save mesh to VTK format after reading

    // Parallel Implementation
    MPIHelper &mpihelper_;
    static const int MASTER_RANK = 0;
    int rank_;
    int size_;

    // Grid Factory
    //Dune::GridFactory<GridType>& factory;
    Dune::CurvilinearGridFactory<GridType> & factory;

    // Reading file
    std::string fileName;
    char buf_[512];

    // Boundary element indexing
    bool insertBoundarySegment;
    std::vector<int> boundaryElement2PhysicalEntityIndex;
    std::vector<int> internalElement2PhysicalEntityIndex;

    // Total data about the mesh
    int nVertexTotal_ = 0;
    int nElementTotal_ = 0;
    int nInternalElementTotal_ = 0;
    int nBoundaryElementTotal_ = 0;

    // A map from GMSH -> Dune for indexing interpolatory points
    std::vector< std::vector< int > > triangularInterpolatoryVertexGmsh2DuneMap;
    std::vector< std::vector< int > > tetrahedralInterpolatoryVertexGmsh2DuneMap;

    // Testing capabilities for writing to VTK.
    CurvilinearVTKWriter<dimWorld_> vtkCurvWriter_;



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
    static void read (Dune::CurvilinearGridFactory<Grid> & factory,
                      const std::string& fileName,
                      MPIHelper &mpihelper,
                      std::vector<int>& boundaryElement2PhysicalEntityIndex,
                      std::vector<int>& internalElement2PhysicalEntityIndex,
                      int & nVertexTotal,
                      int & nElementTotal,
                      bool verbose,
                      bool processVerbose,
                      bool writeVTKFile,
                      bool insertBoundarySegment=true
                     )
    {
        // create parse object
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose, processVerbose, insertBoundarySegment, writeVTKFile, mpihelper);
        parser.read(fileName);

        boundaryElement2PhysicalEntityIndex.swap(parser.boundaryIdMap());
        internalElement2PhysicalEntityIndex.swap(parser.elementIndexMap());

        nVertexTotal = parser.totalVertex();
        nElementTotal = parser.totalInternalElement();
    }
  };


} // namespace Dune

#endif /** DUNE_CURVILINEARGMSHREADER_HH **/
