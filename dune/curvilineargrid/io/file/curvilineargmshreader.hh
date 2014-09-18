// curvilinear gmsh reader functionality - in this implementation of the curvilinear gmsh reader we ONLY read higher order triangles and tetrahedra


#ifndef DUNE_CURVILINEARGMSHREADER_HH
#define DUNE_CURVILINEARGMSHREADER_HH

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/geometry/polynomialinterpolation/curvilinearelementinterpolator.hh>

namespace Dune
{

  namespace {

    // arbitrary dimension, implementation is in specialization
    template< int dimension, int dimWorld = dimension >
    class CurvilinearGmshReaderQuadraticBoundarySegment
    {};

    // quadratic boundary segments in 1d
    /*
       Note the points

       (0)   (alpha)   (1)

       are mapped to the points in global coordinates

       p0 p2 p1

       alpha is determined automatically from the given points.
     */
    template< int dimWorld >
    struct CurvilinearGmshReaderQuadraticBoundarySegment< 2, dimWorld >
      : public Dune::BoundarySegment< 2, dimWorld >
    {
      typedef Dune::FieldVector< double, dimWorld > GlobalVector;

      CurvilinearGmshReaderQuadraticBoundarySegment ( const GlobalVector &p0_, const GlobalVector &p1_, const GlobalVector &p2_)
        : p0(p0_), p1(p1_), p2(p2_)
      {
        GlobalVector d1 = p1;
        d1 -= p0;
        GlobalVector d2 = p2;
        d2 -= p1;

        alpha=d1.two_norm()/(d1.two_norm()+d2.two_norm());
        if (alpha<1E-6 || alpha>1-1E-6)
          DUNE_THROW(Dune::IOError, "ration in quadratic boundary segment bad");
      }

      virtual GlobalVector operator() ( const Dune::FieldVector<double,1> &local ) const
      {
        GlobalVector y;
        y = 0.0;
        y.axpy((local[0]-alpha)*(local[0]-1.0)/alpha,p0);
        y.axpy(local[0]*(local[0]-1.0)/(alpha*(alpha-1.0)),p1);
        y.axpy(local[0]*(local[0]-alpha)/(1.0-alpha),p2);
        return y;
      }

    private:
      GlobalVector p0,p1,p2;
      double alpha;
    };


    // quadratic boundary segments in 2d
    /* numbering of points corresponding to gmsh:

       2

       5  4

       0  3  1

       Note: The vertices 3, 4, 5 are not necessarily at the edge midpoints but can
       be placed with parameters alpha, beta , gamma at the following positions
       in local coordinates:


       2 = (0,1)

       5 = (0,beta)   4 = (1-gamma/sqrt(2),gamma/sqrt(2))

       0 = (0,0)      3 = (alpha,0)                        1 = (1,0)

       The parameters alpha, beta, gamma are determined from the given vertices in
       global coordinates.
     */
    template<>
    class CurvilinearGmshReaderQuadraticBoundarySegment< 3, 3 >
      : public Dune::BoundarySegment< 3 >
    {
    public:
    CurvilinearGmshReaderQuadraticBoundarySegment (Dune::FieldVector<double,3> p0_, Dune::FieldVector<double,3> p1_,
                                          Dune::FieldVector<double,3> p2_, Dune::FieldVector<double,3> p3_,
                                          Dune::FieldVector<double,3> p4_, Dune::FieldVector<double,3> p5_)
        : p0(p0_), p1(p1_), p2(p2_), p3(p3_), p4(p4_), p5(p5_)
      {
        sqrt2 = sqrt(2.0);
        Dune::FieldVector<double,3> d1,d2;

        d1 = p3; d1 -= p0;
        d2 = p1; d2 -= p3;
        alpha=d1.two_norm()/(d1.two_norm()+d2.two_norm());
        if (alpha<1E-6 || alpha>1-1E-6)
          DUNE_THROW(Dune::IOError, "alpha in quadratic boundary segment bad");

        d1 = p5; d1 -= p0;
        d2 = p2; d2 -= p5;
        beta=d1.two_norm()/(d1.two_norm()+d2.two_norm());
        if (beta<1E-6 || beta>1-1E-6)
          DUNE_THROW(Dune::IOError, "beta in quadratic boundary segment bad");

        d1 = p4; d1 -= p1;
        d2 = p2; d2 -= p4;
        gamma=sqrt2*(d1.two_norm()/(d1.two_norm()+d2.two_norm()));
        if (gamma<1E-6 || gamma>1-1E-6)
          DUNE_THROW(Dune::IOError, "gamma in quadratic boundary segment bad");
      }

      virtual Dune::FieldVector<double,3> operator() (const Dune::FieldVector<double,2>& local) const
      {
        Dune::FieldVector<double,3> y;
        y = 0.0;
        y.axpy(phi0(local),p0);
        y.axpy(phi1(local),p1);
        y.axpy(phi2(local),p2);
        y.axpy(phi3(local),p3);
        y.axpy(phi4(local),p4);
        y.axpy(phi5(local),p5);
        return y;
      }

    private:
      // The six Lagrange basis function on the reference element
      // for the points given above

      double phi0 (const Dune::FieldVector<double,2>& local) const
      {
        return (alpha*beta-beta*local[0]-alpha*local[1])*(1-local[0]-local[1])/(alpha*beta);
      }
      double phi3 (const Dune::FieldVector<double,2>& local) const
      {
        return local[0]*(1-local[0]-local[1])/(alpha*(1-alpha));
      }
      double phi1 (const Dune::FieldVector<double,2>& local) const
      {
        return local[0]*(gamma*local[0]-(sqrt2-gamma-sqrt2*alpha)*local[1]-alpha*gamma)/(gamma*(1-alpha));
      }
      double phi5 (const Dune::FieldVector<double,2>& local) const
      {
        return local[1]*(1-local[0]-local[1])/(beta*(1-beta));
      }
      double phi4 (const Dune::FieldVector<double,2>& local) const
      {
        return local[0]*local[1]/((1-gamma/sqrt2)*gamma/sqrt2);
      }
      double phi2 (const Dune::FieldVector<double,2>& local) const
      {
        return local[1]*(beta*(1-gamma/sqrt2)-local[0]*(beta-gamma/sqrt2)-local[1]*(1-gamma/sqrt2))/((1-gamma/sqrt2)*(beta-1));
      }

      Dune::FieldVector<double,3> p0,p1,p2,p3,p4,p5;
      double alpha,beta,gamma,sqrt2;
    };

  }   // end empty namespace















  //! dimension independent parts for CurvilinearGmshReaderParser
  template<typename GridType>
  class CurvilinearGmshReaderParser
  {
  protected:
    // private data
    Dune::GridFactory<GridType>& factory;
    bool verbose = false;
    bool insert_boundary_segments;
    unsigned int number_of_real_vertices;
    int boundary_element_count;
    int element_count;
    // read buffer
    char buf[512];
    std::string fileName;
    // exported data
    std::vector<int> boundary_id_to_physical_entity;
    std::vector<int> element_index_to_physical_entity;


    std::vector<std::vector<int>> triangularPointRenumberings_;
    std::vector<std::vector<int>> tetrahedralPointRenumberings_;


    std::map<int, unsigned int> nodeofelement; /** \brief map stores the number of vertex indices for a specific element type **/

    // static data
    static const int dim = GridType::dimension;
    static const int dimWorld = GridType::dimensionworld;
    static_assert( (dimWorld <= 3), "GmshReader requires dimWorld <= 3." );

    // typedefs
    typedef FieldVector< double, dimWorld > GlobalVector;



    // TODO: Possibly merge this element name enum with internal DUNE one.
    enum{
         GMSH_NODE            = 1,
         GMSH_EDGE            = 2,
         GMSH_TRIANGLE        = 3,
         GMSH_QUADRANGLE      = 4,
         GMSH_TETRAHEDRON     = 5,
         GMSH_HEXAHEDRON      = 6,
         GMSH_PRISM           = 7,
         GMSH_PYRAMID         = 8
    };

    // Returns the dimension corresponding to a given element type name
    // TODO: Throw error if wrong element name or order
    int GMSH_ElementDim(int elemName) {

        switch (elemName) {
            case GMSH_NODE             : return 0;   break;
            case GMSH_EDGE             : return 1;   break;
            case GMSH_TRIANGLE         : return 2;   break;
            case GMSH_QUADRANGLE       : return 2;   break;
            case GMSH_TETRAHEDRON      : return 3;   break;
            case GMSH_HEXAHEDRON       : return 3;   break;
            case GMSH_PRISM            : return 3;   break;
            case GMSH_PYRAMID          : return 3;   break;
        }

        return 0;
    }

    // Returns the corner number corresponding to a given element type name
    // TODO: Throw error if wrong element name
    int GMSH_ElementCorners(int elemName) {

        switch (elemName) {
            case GMSH_NODE             : return 1;   break;
            case GMSH_EDGE             : return 2;   break;
            case GMSH_TRIANGLE         : return 3;   break;
            case GMSH_QUADRANGLE       : return 4;   break;
            case GMSH_TETRAHEDRON      : return 4;   break;
            case GMSH_HEXAHEDRON       : return 8;   break;
            case GMSH_PRISM            : return 6;   break;
            case GMSH_PYRAMID          : return 5;   break;
        }

        return 0;
    }

    // TODO: Throw error if wrong element name or order
    int GMSH_ElementDofNumber(int elemName, int order)
    {
        switch (elemName) {
            case GMSH_EDGE             : return order + 1;                                   break;    // {2,3,4,5,6}
            case GMSH_TRIANGLE         : return (order + 1) * (order + 2) / 2;               break;    // {3,6,10,15,21}
            case GMSH_TETRAHEDRON      : return (order + 1)*(order + 2)*(order + 3) / 6  ;   break;    // {4,10,20,35,56}
        }
    	return 0;
    }

    // Constructs a DUNE geometry type based on GMSH element index
    Dune::GeometryType GMSH_GeometryType(int gmsh_index)
    {
    	Dune::GeometryType rez;
    	int gi = gmsh_index;

    	     if (gi == 15)																									   { rez.makeVertex(); }
    	else if ((gi == 1) || (gi == 8) || (gi == 26) || (gi == 27) || (gi == 28))											   { rez.makeLine(); }
    	else if ((gi == 2) || (gi == 9) || (gi == 20) || (gi == 21) || (gi == 22) || (gi == 23) || (gi == 24) || (gi == 25))   { rez.makeTriangle(); }
    	else if ((gi == 3) || (gi == 10) || (gi == 16))											 							   { rez.makeQuadrilateral(); }
    	else if ((gi == 4) || (gi == 11) || (gi == 29) || (gi == 30) || (gi == 31))											   { rez.makeTetrahedron(); }
    	else if ((gi == 5) || (gi == 12) || (gi == 17) || (gi == 92) || (gi == 93))											   { rez.makeHexahedron(); }
    	else if ((gi == 6) || (gi == 13) || (gi == 18))											 							   { rez.makePrism(); }
    	else if ((gi == 7) || (gi == 14) || (gi == 19))											 							   { rez.makePyramid(); }
    	else  { DUNE_THROW(Dune::IOError, "Unexpected geometry type");  }

    	return rez;
    }

    // Returns the type name of the element given its GMSH_index
    int GMSH_ElementName(int gmsh_index)
    {
    	// Hexahedra of high dimension have funny index
    	if ((gmsh_index == 92) || (gmsh_index == 93)) { return GMSH_HEXAHEDRON; }

    	// Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    	const int elemName[32]        = { -1,
    	            GMSH_EDGE, GMSH_TRIANGLE, GMSH_QUADRANGLE, GMSH_TETRAHEDRON, GMSH_HEXAHEDRON, GMSH_PRISM, GMSH_PYRAMID, GMSH_EDGE, GMSH_TRIANGLE, GMSH_QUADRANGLE,
    	            GMSH_TETRAHEDRON, GMSH_HEXAHEDRON, GMSH_PRISM, GMSH_PYRAMID, GMSH_NODE, GMSH_QUADRANGLE, GMSH_HEXAHEDRON, GMSH_PRISM, GMSH_PYRAMID, GMSH_TRIANGLE,
    	            GMSH_TRIANGLE, GMSH_TRIANGLE, GMSH_TRIANGLE, GMSH_TRIANGLE, GMSH_TRIANGLE, GMSH_EDGE, GMSH_EDGE, GMSH_EDGE, GMSH_TETRAHEDRON, GMSH_TETRAHEDRON, GMSH_TETRAHEDRON };

    	return elemName[gmsh_index];
    }

    // Returns the type name of the element given its GMSH_index
    int GMSH_ElementOrder(int gmsh_index)
    {
    	// Hexahedra of high dimension have funny index
    	if (gmsh_index == 92) { return 3; }
    	if (gmsh_index == 93) { return 4; }

    	// Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    	const int elemOrder[32]          = {-1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 3, 4, 5, 3, 4, 5};

    	return elemOrder[gmsh_index];
    }

    bool GMSH_ElementIsIncomplete(int gmsh_index)
    {
    	// Both high-order hexahedrons are complete
    	if ((gmsh_index == 92) || (gmsh_index == 93)) { return 0; }

    	// Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    	const bool elemIncomplete[32]          = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};

    	return elemIncomplete[gmsh_index];
    }

    // Returns the number of degrees of freedom of the element given its GMSH_index
    // note: This info can not simply be obtained from referenceElement, because some of the elements in GMSH have incomplete order, so less DoF than expected
    int GMSH_ElementDofNumber(int gmsh_index)
    {
    	// Array copy-pasted from GMSH Brute-Force because it does not seem to have any pattern :)
    	const int nDofs[32]              = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56};

    	return nDofs[gmsh_index];
    }

    // Returns the total number of DoF associated with all subentities of a given dimension for this element, subtracting the ones that come from the corners
    int GMSH_ElementSubentityExtraDim(int gmsh_index, int dim)
    {
    	const int nDofsExtraEdges[32] = {-1, 0, 0, 0, 0, 0, 0, 0, 1, 3, 4, 6, 12, 9, 8, 0, 4, 12, 9, 8, 6, 6, 9, 9, 12, 12, 2, 3, 4, 12, 18, 24};
    	const int nDofsExtraFaces[32] = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 6, 3, 1, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 6, 0, 0, 0, 4, 12, 24};
    	const int nDofsExtraElems[32] = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4};

    	switch(dim)
    	{
    	case 1: return nDofsExtraEdges[gmsh_index];  break;
    	case 2: return nDofsExtraFaces[gmsh_index];  break;
    	case 3: return nDofsExtraElems[gmsh_index];  break;
    	}

    	return -1;
    }

    // correct differences between gmsh and Dune in the local vertex numbering
    // NOTE: THIS METHOD DOES NOT WORK WITH INCOMPLETE ORDER GMSH ELEMENTS AT THE MOMENT
    void Gmsh2DuneElementDofNumbering(std::vector<int> &elementDofs, int thisElmName, int thisElmOrder) {
        int thisElmDofNo = GMSH_ElementDofNumber(thisElmName, thisElmOrder);
        std::vector<int> tmp;

        switch (thisElmName)
        {
        case GMSH_TRIANGLE :
            for (int i = 0; i < thisElmDofNo; i++) { tmp.push_back(elementDofs[triangularPointRenumberings_[thisElmOrder - 1][i]]); }
            for (int i = 0; i < thisElmDofNo; i++) { elementDofs[i] = tmp[i]; }

            break;
        case GMSH_TETRAHEDRON :
            for (int i = 0; i < thisElmDofNo; i++) { tmp.push_back(elementDofs[tetrahedralPointRenumberings_[thisElmOrder - 1][i]]); }
            for (int i = 0; i < thisElmDofNo; i++) { elementDofs[i] = tmp[i]; }

            break;
        }
    }

  public:

    CurvilinearGmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      factory(_factory), verbose(v), insert_boundary_segments(i)
   {
        element_count = 0;
        boundary_element_count = 0;
        number_of_real_vertices = 0;


        // Initialize triangular point renumberings for GMSH->DUNE convention
        triangularPointRenumberings_.push_back( std::vector<int> {0, 1, 2} );
        triangularPointRenumberings_.push_back( std::vector<int> {0, 3, 1, 5, 4, 2} );
        triangularPointRenumberings_.push_back( std::vector<int> {0, 3, 4, 1, 8, 9, 5, 7, 6, 2} );
        triangularPointRenumberings_.push_back( std::vector<int> {0, 3, 4, 5, 1, 11, 12, 13, 6, 10, 14, 7, 9, 8, 2} );
        triangularPointRenumberings_.push_back( std::vector<int> {0, 3, 4, 5, 6, 1, 14, 15, 18, 16, 7, 13, 20, 19, 8, 12, 17, 9, 11, 10, 2} );

        // Initialize tetrahedral point renumberings for GMSH->DUNE convention
        tetrahedralPointRenumberings_.push_back( std::vector<int> {0, 3, 1, 2} );
        tetrahedralPointRenumberings_.push_back( std::vector<int> {0, 7, 3, 4, 9, 1, 6, 8, 5, 2} );
        tetrahedralPointRenumberings_.push_back( std::vector<int> {0, 11, 10, 3, 4, 17, 14, 5, 15, 1, 9, 18, 12, 16, 19, 6, 8, 13, 7, 2} );
        tetrahedralPointRenumberings_.push_back( std::vector<int> {0, 15, 14, 13, 3, 4, 25, 27, 19, 5, 26, 20, 6, 21, 1, 12, 28, 29, 16, 22, 34, 31, 24, 32, 7, 11, 30, 17, 23, 33, 8, 10, 18, 9, 2} );
        tetrahedralPointRenumberings_.push_back( std::vector<int> {0, 19, 18, 17, 16, 3, 4, 34, 39, 36, 24, 5, 37, 38, 25, 6, 35, 26, 7, 27, 1, 15, 40, 43, 41, 20, 28, 52, 55, 46, 33, 53, 49, 30, 47, 8, 14, 45, 44, 21, 31, 54, 51, 32, 50, 9, 13, 42, 22, 29, 48, 10, 12, 23, 11, 2} );
   }

    std::vector<int> & boundaryIdMap()   { return boundary_id_to_physical_entity; }

    std::vector<int> & elementIndexMap() { return element_index_to_physical_entity; }

    // This reads the GMSH format to parse the node and element structure
    // TODO: May want to initialise the GMSH info arrays outside the loop, perhaps even outside read function
    void read (const std::string& f)
    {
      verbose = false;

      if (verbose) std::cout << "using file " << fileName.c_str() << std::endl;
      if (verbose) std::cout << "::: reading" << dim << "d curvilinear gmsh grid..." << std::endl;

      // open file name, we use C I/O
      fileName = f;
      FILE* file = fopen(fileName.c_str(),"r");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);



      //=========================================
      // Header: Read vertices into vector
      //         Check vertices that are needed
      //=========================================

      number_of_real_vertices = 0;
      boundary_element_count = 0;
      element_count = 0;

      // process header
      double version_number;
      int file_type, data_size;

      // Reading MeshFormat Header
      // ***********************************************
      fscanf(file, "%s\n", buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      fscanf(file, "%lg %d %d ", &version_number, &file_type, &data_size);
      if( (version_number < 2.0) || (version_number > 2.3) )
        DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files");
      if (verbose) std::cout << "version " << version_number << " Gmsh file detected" << std::endl;
      fscanf(file, "%s\n", buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // Reading Node data
      // ***********************************************
      int number_of_nodes;

      // Number of nodes
      fscanf(file, "%s\n", buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      fscanf(file, "%d\n", &number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;

      // Actual nodes
      std::vector< GlobalVector > nodes( number_of_nodes+1 );       // store positions
      {
        int id;
        double x[3];
        for( int i = 1; i <= number_of_nodes; ++i )
        {
          fscanf(file, "%d %lg %lg %lg\n", &id, &x[0], &x[1], &x[2] );
          if( id != i )
            DUNE_THROW( Dune::IOError, "Expected id " << i << "(got id " << id << "." );

            nodes[i][0] = x[0];
            nodes[i][1] = x[1];
            nodes[i][2] = x[2];
        }

        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$EndNodes")!=0)
          DUNE_THROW(Dune::IOError, "expected $EndNodes");
      }


      // Reading Element Data
      // *************************************************
      int number_of_elements;

      // Number of elements
      fscanf(file, "%s\n", buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      fscanf(file, "%d\n", &number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl;

      long section_element_offset = ftell(file);
      std::map<int,unsigned int> renumber;
      std::vector<int> mesh_partitions;

      // For testing we will create a vector of all possible triangles
      std::vector<GlobalVector> vtk_points;
      std::vector<std::vector<int>> vtk_triangles;

      // For testing we may choose to create a set of edges
      std::vector<std::vector<int>> vtk_edges;


      /** \brief Handle one line in the gmsh file, i.e. handle an element
       *
       * We would like to read:
       *  1) Physical tag associated with each element, associated to, for example, material properties of the element
       *  2) Type of the element, which would also give us its curvilinear order
       *  3) All nodes associated with the element, including the
       *
       * Possible tags
       *   k == 0: physical entity - very important tag
       *   k == 1: elementary entity (not used here)
       *   if version_number < 2.2:
       *     k == 2: mesh partition 0
       *   else
       *     k == 3: number of mesh partitions
       *     k => 3: mesh partition k-3
       *
       * **/


      for (int i=1; i<=number_of_elements; i++)
      {
          int id, gmsh_element_index, number_of_tags;

          if (verbose) std::cout << "Started reading elements" << std::endl;
          fscanf(file, "%d %d %d ", &id, &gmsh_element_index, &number_of_tags);

          if (verbose) std::cout << "element " << i << " has " << number_of_tags << " tags" << std::endl;

          // Read all tags
          std::vector<int> elm_tags;
          for (int k = 0; k < number_of_tags; k++)
          {
              int tmp_tag;
              fscanf(file, "%d ", &tmp_tag);
              elm_tags.push_back(tmp_tag);
          }

          // Possible functionality for mesh partitioning
          if (number_of_tags > 3)
          {
              mesh_partitions.resize(elm_tags[2]);
              for (int j = 1; j <= elm_tags[2]; j++)  {    mesh_partitions[j] = elm_tags[2 + j];  }
          }

          // Obtain all necessary info not to use gmsh_element_index in the following steps
          Dune::GeometryType thisElmType = GMSH_GeometryType(gmsh_element_index);
          int thisElmName                = GMSH_ElementName(gmsh_element_index);
          int thisElmOrder               = GMSH_ElementOrder(gmsh_element_index);
          int thisElmDofNo               = GMSH_ElementDofNumber(gmsh_element_index);

          int thisElmDim                 = GMSH_ElementDim(thisElmName);
          int thisElmCornerNo            = GMSH_ElementCorners(thisElmName);
          int thisElmPhysicalEntity      = elm_tags[0];


          if (verbose) std::cout << "element " << i << " has dimension " << thisElmDim << " and vertex number " << thisElmCornerNo << " and physical entity number " << thisElmPhysicalEntity << std::endl;



          // Testing if the current element type can be handled by DUNE
          // *****************************************************************
          bool isAllowedElement = true;

          // Check if element is polynomial-complete (ask GMSH what that even means I dont know :) )
          isAllowedElement &= !GMSH_ElementIsIncomplete(gmsh_element_index);

          // in 2D only allow lines and triangles
          isAllowedElement &= (dim != 2) || ((thisElmName == GMSH_EDGE)||(thisElmName == GMSH_TRIANGLE));

          // in 3D only allow triangles and tetrahedrons
          isAllowedElement &= (dim != 3) || ((thisElmName == GMSH_TRIANGLE)||(thisElmName == GMSH_TETRAHEDRON));

          // test whether we support the element type at the moment
          if (!isAllowedElement)
              DUNE_THROW(Dune::IOError, "GMSH Reader: Have read an element of unexpected type " << gmsh_element_index);
          if (verbose) std::cout << "element " << i << " can be treated by Dune grid " << std::endl;



          // Reading DoF's
          // *************************************************
          std::vector<int> elementDofs;

          for (int iDof = 0; iDof < thisElmDofNo; iDof++) {
              int tmp;
              fscanf(file, "%d", &tmp);
              if (verbose) std::cout << "  --- have read DoF " << tmp << std::endl;
              elementDofs.push_back(tmp);
          }
          fscanf(file, "\n");

            // correct differences between gmsh and Dune in the local vertex numbering
          Gmsh2DuneElementDofNumbering(elementDofs, thisElmName, thisElmOrder);


          // TESTING SECTION FOR TETRAHEDRA
          // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          switch (thisElmName)
          {
          case GMSH_TRIANGLE    : addTestTriangleFromSurface(elementDofs, thisElmType, thisElmOrder, thisElmDofNo, nodes, vtk_points, vtk_triangles);         break;
          //case GMSH_TETRAHEDRON : addTestTrianglesFromTetrahedron_NoInterp(elementDofs, thisElmOrder, thisElmDofNo, nodes, vtk_points, vtk_triangles);        break;
          case GMSH_TETRAHEDRON : addTestTrianglesFromTetrahedron(elementDofs, thisElmType, thisElmOrder, thisElmDofNo, nodes, vtk_points, vtk_triangles);    break;
          //case GMSH_TETRAHEDRON : addTestEdgesFromTetrahedron(elementDofs, thisElmType, thisElmOrder, thisElmDofNo, nodes, vtk_points, vtk_edges);            break;
          }
          if (verbose) std::cout << "element " << i << " has been added to the VTK triangles  " << std::endl;
          // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



          // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          // AT THE MOMENT DUNE GRID CAN NOT HANDLE CURVILINEAR ELEMENTS
          // THE FACTORY WILL ONLY BE INITIALIZED WITH NODES THAT ARE VERTICES
          // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          // insert each vertex if it hasn't been inserted already
          for (int i=0; i < thisElmCornerNo; i++)
          {
              if (renumber.find(elementDofs[i])==renumber.end())
              {
                  renumber[elementDofs[i]] = number_of_real_vertices++;
                  if (verbose) std::cout << " --- adding " << i + 1 << "th vertex with number " << elementDofs[i] << " to factory" << std::endl;
                  factory.insertVertex(nodes[elementDofs[i]]);
              }
          }

          // renumber corners to account for the explicitly given vertex
          // numbering in the file
          std::vector<unsigned int> corners(thisElmCornerNo);
          for (int i=0; i < thisElmCornerNo; i++)    { corners[i] = renumber[elementDofs[i]]; }


          //Insert boundary segments and elements
          //****************************************************
          addElementToFactory(elementDofs, corners, nodes, thisElmName, thisElmDim, thisElmOrder, thisElmPhysicalEntity);
          if (verbose) std::cout << "element " << i << " has been added to the Geometry Factory " << std::endl;
        }


      // TESTING SECTION - WRITES TEST ELEMENTS TO .    VTK FILE
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      VTK_Write(vtk_points, vtk_edges, vtk_triangles);
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      if (verbose) std::cout << "number of boundary elements = " << boundary_element_count << std::endl;
      if (verbose) std::cout << "number of elements = " << element_count << std::endl;
      fscanf(file, "%s\n", buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");
      boundary_id_to_physical_entity.resize(boundary_element_count);
      element_index_to_physical_entity.resize(element_count);

      fclose(file);
    }

    // generic-case: This is not supposed to be used at runtime.
    template <class E, class V, class V2>
    void boundarysegment_insert(
      const V& nodes,
      const E& elementDofs,
      const V2& vertices
      )
    {
      DUNE_THROW(Dune::IOError, "tried to create a 3D boundary segment in a non-3D Grid");
    }

    // 3d-case:
    template <class E, class V>
    void boundarysegment_insert(
      const std::vector<FieldVector<double, 3> >& nodes,
      const E& elementDofs,
      const V& vertices
      )
    {
      array<FieldVector<double,dimWorld>, 6> v;
      for (int i=0; i<6; i++)
        for (int j=0; j<dimWorld; j++)
          v[i][j] = nodes[elementDofs[i]][j];

      BoundarySegment<dim,dimWorld>* newBoundarySegment
        = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 3, 3 >( v[0], v[1], v[2],
                                                                                           v[3], v[4], v[5] );

      factory.insertBoundarySegment( vertices,
                                     shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment) );
    }


    // Test Function - Takes a triangle with all its interpolatory points, splits it into visualisation triangles
    // and adds all of them to vtk_triangle
    void addTestTriangleFromSurface(
            std::vector<int> & elementDofs,
            Dune::GeometryType & thisElmType,
            int thisElmOrder,
            int thisElmDofNo,
            std::vector< GlobalVector > & nodes,
            std::vector<GlobalVector> & vtk_points,
            std::vector<std::vector<int>> & vtk_triangles
            )
    {
        // Get nodes associated with this element
        std::vector<GlobalVector> elementNodes;
        for (int iDof = 0; iDof < thisElmDofNo; iDof++) { elementNodes.push_back(nodes[elementDofs[iDof]]); }

        CurvilinearElementInterpolator<double, 2, 3> thisTriangleInt(thisElmType, elementNodes, thisElmOrder);

          // Shrink triangle
        // 0) Find coordinates of the corners
        int corner1 = 0;
        int corner2 = thisElmOrder;
        int corner3 = thisElmDofNo - 1;

          // 1) Compute Center of Mass
          double CoM[3];
          double shrink_magnitude = 0.2;                // 0.0 - no shrinking, 0.99 - very small element (must be < 1)
          CoM[0] = (elementNodes[corner1][0] + elementNodes[corner2][0] + elementNodes[corner3][0]) / 3.0;
          CoM[1] = (elementNodes[corner1][1] + elementNodes[corner2][1] + elementNodes[corner3][1]) / 3.0;
          CoM[2] = (elementNodes[corner1][2] + elementNodes[corner2][2] + elementNodes[corner3][2]) / 3.0;


          //std::cout << "Element with init dof " << elementDofs[0] << std::endl;


          // 2) Construct sampling coordinates
          std::map<std::vector<int>, int> parametricToIndex;
          int nDiscretizationPoints = 16;
          int nIntervals = nDiscretizationPoints - 1;

          // Step 1: Generate all points on this face. Generic triangle looks like this
          //
          //  (i,j)=(n,0) -> *
          //                 **
          //                 ***
          //                 ****
          //  (i,j)=(0,0) -> ****** <- (i,j) = (0,n)
          //
          for (int i = 0; i < nDiscretizationPoints; i++) {
              for (int j = 0; j < nDiscretizationPoints - i; j++) {
                  // Map local triangle parametric coordinates to local tetrahedral parametric coordinates
                  std::vector<int> parUV {j, i};

                  // Check if this point was already calculated
                  if (parametricToIndex.find(parUV) == parametricToIndex.end())
                  {
                      // Calculate the real coordinates of the point
                      FieldVector<double, 2> localVector;
                      localVector[0] = double(parUV[0]) / nIntervals;
                      localVector[1] = double(parUV[1]) / nIntervals;
                      GlobalVector tmpPoint = thisTriangleInt.realCoordinate(localVector);

                      // Shrink each coordinate of each point towards the CoM
                      tmpPoint[0] += shrink_magnitude * (CoM[0] - tmpPoint[0]);
                      tmpPoint[1] += shrink_magnitude * (CoM[1] - tmpPoint[1]);
                      tmpPoint[2] += shrink_magnitude * (CoM[2] - tmpPoint[2]);

                      // Expand the surface a little bit so that it does not interlay with the volume
                      tmpPoint[0] *= 1.2;
                      tmpPoint[1] *= 1.2;
                      tmpPoint[2] *= 1.2;

                      // Add coordinates to the coordinates array
                      vtk_points.push_back(tmpPoint);

                      // Add point to the point map
                      parametricToIndex[parUV] = vtk_points.size() - 1;

                      //std::cout << "* coords " << parUV[0] << ", " << parUV[1] << ", in the map cooresponds to " << parametricToIndex[parUV] << std::endl;
                  }
              }
          }




          // Step 2: Split this face into tiny triangles and add them to the triangle array
          for (int i = 0; i < nDiscretizationPoints - 1; i++) {
              for (int j = 0; j < nDiscretizationPoints - 1 - i; j++) {
                  // Construct two triangles (123) and (234) from 4 points, where 1 is the index point
                  //   3--4
                  //   |\ |
                  //   | \|
                  //   1--2
                  // First triangle always exists, because we iterate over (i,j) in such a way that there is 1 free point at the edge
                  // Second triangle we construct only if point 4 is still in the triangle

                  // Construct triangle (123)
                  std::vector<int> parUV_1 {j,    i};
                  std::vector<int> parUV_2 {j+1,    i};
                  std::vector<int> parUV_3 {j,    i+1};

                  std::vector<int> triangle123;
                  triangle123.push_back(parametricToIndex[parUV_1]);
                  triangle123.push_back(parametricToIndex[parUV_2]);
                  triangle123.push_back(parametricToIndex[parUV_3]);

                  //std::cout << "- coords " << parUV_1[0] << ", " << parUV_1[1] << ", in the map cooresponds to " << parametricToIndex[parUV_1] << std::endl;
                  //std::cout << "- coords " << parUV_2[0] << ", " << parUV_2[1] << ", in the map cooresponds to " << parametricToIndex[parUV_2] << std::endl;
                  //std::cout << "- coords " << parUV_3[0] << ", " << parUV_3[1] << ", in the map cooresponds to " << parametricToIndex[parUV_3] << std::endl;

                  // Add triangle (123)
                  vtk_triangles.push_back(triangle123);

                  // Check if point 4 is still in the triangle
                  if (i + j + 2 < nDiscretizationPoints)
                  {
                      std::vector<int> parUV_4 {j+1,    i+1};

                      std::vector<int> triangle234;
                      triangle234.push_back(parametricToIndex[parUV_2]);
                      triangle234.push_back(parametricToIndex[parUV_3]);
                      triangle234.push_back(parametricToIndex[parUV_4]);

                      // Add triangle (234)
                      vtk_triangles.push_back(triangle234);
                  }
              }
          }
    }

    // Test Function
    // Takes Tetrahedron with all its interpolatory points, shrinks it by 5%, splits all DoF points into triangles, adds all of them to vtk_triangles
    void addTestTrianglesFromTetrahedron_NoInterp(
            std::vector<int> & elementDofs,
            int thisElmOrder,
            int thisElmDofNo,
            std::vector< GlobalVector > & nodes,
            std::vector<GlobalVector> & vtk_points,
            std::vector<std::vector<int>> & vtk_triangles

            )
    {

        // Get nodes associated with this element
        std::vector<GlobalVector> elementNodes;
        for (int iDof = 0; iDof < thisElmDofNo; iDof++) { elementNodes.push_back(nodes[elementDofs[iDof]]); }

          // Shrink tet
        // 0) Find coordinates of the corners
        int corner1 = 0;
        int corner2 = thisElmOrder;
        int corner3 = thisElmOrder*(thisElmOrder + 3) / 2;
        int corner4 = thisElmDofNo - 1;


          // 1) Compute Center of Mass
          double CoM[3];
          double shrink_magnitude = 0.2;                // 0.0 - no shrinking, 0.99 - very small element (must be < 1)
          CoM[0] = (elementNodes[corner1][0] + elementNodes[corner2][0] + elementNodes[corner3][0] + elementNodes[corner4][0]) / 4.0;
          CoM[1] = (elementNodes[corner1][1] + elementNodes[corner2][1] + elementNodes[corner3][1] + elementNodes[corner4][1]) / 4.0;
          CoM[2] = (elementNodes[corner1][2] + elementNodes[corner2][2] + elementNodes[corner3][2] + elementNodes[corner4][2]) / 4.0;

          // 3) Construct sampling coordinates
          std::map<std::vector<int>, int> parametricToIndex;
          int nDiscretizationPoints = thisElmOrder + 1;        // Important - we only use already existing points for triangles here
          int nIntervals = nDiscretizationPoints - 1;

          int iNodes = 0;

          // Generate all points
          for (int i = 0; i < nDiscretizationPoints; i++) {
              for (int j = 0; j < nDiscretizationPoints - i; j++) {
                  for (int k = 0; k < nDiscretizationPoints - i - j; k++) {
                      // Iterate over all DoF and the corresponding coordinates in the reference tetrahedron
                      std::vector<int> parUVW {k, j, i};

                      // If the point is on the surface, add it
                      if ((i == 0)||(j==0)||(k==0)||(i+j+k == nDiscretizationPoints - 1))
                      {
                          GlobalVector tmpPoint = nodes[elementDofs[iNodes]];

                          // Shrink each coordinate of each point towards the CoM
                          tmpPoint[0] += shrink_magnitude * (CoM[0] - tmpPoint[0]);
                          tmpPoint[1] += shrink_magnitude * (CoM[1] - tmpPoint[1]);
                          tmpPoint[2] += shrink_magnitude * (CoM[2] - tmpPoint[2]);


                          // Add coordinates to the coordinates array
                          vtk_points.push_back(tmpPoint);

                          // Add point to the point map
                          parametricToIndex[parUVW] = vtk_points.size() - 1;
                      }
                      iNodes++;
                  }
              }
          }

          for (int iFace = 0; iFace < 4; iFace++)
          {
              // Step 2: Split this face into tiny triangles and add them to the triangle array
              for (int i = 0; i < nDiscretizationPoints - 1; i++) {
                  for (int j = 0; j < nDiscretizationPoints - 1 - i; j++) {
                      // Construct two triangles (123) and (234) from 4 points, where 1 is the index point
                      //   3--4
                      //   |\ |
                      //   | \|
                      //   1--2
                      // First triangle always exists, because we iterate over (i,j) in such a way that there is 1 free point at the edge
                      // Second triangle we construct only if point 4 is still in the triangle

                      // Construct triangle (123)
                      std::vector<int> parUVW_1 = referenceTetrahedronFaceParameterization(iFace, i, j, nIntervals);
                      std::vector<int> parUVW_2 = referenceTetrahedronFaceParameterization(iFace, i, j+1, nIntervals);
                      std::vector<int> parUVW_3 = referenceTetrahedronFaceParameterization(iFace, i+1, j, nIntervals);

                      std::vector<int> triangle123;
                      triangle123.push_back(parametricToIndex[parUVW_1]);
                      triangle123.push_back(parametricToIndex[parUVW_2]);
                      triangle123.push_back(parametricToIndex[parUVW_3]);

                      // Add triangle (123)
                      vtk_triangles.push_back(triangle123);

                      // Check if point 4 is still in the triangle
                      if (i + j + 2 < nDiscretizationPoints)
                      {
                          std::vector<int> parUVW_4 = referenceTetrahedronFaceParameterization(iFace, i+1, j+1, nIntervals);

                            std::vector<int> triangle234;
                            triangle234.push_back(parametricToIndex[parUVW_2]);
                            triangle234.push_back(parametricToIndex[parUVW_3]);
                            triangle234.push_back(parametricToIndex[parUVW_4]);

                            // Add triangle (234)
                            vtk_triangles.push_back(triangle234);
                      }
                  }
              }
          }
    }

    // Test Function
    // Takes Tetrahedron - shrinks it by 5%, interpolates over it, iterates over all surfaces,
    // for each surface uniformly generates points, splits them into triangles, adds all of them to vtk_triangles
    void addTestTrianglesFromTetrahedron(
            std::vector<int> & elementDofs,
            Dune::GeometryType & thisElmType,
            int thisElmOrder,
            int thisElmDofNo,
            std::vector< GlobalVector > & nodes,
            std::vector<GlobalVector> & vtk_points,
            std::vector<std::vector<int>> & vtk_triangles

            )
    {

        // Get nodes associated with this element
        std::vector<GlobalVector> elementNodes;
        for (int iDof = 0; iDof < thisElmDofNo; iDof++) { elementNodes.push_back(nodes[elementDofs[iDof]]); }

        // Construct Tetrahedral Interpolator used to map parametric coordinates to real interpolated ones
        CurvilinearElementInterpolator<double, 3, 3> thisTetrahedronInt(thisElmType, elementNodes, thisElmOrder);

          // Shrink tet
        // 0) Find coordinates of the corners
        int corner1 = 0;
        int corner2 = thisElmOrder;
        int corner3 = thisElmOrder*(thisElmOrder + 3) / 2;
        int corner4 = thisElmDofNo - 1;


          // 1) Compute Center of Mass
          double CoM[3];
          double shrink_magnitude = 0.2;                // 0.0 - no shrinking, 0.99 - very small element (must be < 1)
          CoM[0] = (elementNodes[corner1][0] + elementNodes[corner2][0] + elementNodes[corner3][0] + elementNodes[corner4][0]) / 4.0;
          CoM[1] = (elementNodes[corner1][1] + elementNodes[corner2][1] + elementNodes[corner3][1] + elementNodes[corner4][1]) / 4.0;
          CoM[2] = (elementNodes[corner1][2] + elementNodes[corner2][2] + elementNodes[corner3][2] + elementNodes[corner4][2]) / 4.0;

          // 3) Construct sampling coordinates
          std::map<std::vector<int>, int> parametricToIndex;
          //int nDiscretizationPoints = 2;
          int nDiscretizationPoints = 16;
          int nIntervals = nDiscretizationPoints - 1;

          for (int iFace = 0; iFace < 4; iFace++)
          {
              // Step 1: Generate all points on this face. Generic triangle looks like this
              //
              //  (i,j)=(n,0) -> *
              //                 **
              //                 ***
              //                 ****
              //  (i,j)=(0,0) -> ****** <- (i,j) = (0,n)
              //
              for (int i = 0; i < nDiscretizationPoints; i++) {
                  for (int j = 0; j < nDiscretizationPoints - i; j++) {
                      // Map local triangle parametric coordinates to local tetrahedral parametric coordinates
                      std::vector<int> parUVW = referenceTetrahedronFaceParameterization(iFace, i, j, nIntervals);

                      // Check if this point was already calculated
                      if (parametricToIndex.find(parUVW) == parametricToIndex.end())
                      {
                          // Calculate the real coordinates of the point
                          FieldVector<double, 3> localVector;
                          localVector[0] = double(parUVW[0]) / nIntervals;
                          localVector[1] = double(parUVW[1]) / nIntervals;
                          localVector[2] = double(parUVW[2]) / nIntervals;
                          GlobalVector tmpPoint = thisTetrahedronInt.realCoordinate(localVector);

                          // Shrink each coordinate of each point towards the CoM
                          tmpPoint[0] += shrink_magnitude * (CoM[0] - tmpPoint[0]);
                          tmpPoint[1] += shrink_magnitude * (CoM[1] - tmpPoint[1]);
                          tmpPoint[2] += shrink_magnitude * (CoM[2] - tmpPoint[2]);

                          // Add coordinates to the coordinates array
                          vtk_points.push_back(tmpPoint);

                          // Add point to the point map
                          parametricToIndex[parUVW] = vtk_points.size() - 1;
                      }


                  }
              }


              // Step 2: Split this face into tiny triangles and add them to the triangle array
              for (int i = 0; i < nDiscretizationPoints - 1; i++) {
                  for (int j = 0; j < nDiscretizationPoints - 1 - i; j++) {
                      // Construct two triangles (123) and (234) from 4 points, where 1 is the index point
                      //   3--4
                      //   |\ |
                      //   | \|
                      //   1--2
                      // First triangle always exists, because we iterate over (i,j) in such a way that there is 1 free point at the edge
                      // Second triangle we construct only if point 4 is still in the triangle

                      // Construct triangle (123)
                      std::vector<int> parUVW_1 = referenceTetrahedronFaceParameterization(iFace, i, j, nIntervals);
                      std::vector<int> parUVW_2 = referenceTetrahedronFaceParameterization(iFace, i, j+1, nIntervals);
                      std::vector<int> parUVW_3 = referenceTetrahedronFaceParameterization(iFace, i+1, j, nIntervals);

                      std::vector<int> triangle123;
                      triangle123.push_back(parametricToIndex[parUVW_1]);
                      triangle123.push_back(parametricToIndex[parUVW_2]);
                      triangle123.push_back(parametricToIndex[parUVW_3]);

                      // Add triangle (123)
                      vtk_triangles.push_back(triangle123);

                      // Check if point 4 is still in the triangle
                      if (i + j + 2 < nDiscretizationPoints)
                      {
                          std::vector<int> parUVW_4 = referenceTetrahedronFaceParameterization(iFace, i+1, j+1, nIntervals);

                            std::vector<int> triangle234;
                            triangle234.push_back(parametricToIndex[parUVW_2]);
                            triangle234.push_back(parametricToIndex[parUVW_3]);
                            triangle234.push_back(parametricToIndex[parUVW_4]);

                            // Add triangle (234)
                            vtk_triangles.push_back(triangle234);
                      }
                  }
              }
          }
    }

    // Test Function
    // Takes Tetrahedron - splits it into points uniformly, makes some edges between those points
    // adds all new points to the vtk_points and all new edges to the vtk_edges
    void addTestEdgesFromTetrahedron(
            std::vector<int> & elementDofs,
            Dune::GeometryType & thisElmType,
            int thisElmName,
            int thisElmOrder,
            int thisElmDofNo,
            std::vector< GlobalVector > & nodes,
            std::vector<GlobalVector> & vtk_points,
            std::vector<std::vector<int>> & vtk_edges

            )
    {

        // Get nodes associated with this element
        std::vector<GlobalVector> elementNodes;
        for (int iDof = 0; iDof < thisElmDofNo; iDof++) { elementNodes.push_back(nodes[elementDofs[iDof]]); }

        // Construct Tetrahedral Interpolator used to map parametric coordinates to real interpolated ones
        CurvilinearElementInterpolator <double, 3, 3> thisTetrahedronInt(thisElmType, elementNodes, thisElmOrder);


          // 3) Construct sampling coordinates
          std::map<std::vector<int>, int> parametricToIndex;
          int nDiscretizationPoints = 16;
          int nIntervals = nDiscretizationPoints - 1;

          // Generate all points
          for (int i = 0; i < nDiscretizationPoints; i++) {
              for (int j = 0; j < nDiscretizationPoints - i; j++) {
                  for (int k = 0; k < nDiscretizationPoints - i - j; k++) {
                      // Map local triangle parametric coordinates to local tetrahedral parametric coordinates
                      std::vector<int> parUVW {k, j, i};

                      // Check if this point was already calculated
                      if (parametricToIndex.find(parUVW) == parametricToIndex.end())
                      {
                          // Calculate the real coordinates of the point
                          FieldVector<double, 3> localVector;
                          localVector[0] = double(parUVW[0]) / nIntervals;
                          localVector[1] = double(parUVW[1]) / nIntervals;
                          localVector[2] = double(parUVW[2]) / nIntervals;
                          GlobalVector tmpPoint = thisTetrahedronInt.realCoordinate(localVector);

                          // Add coordinates to the coordinates array
                          vtk_points.push_back(tmpPoint);

                          // Add point to the point map
                          parametricToIndex[parUVW] = vtk_points.size() - 1;
                      }
                  }
              }
          }

          // Construct all edges and add them to the edge array
          for (int i = 0; i < nDiscretizationPoints - 1; i++) {
              for (int j = 0; j < nDiscretizationPoints - 1 - i; j++) {
                  for (int k = 0; k < nDiscretizationPoints - 1 - i - j; k++) {

                      // Construct triangle (123)
                      std::vector<int> parUVW_0 {k,    j,      i};
                      std::vector<int> parUVW_1 {k+1,  j,      i};
                      std::vector<int> parUVW_2 {k,    j+1,    i};
                      std::vector<int> parUVW_3 {k,    j,      i+1};

                      std::vector<int> edge01;    edge01.push_back(parametricToIndex[parUVW_0]);    edge01.push_back(parametricToIndex[parUVW_1]);
                      std::vector<int> edge12;    edge12.push_back(parametricToIndex[parUVW_1]);    edge12.push_back(parametricToIndex[parUVW_2]);
                      std::vector<int> edge20;    edge20.push_back(parametricToIndex[parUVW_2]);    edge20.push_back(parametricToIndex[parUVW_0]);
                      std::vector<int> edge30;    edge30.push_back(parametricToIndex[parUVW_3]);    edge30.push_back(parametricToIndex[parUVW_0]);
                      std::vector<int> edge31;    edge31.push_back(parametricToIndex[parUVW_3]);    edge31.push_back(parametricToIndex[parUVW_1]);
                      std::vector<int> edge32;    edge32.push_back(parametricToIndex[parUVW_3]);    edge32.push_back(parametricToIndex[parUVW_2]);

                      // Add all edges
                      vtk_edges.push_back(edge01);
                      vtk_edges.push_back(edge12);
                      vtk_edges.push_back(edge20);
                      vtk_edges.push_back(edge30);
                      vtk_edges.push_back(edge31);
                      vtk_edges.push_back(edge32);
                  }
              }
          }
    }

    // Given face of the tetrahedron and parametric local coordinates of face gives tetrahedral 3D parametric coordinates
    std::vector<int> referenceTetrahedronFaceParameterization(int iFace, int i, int j, int maxN)
    {
        std::vector<int> rez;
        switch (iFace)
        {
            case 0 : rez.push_back(j);            rez.push_back(i);    rez.push_back(0);    break;
            case 1 : rez.push_back(j);            rez.push_back(0);    rez.push_back(i);    break;
            case 2 : rez.push_back(0);            rez.push_back(j);    rez.push_back(i);    break;
            case 3 : rez.push_back(maxN-i-j);     rez.push_back(j);    rez.push_back(i);    break;
        }
        return rez;
    }

    // Test Function
    // Takes a vector of nodes, a vector of edges and a vector of triangles, writes them all to a .VTK file
    void VTK_Write(
            std::vector<GlobalVector> vtk_points,
            std::vector<std::vector<int>> & vtk_edges,
            std::vector<std::vector<int>> & vtk_triangles
            )
    {
        FILE* vtk_file = fopen("./triangle_test_output.vtk", "w");

        fprintf(vtk_file, "# vtk DataFile Version 2.0\n");
        fprintf(vtk_file, "triangle_test_output, Created by CurvilinearGmshReaderTest\n");
        fprintf(vtk_file, "ASCII\n");
        fprintf(vtk_file, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(vtk_file, "POINTS %d double\n", vtk_points.size() );

        // Write all points
        for (int i = 0; i < vtk_points.size(); i++ ) {
            fprintf(vtk_file, "%lg %lg %lg\n", vtk_points[i][0], vtk_points[i][1], vtk_points[i][2] );
        }

        fprintf(vtk_file, "\n");
        fprintf(vtk_file, "CELLS %d %d\n", vtk_edges.size() + vtk_triangles.size(), vtk_edges.size()*3 + vtk_triangles.size()*4 );

        // Write all elements
        for (int i = 0; i < vtk_edges.size(); i++ )        {    fprintf(vtk_file, "2 %d %d \n", vtk_edges[i][0], vtk_edges[i][1]); }
        for (int i = 0; i < vtk_triangles.size(); i++ )    {    fprintf(vtk_file, "3 %d %d %d\n", vtk_triangles[i][0], vtk_triangles[i][1], vtk_triangles[i][2] ); }

        fprintf(vtk_file, "\n");
        fprintf(vtk_file, "CELL_TYPES %d\n", vtk_edges.size() + vtk_triangles.size() );

        // Write edge codes and triangle codes
        for (int i = 0; i < vtk_edges.size(); i++ )        { fprintf(vtk_file, "3\n"); }
        for (int i = 0; i < vtk_triangles.size(); i++ )    { fprintf(vtk_file, "5\n"); }

        fprintf(vtk_file, "\n");
        fclose(vtk_file);
    }

    // After the element has been read it needs to be added to the factory
    void addElementToFactory(
            std::vector< int > & elementDofs,
            std::vector< unsigned int > &corners,
            std::vector< GlobalVector > &nodes,
            int thisElmName,
            int thisElmDim,
            int thisElmOrder,
            int thisElmPhysicalEntity)
    {
        if (verbose) std::cout << " --- adding element to factory" << std::endl;

        // If it is an element, insert it as such
        if (thisElmDim == dim) {

            if (verbose) std::cout << " --- element dimension is equal to space dimension" << std::endl;

            // Note that hexahedron GeometryType is missing
            switch (thisElmName)
            {
                   case GMSH_EDGE           : factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), corners);   break;
                   case GMSH_TRIANGLE       : factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), corners);   break;
                   case GMSH_QUADRANGLE     : factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim),    corners);   break;
                   case GMSH_TETRAHEDRON    : factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), corners);   break;
                   case GMSH_HEXAHEDRON     : factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim), corners);   break;
                   case GMSH_PRISM          : factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism,dim),   corners);   break;
                   case GMSH_PYRAMID        : factory.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid,dim), corners);   break;
            }
        } else
        {
            if (verbose) std::cout << " --- element is a boundary element" << std::endl;
            // it must be a boundary segment then
            if (insert_boundary_segments)
            {
            	switch (thisElmName)
            	{
            	case GMSH_EDGE :
            		switch (thisElmOrder)
            		{
            		// 2-node line
            		case 1:  factory.insertBoundarySegment(corners);  break;

            		// 3-node line
            		case 2: {
                         array<FieldVector<double,dimWorld>, 3> v;
                         for (int i=0; i<dimWorld; i++) {
                             v[0][i] = nodes[elementDofs[0]][i];
                             v[1][i] = nodes[elementDofs[2]][i];                    // yes, the renumbering is intended!
                             v[2][i] = nodes[elementDofs[1]][i];
                         }

                         BoundarySegment<dim,dimWorld>* newBoundarySegment = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 2, dimWorld >(v[0], v[1], v[2]);
                         factory.insertBoundarySegment(corners, shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment));
                       } break;
            		} break;

            	case GMSH_TRIANGLE :
            		switch (thisElmOrder)
            		{
            		// 3-node triangle
            		case 1: factory.insertBoundarySegment(corners);  break;

            		// 6-node triangle
            		case 2: /* boundarysegment_insert(nodes, elementDofs, vertices); */  break;
            		}

            	}
            }
        }

        // count elements and boundary elements
        if (thisElmDim == dim) {
            if (verbose) std::cout << " --- increasing element count" << std::endl;
            element_index_to_physical_entity.push_back(thisElmPhysicalEntity);
            element_count++;
        } else {
            if (verbose) std::cout << " --- incresing boundary element count" << std::endl;
            boundary_id_to_physical_entity.push_back(thisElmPhysicalEntity);
            boundary_element_count++;
        }
        if (verbose) std::cout << " --- added" << std::endl;
    }





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
    /** \todo doc me */
    static Grid* read (const std::string& fileName, bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<Grid> factory;

      // create parse object
      CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      return factory.createGrid();
    }

    /** \todo doc me */
    static Grid* read (const std::string& fileName,
                       std::vector<int>& boundary_id_to_physical_entity,
                       std::vector<int>& element_index_to_physical_entity,
                       bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<Grid> factory;

      // create parse object
      CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      boundary_id_to_physical_entity.swap(parser.boundaryIdMap());
      element_index_to_physical_entity.swap(parser.elementIndexMap());

      return factory.createGrid();
    }

    /** \todo doc me */
    static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
                      bool verbose = true, bool insert_boundary_segments=true)
    {
      // create parse object
     CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);
    }

    /** \todo doc me */
    static void read (Dune::GridFactory<Grid>& factory,
                      const std::string& fileName,
                      std::vector<int>& boundary_id_to_physical_entity,
                      std::vector<int>& element_index_to_physical_entity,
                      bool verbose = true, bool insert_boundary_segments=true)
 {
  // create parse object
  CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
  parser.read(fileName);

  boundary_id_to_physical_entity.swap(parser.boundaryIdMap());
  element_index_to_physical_entity.swap(parser.elementIndexMap());
 }

  };

  /** \} */

} // namespace Dune

#endif /** DUNE_CURVILINEARGMSHREADER_HH **/
