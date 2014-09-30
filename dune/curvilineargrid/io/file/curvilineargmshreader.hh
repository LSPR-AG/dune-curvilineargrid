// curvilinear gmsh reader functionality - in this implementation of the curvilinear gmsh reader we ONLY read higher order triangles and tetrahedra


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

#include <dune/geometry/type.hh>
#include <dune/geometry/polynomialinterpolation/curvilinearelementinterpolator.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/curvilineargrid/io/file/curvilinearvtkwriter.hh>

namespace Dune
{

  // Stores all info associated with an element, except explicit vertex coordinates
  struct ElementData
  {
      int element_id;
      int gmsh_index;
      int physical_entity_tag;
      int element_tag;
      std::vector<int> process_tags;
      std::vector<int> elementDofs;
  };



  //! dimension independent parts for CurvilinearGmshReaderParser
  template<typename GridType>
  class CurvilinearGmshReaderParser
  {
  protected:
    // If to show debug output
    bool verbose = false;

    // If to save mesh to VTK format after reading
    bool vtk_output = true;

    // Grid Factory
    Dune::GridFactory<GridType>& factory;

    // static data
    static const int dim = GridType::dimension;
    static const int dimWorld = GridType::dimensionworld;
    static_assert( (dimWorld <= 3), "GmshReader requires dimWorld <= 3." );

    // Reading file
    std::string fileName;
    char buf[512];

    // Boundary element indexing
    bool insert_boundary_segments;
    std::vector<int> boundary_element_index_to_physical_entity;
    std::vector<int> internal_element_index_to_physical_entity;

    // Parallel Implementation
    int rank_;
    int size_;

    // A map from GMSH -> Dune for indexing interpolatory points
    std::vector< std::vector< int > > triangularPointRenumberings_;
    std::vector< std::vector< int > > tetrahedralPointRenumberings_;


    std::map<int, unsigned int> nodeofelement; /** \brief map stores the number of vertex indices for a specific element type **/



    // typedefs
    typedef Dune::FieldVector< double, dimWorld > GlobalVector;



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
    GeometryType GMSH_GeometryType(int gmsh_index)
    {
        GeometryType rez;
        int gi = gmsh_index;

             if (gi == 15)                                                                                                       { rez.makeVertex(); }
        else if ((gi == 1) || (gi == 8) || (gi == 26) || (gi == 27) || (gi == 28))                                               { rez.makeLine(); }
        else if ((gi == 2) || (gi == 9) || (gi == 20) || (gi == 21) || (gi == 22) || (gi == 23) || (gi == 24) || (gi == 25))   { rez.makeTriangle(); }
        else if ((gi == 3) || (gi == 10) || (gi == 16))                                                                            { rez.makeQuadrilateral(); }
        else if ((gi == 4) || (gi == 11) || (gi == 29) || (gi == 30) || (gi == 31))                                               { rez.makeTetrahedron(); }
        else if ((gi == 5) || (gi == 12) || (gi == 17) || (gi == 92) || (gi == 93))                                               { rez.makeHexahedron(); }
        else if ((gi == 6) || (gi == 13) || (gi == 18))                                                                            { rez.makePrism(); }
        else if ((gi == 7) || (gi == 14) || (gi == 19))                                                                            { rez.makePyramid(); }
        else  { DUNE_THROW(Dune::IOError, "Unexpected geometry type");  }

/*             // Note that hexahedron GeometryType is missing
             switch (thisElmName)
             {
             case GMSH_EDGE           : factory.insertElement(GeometryType(GeometryType::simplex,dim), corners);   break;
             case GMSH_TRIANGLE       : factory.insertElement(GeometryType(GeometryType::simplex,dim), corners);   break;
             case GMSH_QUADRANGLE     : factory.insertElement(GeometryType(GeometryType::cube,dim),    corners);   break;
             case GMSH_TETRAHEDRON    : factory.insertElement(GeometryType(GeometryType::simplex,dim), corners);   break;
             case GMSH_HEXAHEDRON     : factory.insertElement(GeometryType(GeometryType::simplex,dim), corners);   break;
             case GMSH_PRISM          : factory.insertElement(GeometryType(GeometryType::prism,dim),   corners);   break;
             case GMSH_PYRAMID        : factory.insertElement(GeometryType(GeometryType::pyramid,dim), corners);   break;
             }*/

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

    // Finds the corners in the DUNE-convention DoF vector
    int DUNE_CornerIndex(int cornerNo, GeometryType gt, int thisElmOrder, int thisElmDoF) {

        if (!gt.isSimplex())  { DUNE_THROW(Dune::IOError, "CornerIndex only implemented for Simplex geometries at the moment"); }

        int ind = 0;

        switch(gt.dim()) {
        case 1 :  // EDGE
        {
            switch (cornerNo)
            {
            case 0 : ind = 0;                   break;
            case 1 : ind = thisElmOrder;              break;
            }
        } break;
        case 2 : // TRIANGLE
        {
            switch (cornerNo)
            {
            case 0 : ind = 0;                  break;
            case 1 : ind = thisElmOrder;             break;
            case 2 : ind = thisElmDoF - 1;  break;
            }
        } break;

        case 3 : // TETRAHEDRON
        {
            switch (cornerNo)
            {
            case 0 : ind = 0;                         break;
            case 1 : ind = thisElmOrder;                    break;
            case 2 : ind = thisElmOrder*(thisElmOrder + 3) / 2;   break;
            case 3 : ind = thisElmDoF - 1;         break;
            }
        } break;
        }

        return ind;
    }

    // Says how many d-1 subentities does an element with this name have
    int GEOMETRY_ElementSubentityNo(int elemName)
    {
    	int rez;
    	switch (elemName)
    	{
    		case GMSH_TRIANGLE :     rez = 3;  break;
    		case GMSH_TETRAHEDRON :  rez = 4;  break;
    		default:  DUNE_THROW(Dune::IOError, "GMSH Reader: Not implemented element subentityNo for this element type " );
    	}
    	return rez;
    }

    // Returns d-1 subentity corner local indices, sorted
    std::vector<int> GEOMETRY_ElementSubentityCorners(int elemName, int ind)
    {
    	std::vector<int> rez;
    	switch (elemName)
    	{
    		case GMSH_TRIANGLE :
    		{
    			switch (ind)
    			{
    			case 0 :  rez = std::vector<int> {0, 1};  break;
    			case 1 :  rez = std::vector<int> {1, 2};  break;
    			case 2 :  rez = std::vector<int> {0, 2};  break;
    			default : DUNE_THROW(Dune::IOError, "GMSH Reader: Wrong input arguments for SubentityCorners " );
    			}
    		}  break;
    		case GMSH_TETRAHEDRON :
    		{
    			switch (ind)
    			{
    			case 0 :  rez = std::vector<int> {0, 1, 2};  break;
				case 1 :  rez = std::vector<int> {0, 1, 3};  break;
				case 2 :  rez = std::vector<int> {0, 2, 3};  break;
				case 3 :  rez = std::vector<int> {1, 2, 3};  break;
				default : DUNE_THROW(Dune::IOError, "GMSH Reader: Wrong input arguments for SubentityCorners " );
    			}
    		}  break;
    		default:  DUNE_THROW(Dune::IOError, "GMSH Reader: Not implemented element subentityNo for this element type " );
    	}
    	return rez;
    }


    // Testing if the current element type can be handled by DUNE
    bool checkElementAllowed(int gmsh_element_index)
    {
        bool isAllowedElement = true;
        int thisElmName = GMSH_ElementName(gmsh_element_index);

        // Check if element is polynomial-complete (ask GMSH what that even means I dont know :) )
        isAllowedElement &= !GMSH_ElementIsIncomplete(gmsh_element_index);

        // in 2D only allow lines and triangles
        isAllowedElement &= (dim != 2) || ((thisElmName == GMSH_EDGE)||(thisElmName == GMSH_TRIANGLE));

        // in 3D only allow triangles and tetrahedrons
        isAllowedElement &= (dim != 3) || ((thisElmName == GMSH_TRIANGLE)||(thisElmName == GMSH_TETRAHEDRON));


        // test whether we support the element type at the moment
        if (!isAllowedElement) { DUNE_THROW(Dune::IOError, "GMSH Reader: Have read an element of unexpected type " << gmsh_element_index); }


        return isAllowedElement;
    }

    // Checks whether this element (or boundary element) belongs on this parallel process
    bool elementOnProcess(int eIndex, int eTotal) {
        int eFirst = (eTotal * rank_) / size_;
        int eLast = (eTotal * (rank_+1)) / size_;

        print_debug(" == checkprocess if " + std::to_string(eIndex) + " in [" + std::to_string(eFirst) + "," + std::to_string(eLast) + "]");

        return ((eIndex >= eFirst)&&(eIndex < eLast));
    }


    /** \brief Reads vertices into a map, given a set of indices of vertices that should be read
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  vertices_total                 the total number of vertices specified in the file
     *  \param[in]  vertices                       the map from vertex globalID to vertex coordinate. Here the vertices will be stored
     *  \param[in]  vertex_index_set               the set of globalID's of all vertices that belong to this process
     *  \param[in]  global_to_local_vertex_id_map  the map from vertex globalID to localID
     *
     *  \return The total number of vertices on this process
     *  \note Assumes it is in the correct position in the file
     *
     *  TODO: Make factory.insertVertex() work with globalId
     *
     */
    int readVertices(
            FILE* file,
            int vertices_total,
            std::map<int, GlobalVector> & vertices,
            std::set<int> & vertex_index_set,
            std::map<int, int> & global_to_local_vertex_id_map
            )
    {
        int id;
        GlobalVector x;

        // Iterator starts from 1 because GMSH numbers vertices [1,n]
        for( int i = 1; i <= vertices_total; ++i )
        {
            // If this vertex does not belong to this process, just skip it
            if (vertex_index_set.count(i) == 0)  { fgets(buf, 512, file ); }
            else
            {
            	fscanf(file, "%d ", &id);
            	std::string tmp_out;

                if( id != i )  { DUNE_THROW( Dune::IOError, "Expected id " << i << "(got id " << id << "." ); }
                for (int d = 0; d < dimWorld; d++) {
                	double tmp_coord;
                	fscanf(file, "%lg", &tmp_coord);
                	x[d] = tmp_coord;

                	tmp_out += std::to_string(tmp_coord) + " ";
                }

                // The local id of this vertex is equalt to the number of vertices added so far
                global_to_local_vertex_id_map[i] = vertices.size();

                print_debug("  * Have read vertex " + tmp_out);

                // Maps global id to global coordinate
                vertices[i] = x;

                // Insert vertex into a factory, noting its globalId.
                // Its localId in the factory should be given by the order the vertices are added to the factory
                //factory.insertVertex(x, i);

                fscanf(file, "\n");
            }
        }

        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$EndNodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndNodes"); }

        return vertices.size();
    }

    // Reads the data about this element, everything except the interpolation vertex id's
    // Assumes it is in the correct position in the file
    ElementData readElementSpec(FILE* file)
    {
    	int number_of_tags;
    	ElementData thisElement;

        fscanf(file, "%d %d %d ", &thisElement.element_id, &thisElement.gmsh_index, &number_of_tags);
        print_debug("    * element " + std::to_string(thisElement.element_id) + " has " + std::to_string(number_of_tags) + " tags");

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
        fscanf(file, "%d ", &thisElement.physical_entity_tag);
        fscanf(file, "%d ", &thisElement.element_tag);

        // Possible functionality for mesh partitioning using GMSH
        // TODO: Currently this functionality not available
        for (int k = 2; k < number_of_tags; k++)
        {
            int tmp_tag;
            fscanf(file, "%d ", &tmp_tag);
            thisElement.process_tags.push_back(tmp_tag);
        }

        return thisElement;

    }


    /** \brief Reads all internal element data into a vector. Only reads elements which should be on this process
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  elements_total                 the total number of elements specified in the file
     *  \param[in]  internal_elements_total        the number of internal elements specified in the file
     *  \param[in]  internal_elements              A vector in which the globalID's of internal elements will be stored
     *  \param[in]  vertex_index_set               the set of globalID's of all vertices that belong to this process
     *  \param[in]  element_boundary_linker        A map from a sorted array of globalID's of vertices that make up a boundary to an array of localID's of internal elements to whom this boundary belongs.
     *
     *  \note Assumes it is in the correct position in the file
     *
     *  TODO: Make factory.insertVertex() work with globalId
     */
    void readInternalElements(
    		FILE* file,
    		int elements_total,
    		int internal_elements_total,
    		std::vector<ElementData> & internal_elements,
    		std::set<int> & vertex_index_set,
    		std::map< std::vector<int>, std::vector<int> > & element_boundary_linker
    		)
    {
    	int iSelectElem = 0;

        // Reading element info - tag information and vertex global indices
        for (int i = 0; i < elements_total; i++)
        {
        	// Read the first part of the element info
        	ElementData thisElement = readElementSpec(file);

            // Find if this is a boundary element
            bool on_boundary = (GMSH_ElementDim(GMSH_ElementName(thisElement.gmsh_index)) < dimWorld);

            // Check if we want to read this element at this stage
            if (on_boundary || !elementOnProcess(iSelectElem++, internal_elements_total)) { fgets(buf, 512, file ); }
            else
            {
                // Testing if the current element type can be handled by DUNE
                // *****************************************************************
                checkElementAllowed(thisElement.gmsh_index);
                print_debug("    * element " + std::to_string(thisElement.element_id) + " can be treated by Dune grid ");

                // Obtain all necessary info not to use gmsh_element_index in the following steps
                // *****************************************************
                int thisElmName                = GMSH_ElementName(thisElement.gmsh_index);
                int thisElmOrder               = GMSH_ElementOrder(thisElement.gmsh_index);
                int thisElmDofNo               = GMSH_ElementDofNumber(thisElement.gmsh_index);
                int thisElmCorners             = GMSH_ElementCorners(thisElmName);
                int thisElmSubentities		   = GEOMETRY_ElementSubentityNo(thisElmName);

                // Reading DoF's
                // *************************************************
                for (int iDof = 0; iDof < thisElmDofNo; iDof++) {
                    int tmp_vertex_global_id;
                    fscanf(file, "%d", &tmp_vertex_global_id);
                    print_debug("  --- have read DoF " + std::to_string(tmp_vertex_global_id) );
                    thisElement.elementDofs.push_back(tmp_vertex_global_id);

                    // Insert all used global vertex indexes into set
                    // Note: set ignores request to add an already existing element
                    vertex_index_set.insert(tmp_vertex_global_id);
                }
                fscanf(file, "\n");


                // Add all subentities of this element to a map, such that it is easy afterwards to find
                // boundary elements associated with this process
                // *************************************************************************************

                // 1) get all corners of this process
                std::vector<int> corners;
                for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) { corners.push_back(thisElement.elementDofs[iCorner]); }

                // 2) sort by increasing globalID
                std::sort(corners.begin(), corners.end());

                // 3) Calculate the localID of this element
                int localID = internal_elements.size();

                for (int iSub = 0; iSub < thisElmSubentities; iSub++)
                {
                	// 4) get all subsets associated with this element type
                	std::vector<int> this_subentity_ind = GEOMETRY_ElementSubentityCorners(thisElmName, iSub);
                	std::vector<int> key;

                	print_debug(" for iSub = " + std::to_string(iSub) + " out of " + std::to_string(thisElmSubentities) + " have sub_size = " + std::to_string(this_subentity_ind.size())  );

                	for (int iCoord = 0; iCoord < this_subentity_ind.size(); iCoord++) { key.push_back(corners[this_subentity_ind[iCoord]]); }

                	// 5) Check if map empty for this entry, then add to the map
                	if (element_boundary_linker.find(key) == element_boundary_linker.end()) {
                		std::cout << " -- element " << localID << " has added boundary ";
                		for (int k = 0; k < key.size(); k++) { std::cout << key[k] << " "; }
                		std::cout << std::endl;


                		element_boundary_linker[key] = std::vector<int> (1, localID);
                	} else
                	{
                    	// 6) In this case this is the 2nd element sharing this boundary
                		std::vector<int>  tmp = element_boundary_linker[key];
                		tmp.push_back(localID);
                		element_boundary_linker[key] = tmp;   // Should overwrite previous value

                		std::cout << " -- element " << localID << " shares boundary ";
                		for (int k = 0; k < key.size(); k++) { std::cout << key[k] << " "; }
                		std::cout << " with element " << tmp[0] << std::endl;
                	}
                }

                // correct differences between gmsh and Dune in the local vertex numbering
                // *************************************************
                Gmsh2DuneElementDofNumbering(thisElement.elementDofs, thisElmName, thisElmOrder);

                internal_elements.push_back(thisElement);
            }
        }

        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    /** \brief Reads all internal element data into a vector. Only reads elements which should be on this process
     *
     *  \param[in]  file                           file pointer to read from
     *  \param[in]  elements_total                 the total number of elements specified in the file
     *  \param[in]  boundary_elements              A vector in which the globalID's of boundary elements will be stored
     *  \param[in]  element_boundary_linker        A map from a sorted array of globalID's of vertices that make up a boundary to an array of localID's of internal elements to whom this boundary belongs.
     *  \param[in]  internal_element_boundaries    A vector that stores a vector of localID's of all boundaries linked this element, for each element
     *
     *  \note Assumes it is in the correct position in the file
     *
     *  TODO: Make factory.insertVertex() work with globalId
     */
    void readBoundaryElements(
    		FILE* file,
    		int elements_total,
    		std::vector<ElementData> & boundary_elements,
    		std::map< std::vector<int>, std::vector<int> > & element_boundary_linker,
    		std::vector< std::vector<int> > & internal_element_boundaries
    		)
    {
    	int iSelectElem = 0;

        // Reading element info - tag information and vertex global indices
        for (int i = 0; i < elements_total; i++)
        {
        	// Read the first part of the element info
        	ElementData thisElement = readElementSpec(file);

            // Find if this is a boundary element
            bool on_boundary = (GMSH_ElementDim(GMSH_ElementName(thisElement.gmsh_index)) < dimWorld);

            // If this element is not on the boundary just skip the rest of info on this line
            if (!on_boundary) { fgets(buf, 512, file ); }
            else
            {
                // Testing if the current element type can be handled by DUNE
                // *****************************************************************
                checkElementAllowed(thisElement.gmsh_index);
                print_debug("    * element " + std::to_string(i) + " can be treated by Dune grid ");

                // Obtain all necessary info not to use gmsh_element_index in the following steps
                // *****************************************************
                int thisElmName                = GMSH_ElementName(thisElement.gmsh_index);
                int thisElmOrder               = GMSH_ElementOrder(thisElement.gmsh_index);
                int thisElmDofNo               = GMSH_ElementDofNumber(thisElement.gmsh_index);
                int thisElmCorners             = GMSH_ElementCorners(thisElmName);

                // Reading DoF's
                // *************************************************
                for (int iDof = 0; iDof < thisElmDofNo; iDof++) {
                    int tmp_vertex_global_id;
                    fscanf(file, "%d", &tmp_vertex_global_id);
                    print_debug("  --- have read DoF " + std::to_string(tmp_vertex_global_id) );
                    thisElement.elementDofs.push_back(tmp_vertex_global_id);
                }
                fscanf(file, "\n");



                // Check if associated with any internal element
                // ****************************************************************

                // 1) Obtain corners
                std::vector<int> corners;
                for (int iCorner = 0; iCorner < thisElmCorners; iCorner++) { corners.push_back(thisElement.elementDofs[iCorner]); }

                // 2) sort by increasing globalID
                std::sort(corners.begin(), corners.end());

                // 3) Obtain this boundary element's localID
                int localID = boundary_elements.size();

        		std::cout << " -- boundary " << localID << " checking key ";
        		for (int k = 0; k < corners.size(); k++) { std::cout << corners[k] << " "; }
        		std::cout << std::endl;


                // If this boundary already exists in the map, then add this element
                if (element_boundary_linker.find(corners) != element_boundary_linker.end())
                {
                    std::cout << " -found b.e localID = " << element_boundary_linker[corners][0] << std::endl;

                	// correct differences between gmsh and Dune in the local vertex numbering
                    // *************************************************
                    Gmsh2DuneElementDofNumbering(thisElement.elementDofs, thisElmName, thisElmOrder);

                    // Add boundary element
                    boundary_elements.push_back(thisElement);

                    // If this boundary is linked to an element, write its localID to the internal_element_boundaries for corresponding element
                    std::vector<int> thisBoundaryLinkedElements = element_boundary_linker[corners];
                    if (thisBoundaryLinkedElements.size() > 1) { DUNE_THROW(Dune::IOError, "Have 2 boundaries associated with 1 element. Not implemented yet"); }

                    for (int iLinked = 0; iLinked < thisBoundaryLinkedElements.size(); iLinked++)
                    {
                    	std::cout << iLinked << " " << thisBoundaryLinkedElements.size() << " " << thisBoundaryLinkedElements[iLinked] << std::endl;

                    	internal_element_boundaries[thisBoundaryLinkedElements[iLinked]].push_back(localID);
                    }
                }
            }
        }

        // Finish reading file
        // *************************************************************
        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$EndElements")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndElements"); }
    }


    /** \brief Adds all internal elements to factory, also writes them to .vtk file.
     *
     *  \param[in]  vtk_curv_writer                A class that writes debug output to .vtk file(s)
     *  \param[in]  global_to_local_vertex_id_map  the map from vertex globalID to localID
     *  \param[in]  vertices                       the map from vertex globalID to vertex coordinate.
     *  \param[in]  internal_elements              A vector in which the globalID's of internal elements are stored
     *  \param[in]  internal_element_boundaries    A vector that stores a vector of localID's of all boundaries linked this element, for each element
     *
     *  TODO: The vertex and element vectors are stored twice - once inside the read procedure and once in factory
     *  as they are being added. Maybe possible to save space
     */
    void addInternalElements(
            CurvilinearVTKWriter<dimWorld> & vtk_curv_writer,
            std::map<int, int> & global_to_local_vertex_id_map,
            std::map<int, GlobalVector> & vertices,
            std::vector< ElementData > & internal_elements,
            std::vector< std::vector<int> > & internal_element_boundaries
            )
    {
        // Write elements to factory
        for (int i = 0; i < internal_elements.size(); i++)
        {
              // Obtain all necessary info not to use gmsh_element_index in the following steps
            // *****************************************************
            GeometryType elemType = GMSH_GeometryType(internal_elements[i].gmsh_index);
            int elemName                = GMSH_ElementName(internal_elements[i].gmsh_index);
            int elemOrder               = GMSH_ElementOrder(internal_elements[i].gmsh_index);
            int elemDofNo               = GMSH_ElementDofNumber(internal_elements[i].gmsh_index);

            int elemDim                 = GMSH_ElementDim(elemName);
            int elemCornerNo            = GMSH_ElementCorners(elemName);

            // Compute local DoF vector
            std::vector<int> localDofs;
            for (int iDof = 0; iDof < internal_elements[i].elementDofs.size(); iDof++) {
                localDofs.push_back(global_to_local_vertex_id_map[internal_elements[i].elementDofs[iDof]]);
            }

            print_debug("    * internal_element " + std::to_string(i) + " has dimension " + std::to_string(elemDim) + " and vertex number " + std::to_string(elemCornerNo) + " and physical entity number " + std::to_string(internal_elements[i].physical_entity_tag));

            // TESTING SECTION FOR TETRAHEDRA
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (vtk_output)
            {
                switch (elemName)
                {
                case GMSH_TRIANGLE    : vtk_curv_writer.addTestTriangleFromSurface(elemType, internal_elements[i].elementDofs, vertices, elemOrder, elemDofNo);         break;
                //case GMSH_TETRAHEDRON : vtk_curv_writer.addTestTrianglesFromTetrahedron_NoInterp(elements[i].elementDofs, vertices, elemOrder, elemDofNo);              break;
                case GMSH_TETRAHEDRON : vtk_curv_writer.addTestTrianglesFromTetrahedron(elemType, internal_elements[i].elementDofs, vertices, elemOrder, elemDofNo);    break;
                //case GMSH_TETRAHEDRON : vtk_curv_writer.addTestEdgesFromTetrahedron(elemType, elements[i].elementDofs, vertices, elemOrder, elemDofNo);                 break;
                }
              print_debug("    * internal_element " + std::to_string(i) + " has been added to the VTK triangles  ");
            }
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            //Insert internal element
            //****************************************************

            internal_element_index_to_physical_entity.push_back(internal_elements[i].physical_entity_tag);

            if (insert_boundary_segments) {
                // We should pass an array of boundaries that are connected to this element (internal_element_boundaries, if any)
                // Without this it is impossible to correctly loadBalance the mesh

                //factory.insertElement(elemType, internal_elements[i].element_id, localDofs, internal_elements[i].elementDofs, elemOrder, internal_element_boundaries[i]);
            } else {
                // The corners can in be computed at this stage as well, but it is simpler if we just do following in factory
                // 1) Create LagrangeGeometry from all DoF, store in Metagrid
                // 2) Request corners from LG, send them to HostGrid

                //factory.insertElement(elemType, corners);  // Old
                //factory.insertElement(elemType, internal_elements[i].element_id, localDofs, internal_elements[i].elementDofs, elemOrder);
            }

            print_debug("    * internal_element " + std::to_string(i) + " has been added to the Geometry Factory ");
        }
    }


    /** \brief Adds all boundary elements to factory, also writes them to .vtk file.
     *
     *  \param[in]  vtk_curv_writer                A class that writes debug output to .vtk file(s)
     *  \param[in]  global_to_local_vertex_id_map  the map from vertex globalID to localID
     *  \param[in]  vertices                       the map from vertex globalID to vertex coordinate.
     *  \param[in]  boundary_elements              A vector in which the globalID's of boundary elements are stored
     *
     *  TODO: The vertex and element vectors are stored twice - once inside the read procedure and once in factory
     *  as they are being added. Maybe possible to save space
     */
    void addBoundaryElements(
            CurvilinearVTKWriter<dimWorld> & vtk_curv_writer,
            std::map<int, int> & global_to_local_vertex_id_map,
            std::map<int, GlobalVector> & vertices,
            std::vector< ElementData > & boundary_elements
            )
    {
        // Write elements to factory
        for (int i = 0; i < boundary_elements.size(); i++)
        {
            // Obtain all necessary info not to use gmsh_element_index in the following steps
            // *****************************************************
            GeometryType boundaryType = GMSH_GeometryType(boundary_elements[i].gmsh_index);
            int boundaryName                = GMSH_ElementName(boundary_elements[i].gmsh_index);
            int boundaryOrder               = GMSH_ElementOrder(boundary_elements[i].gmsh_index);
            int boundaryDofNo               = GMSH_ElementDofNumber(boundary_elements[i].gmsh_index);

            int boundaryDim                    = GMSH_ElementDim(boundaryName);
            int boundaryCornerNo            = GMSH_ElementCorners(boundaryName);

            print_debug("    * boundary_element " + std::to_string(i) + " has dimension " + std::to_string(boundaryDim) + " and vertex number " + std::to_string(boundaryCornerNo) + " and physical entity number " + std::to_string(boundary_elements[i].physical_entity_tag));

            // Compute local DoF vector
            std::vector<int> localDofs;
            for (int iDof = 0; iDof < boundary_elements[i].elementDofs.size(); iDof++) {
                localDofs.push_back(global_to_local_vertex_id_map[boundary_elements[i].elementDofs[iDof]]);
            }

            // TESTING SECTION FOR TETRAHEDRA
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (vtk_output)
            {
                switch (boundaryName)
                {
                case GMSH_TRIANGLE    : vtk_curv_writer.addTestTriangleFromSurface(boundaryType, boundary_elements[i].elementDofs, vertices, boundaryOrder, boundaryDofNo);         break;
                //case GMSH_TETRAHEDRON : vtk_curv_writer.addTestTrianglesFromTetrahedron_NoInterp(elements[i].elementDofs, vertices, boundaryOrder, boundaryDofNo);        break;
                case GMSH_TETRAHEDRON : vtk_curv_writer.addTestTrianglesFromTetrahedron(boundaryType, boundary_elements[i].elementDofs, vertices, boundaryOrder, boundaryDofNo);    break;
                //case GMSH_TETRAHEDRON : vtk_curv_writer.addTestEdgesFromTetrahedron(boundaryType, elements[i].elementDofs, vertices, boundaryOrder, boundaryDofNo);        break;
                }
              print_debug("    * boundary_element " + std::to_string(i) + " has been added to the VTK triangles  ");
            }
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            //Insert boundary segments and elements
            //****************************************************

            if (insert_boundary_segments)
            {
                // Adding element to factory
                //factory.insertBoundarySegment(corners);
                //factory.insertBoundarySegment(boundaryType, boundary_elements[i].element_id, localDofs, boundary_elements[i].elementDofs, boundaryOrder);

                // Adding physical tag
                boundary_element_index_to_physical_entity.push_back(boundary_elements[i].physical_entity_tag);

                print_debug("    * boundary_element " + std::to_string(i) + " has been added to the Geometry Factory ");
            }
        }
    }


    // Writes debug info to the command line
    // TODO: Use IFDEF to manipulate between no output, all output, or only master process output
    void print_debug(std::string s)
    {
        if (verbose) { std::cout << "Process_" << rank_ << ": " << s << std::endl; }
    }

  public:

    // TODO: Processor number should be passed to the constructor as argument
    CurvilinearGmshReaderParser(
    		Dune::GridFactory<GridType>& _factory,
    		bool v,
    		bool i
#ifdef HAVE_MPI
    		,MPIHelper &mpihelper
#endif
    		) : factory(_factory), verbose(v), insert_boundary_segments(i)
   {
        // Initialize process parameters
#ifdef HAVE_MPI
    rank_=mpihelper.rank();
    size_=mpihelper.size();
#else
    rank_=0;
    size_=1;
#endif

        print_debug("I am process " + std::to_string(rank_) + " with total processes " + std::to_string(size_));

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

    std::vector<int> & boundaryIdMap()   { return boundary_element_index_to_physical_entity; }

    std::vector<int> & elementIndexMap() { return internal_element_index_to_physical_entity; }

    // This reads the GMSH format to parse the node and element structure
    void read (const std::string& f)
    {
        verbose = false;

        int vertices_total = 0;
        int elements_total = 0;
        int internal_elements_total = 0;
        int boundary_elements_total = 0;

        print_debug(":: using file " + fileName);
        print_debug(":: reading" + std::to_string(dim) + "d curvilinear gmsh grid...");

        // open file name, we use C I/O
        // ***********************************************
        fileName = f;
        FILE* file = fopen(fileName.c_str(),"r");
        if (file==0)  { DUNE_THROW(Dune::IOError, "Could not open " << fileName); }

        // Reading MeshFormat Header
        // ***********************************************
        double version_number;        // process header
        int file_type, data_size;

        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$MeshFormat")!=0)   { DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line"); }
        fscanf(file, "%lg %d %d ", &version_number, &file_type, &data_size);
        if( (version_number < 2.0) || (version_number > 2.3) )  { DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files"); }
        print_debug(":: version " + std::to_string(version_number) + " Gmsh file detected");
        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$EndMeshFormat")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndMeshFormat"); }

        // Reading Node data
        // ***********************************************
        print_debug("----- Reading vertex header ------------------------");
        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$Nodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $Nodes"); }
        fscanf(file, "%d\n", &vertices_total);
        print_debug(":: file contains " + std::to_string(vertices_total) + " vertices");


        //==========================================================
        // VERTEX PASS 1: Put file pointer and skip all vertices
        //==========================================================
        print_debug("----- Vertex-Pass 1: skip all vertices, since we need element data first ---");
        long section_vertex_offset = ftell(file);
        for (int i_vertex = 0; i_vertex < vertices_total; i_vertex++ ) { fgets(buf, 512, file );   print_debug(std::string(buf)); }
        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$EndNodes")!=0)  { DUNE_THROW(Dune::IOError, "expected $EndNodes"); }
        //==========================================================
        // VERTEX PASS 1: Finished
        //==========================================================


        // Reading Element Data
        // *************************************************
        print_debug("----- Reading elements-header -----------------------");
        fscanf(file, "%s\n", buf);
        if (strcmp(buf,"$Elements")!=0)  { DUNE_THROW(Dune::IOError, "expected $Elements"); }
        fscanf(file, "%d\n", &elements_total);
        print_debug(":: file contains " + std::to_string(elements_total) + " elements");



        //=========================================================
        // ELEMENT PASS 1: Count the number of boundary segments
        //=========================================================
        print_debug("----- Elements-Pass 1: counting elements on the boundary---");
        long section_element_offset = ftell(file);

        for (int i = 0; i < elements_total; i++)
        {
            int id, gmsh_element_index;
            fscanf(file, "%d %d ", &id, &gmsh_element_index);

            // Ignore the rest of data on this line
            fgets(buf, 512, file );

            print_debug(std::to_string(id) + " " + std::to_string(gmsh_element_index) + " " + std::string(buf));

            // A boundary segment is defined here as any element with dimension less than world dimension
            if (GMSH_ElementDim(GMSH_ElementName(gmsh_element_index)) < dimWorld )     { boundary_elements_total++; }
            else                                                                         { internal_elements_total++; }
        }


        //==========================================================
        // ELEMENT PASS 2: Read all data associated with element.
        // Note all necessary vertex global indices into a map
        // Map all d-1 subentities of all elements to the element localID's.
        //    - Needed to find boundaries corresponding to each element.
        //==========================================================
        print_debug("----- Elements-Pass 2: reading internal elements for this process ---");
        fseek(file, section_element_offset, SEEK_SET);
        std::set<int> this_process_vertex_index_set;
        std::vector< ElementData > internal_elements;
        std::map< std::vector<int>, std::vector<int> > element_boundary_linker;
        readInternalElements(file, elements_total, internal_elements_total, internal_elements, this_process_vertex_index_set, element_boundary_linker);


        //==========================================================
        // ELEMENT PASS 3: Read all data associated with boundary elements.
        // Only read boundary element if it is associated with this process
        //    - Note: Can not read internal and boundary elements at the same time
        //            because need d-1 subentity map from all internal elements first
        //==========================================================
        print_debug("----- Elements-Pass 3: reading boundary elements for this process ---");
        fseek(file, section_element_offset, SEEK_SET);
        std::vector< ElementData > boundary_elements;
        std::vector< std::vector<int> > internal_element_boundaries;
        readBoundaryElements(file, elements_total, boundary_elements, element_boundary_linker, internal_element_boundaries);




        //==========================================================
        // VERTEX PASS 2: Read the vertices
        // But only the ones that correspond to elements on this process
        //==========================================================
        print_debug("----- Vertex-Pass 2: reading all vertices necessary for this process ---");
        fseek(file, section_vertex_offset, SEEK_SET);
        std::map<int, GlobalVector> vertices;                // Only for testing purposes
        std::map<int, int> global_to_local_vertex_id_map;
        int vertices_on_process = readVertices(file, vertices_total, vertices, this_process_vertex_index_set, global_to_local_vertex_id_map);


        //==========================================================
        // Final Step: Insert boundary segments and elements
        //==========================================================
        print_debug("----- Adding internal boundary elements to factory ---");

        // Testing VTK output
        CurvilinearVTKWriter<dimWorld> vtk_curv_writer;

        addInternalElements(vtk_curv_writer, global_to_local_vertex_id_map, vertices, internal_elements, internal_element_boundaries);
        addBoundaryElements(vtk_curv_writer, global_to_local_vertex_id_map, vertices, boundary_elements);




        print_debug(":: total vertices          = " + std::to_string(vertices_total)          + " of which on this process " + std::to_string(vertices_on_process) );
        print_debug(":: total internal elements = " + std::to_string(internal_elements_total) + " of which on this process " + std::to_string(internal_elements.size()) );
        print_debug(":: total boundary elements = " + std::to_string(boundary_elements_total) + " of which on this process " + std::to_string(boundary_elements.size()) );


        // TESTING SECTION - WRITES TEST ELEMENTS TO .VTK FILE
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (vtk_output)
        {
            vtk_curv_writer.VTK_Write("./curvreader_output_process_" + std::to_string(rank_) + ".vtk");
            print_debug( ">>> wrote test output to .vtk " );
        }

        // Close file
        fclose(file);
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
    /** \brief Reads .GMSH grid, constructing its own factory to put it in
     * Returns Grid from own factory
     * */
    static Grid* read (const std::string& fileName,
#ifdef HAVE_MPI
    					MPIHelper &mpihelper,
#endif
    					bool verbose = true,
    					bool insert_boundary_segments=true
                       )
    {
        // make a grid factory
        Dune::GridFactory<Grid> factory;

        // create parse object
#ifdef HAVE_MPI
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments, mpihelper);
#else
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
#endif

        parser.read(fileName);

        return factory.createGrid();
    }

    /** \brief Reads .GMSH grid, constructing its own factory to put it in
     *  Also receives physical_tag vector for both internal and boundary elements
     *  Returns Grid from own factory
     * */
    static Grid* read (const std::string& fileName,
#ifdef HAVE_MPI
    					MPIHelper &mpihelper,
#endif
                       std::vector<int>& boundary_element_index_to_physical_entity,
                       std::vector<int>& internal_element_index_to_physical_entity,
                       bool verbose = true, bool insert_boundary_segments=true
                       )
    {
        // make a grid factory
        Dune::GridFactory<Grid> factory;

        // create parse object
#ifdef HAVE_MPI
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments, mpihelper);
#else
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
#endif

        parser.read(fileName);

        boundary_element_index_to_physical_entity.swap(parser.boundaryIdMap());
        internal_element_index_to_physical_entity.swap(parser.elementIndexMap());

        return factory.createGrid();
    }

    /** \brief Reads .GMSH grid, factory provided as argument */
    static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
#ifdef HAVE_MPI
    					MPIHelper &mpihelper,
#endif
                      bool verbose = true, bool insert_boundary_segments=true
                     )
    {
#ifdef HAVE_MPI
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments, mpihelper);
#else
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
#endif

        parser.read(fileName);
    }

    /** \brief Reads .GMSH grid, factory provided as argument
     *  Also receives physical_tag vector for both internal and boundary elements
     * */
    static void read (Dune::GridFactory<Grid>& factory,
                      const std::string& fileName,
#ifdef HAVE_MPI
    					MPIHelper &mpihelper,
#endif
                      std::vector<int>& boundary_element_index_to_physical_entity,
                      std::vector<int>& internal_element_index_to_physical_entity,
                      bool verbose = true, bool insert_boundary_segments=true
                     )
    {
        // create parse object
#ifdef HAVE_MPI
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments, mpihelper);
#else
        CurvilinearGmshReaderParser<Grid> parser(factory,verbose,insert_boundary_segments);
#endif
          parser.read(fileName);

          boundary_element_index_to_physical_entity.swap(parser.boundaryIdMap());
          internal_element_index_to_physical_entity.swap(parser.elementIndexMap());
    }
  };


} // namespace Dune

#endif /** DUNE_CURVILINEARGMSHREADER_HH **/
