/*******************************************************************
 * Curvilinear VTK writer
 * 
 * author: Aleksejs Fomins
 * date: 01.09.2014 - created
 * 
 * description:
 * Generates VTK, VTU and PVTU files to visualise curvilinear elements.
 * Curvilinear elements added to the writer are automatically discretized using
 * fixed interval linear elements, which are consequently written to a file
 * 
 *******************************************************************/


#ifndef DUNE_CURVILINEARVTKWRITER_HH
#define DUNE_CURVILINEARVTKWRITER_HH

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

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>



namespace Dune
{
  



namespace VTKEntitySubset
{

typedef std::vector<int>                            IndexVector;
typedef std::vector< IndexVector >                  ElemGridEnumerate;

typedef std::map<std::vector<int>, int>             LocalCoordinate2GlobalIdMap;
typedef std::vector< LocalCoordinate2GlobalIdMap >  Coord2GlobalMapVector;

typedef std::vector<std::vector<int> >              SubEntityIndexVector;

// Codimensions of entity types for better code readability
const int   VERTEX_CODIM   = 3;
const int   EDGE_CODIM     = 2;
const int   FACE_CODIM     = 1;
const int   ELEMENT_CODIM  = 0;




/** \brief Connects discretization points of an entity<mydim> using linear subentities<subdim>, adds these subentities for writing to VTK
 *
 *  \param[in]  entityEnumeratorReduced        A vector of keys that enumerate triangular discretization points.
 *                                             Reduced means that enumerator's points-per-edge is one less than the number used to discretize the element.
 *  \param[in]  nInterval                      Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
 *  \param[in]  thisElmPhysTag                 Physical Tag of the element being discretized
 *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
 *
 */
template <int mycodim, int subcodim>
SubEntityIndexVector refineEntitySubset(ElemGridEnumerate & entityEnumeratorReduced, int nInterval, LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	std::cout << "ERROR: called refineErrorSubset generic (non-specialized) form that does nothing!" << std::endl;
	assert(0);
    return SubEntityIndexVector();
}


/** \brief Connects discretization points of a curvilinear edge using linear edges, adds these edges for writing to VTK */
template <>
SubEntityIndexVector refineEntitySubset<EDGE_CODIM, EDGE_CODIM>(
        ElemGridEnumerate & edgeEnumeratorReduced,
        int nInterval,
        LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	SubEntityIndexVector rez;

    // Construct all edges and add them to the edge array
    for (int i = 0; i < edgeEnumeratorReduced.size(); i++)
    {
        int x = edgeEnumeratorReduced[i][0];

        // Construct triangle (123)
        std::vector<int> parUV_0 {x};
        std::vector<int> parUV_1 {x + 1};

        std::vector<int> edge01 {parametricToIndex[parUV_0], parametricToIndex[parUV_1]};

        rez.push_back(edge01);
    }
    return rez;
}


/** \brief Connects discretization points of a curvilinear triangle using linear edges, adds these edges for writing to VTK */
template <>
SubEntityIndexVector refineEntitySubset<FACE_CODIM, EDGE_CODIM> (
        ElemGridEnumerate & triangleEnumeratorReduced,
        int nInterval,
        LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	SubEntityIndexVector rez;

    // Construct all edges and add them to the edge array
    for (int i = 0; i < triangleEnumeratorReduced.size(); i++)
    {
        int x = triangleEnumeratorReduced[i][0];
        int y = triangleEnumeratorReduced[i][1];

        // Construct triangle (123)
        std::vector<int> parUV_0 {x,      y,   };
        std::vector<int> parUV_1 {x + 1,  y,   };
        std::vector<int> parUV_2 {x,      y + 1};

        std::vector<int> edge01 { parametricToIndex[parUV_0], parametricToIndex[parUV_1] };
        std::vector<int> edge12 { parametricToIndex[parUV_1], parametricToIndex[parUV_2] };
        std::vector<int> edge20 { parametricToIndex[parUV_2], parametricToIndex[parUV_0] };

        // Add all edges
        rez.push_back(edge01);
        rez.push_back(edge12);
        rez.push_back(edge20);
    }
    return rez;
}


/** \brief Connects discretization points of a curvilinear tetrahedron using linear edges, adds these edges for writing to VTK */
template<>
SubEntityIndexVector refineEntitySubset<ELEMENT_CODIM, EDGE_CODIM> (
        ElemGridEnumerate & tetrahedronEnumeratorReduced,
        int nInterval,
        LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	SubEntityIndexVector rez;

    // Construct all edges and add them to the edge array
    for (int i = 0; i < tetrahedronEnumeratorReduced.size(); i++)
    {
        int x = tetrahedronEnumeratorReduced[i][0];
        int y = tetrahedronEnumeratorReduced[i][1];
        int z = tetrahedronEnumeratorReduced[i][2];

        // Construct triangle (123)
        std::vector<int> parUVW_0 {x,      y,      z    };
        std::vector<int> parUVW_1 {x + 1,  y,      z    };
        std::vector<int> parUVW_2 {x,      y + 1,  z    };
        std::vector<int> parUVW_3 {x,      y,      z + 1};

        std::vector<int> edge01 { parametricToIndex[parUVW_0], parametricToIndex[parUVW_1] };
        std::vector<int> edge12 { parametricToIndex[parUVW_1], parametricToIndex[parUVW_2] };
        std::vector<int> edge20 { parametricToIndex[parUVW_2], parametricToIndex[parUVW_0] };
        std::vector<int> edge30 { parametricToIndex[parUVW_3], parametricToIndex[parUVW_0] };
        std::vector<int> edge31 { parametricToIndex[parUVW_3], parametricToIndex[parUVW_1] };
        std::vector<int> edge32 { parametricToIndex[parUVW_3], parametricToIndex[parUVW_2] };

        // Add all edges
        rez.push_back(edge01);
        rez.push_back(edge12);
        rez.push_back(edge20);
        rez.push_back(edge30);
        rez.push_back(edge31);
        rez.push_back(edge32);
    }
    return rez;
}


/** \brief Connects discretization points of a curvilinear triangle using linear triangles, adds these edges for writing to VTK */
template<>
SubEntityIndexVector refineEntitySubset<FACE_CODIM, FACE_CODIM> (
        ElemGridEnumerate & triangleEnumeratorReduced,
        int nInterval,
        LocalCoordinate2GlobalIdMap & parametricToIndex
        )
{
	SubEntityIndexVector rez;

    for (int i = 0; i < triangleEnumeratorReduced.size(); i++) {
        // Construct two triangles (123) and (234) from 4 points, where 1 is the index point
        //   3--4
        //   |\ |
        //   | \|
        //   1--2
        // First triangle always exists, because we iterate over (i,j) in such a way that there is 1 free point at the edge
        // Second triangle we construct only if point 4 is still in the triangle

        int x = triangleEnumeratorReduced[i][0];
        int y = triangleEnumeratorReduced[i][1];

        // Construct triangle (123)
        std::vector<int> parUV_1 {x    , y    };
        std::vector<int> parUV_2 {x + 1, y    };
        std::vector<int> parUV_3 {x    , y + 1};

        std::vector<int> triangle123 { parametricToIndex[parUV_1], parametricToIndex[parUV_2], parametricToIndex[parUV_3] };

        // Add triangle (123)
        rez.push_back(triangle123);

        // Check if point 4 is still within the original triangle
        // If yes then add its symmetric triangle too
        if (x + y + 1 < nInterval)
        {
            std::vector<int> parUV_4 {x + 1, y + 1};
            std::vector<int> triangle234 { parametricToIndex[parUV_2], parametricToIndex[parUV_3], parametricToIndex[parUV_4] };

            // Add triangle (234)
            rez.push_back(triangle234);
        }
    }
    return rez;
}


/** \brief Connects discretization points of a curvilinear tetrahedron using linear triangles, adds these edges for writing to VTK */
template<>
SubEntityIndexVector refineEntitySubset<ELEMENT_CODIM, FACE_CODIM> (
        ElemGridEnumerate & tetrahedralEnumeratorReduced,
        int nInterval,
        LocalCoordinate2GlobalIdMap & parametricToIndex
        )
{
	SubEntityIndexVector rez;

    ElemGridEnumerate triangleEnumerator = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<2>(nInterval);
    ElemGridEnumerate triangleEnumeratorReduced = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<2>(nInterval - 1);

    // Takes tetrahedral discretization point key to globalID 3D map, and splits it into 4 triangle 2D maps, corresponding to the faces of the tetrahedron
    Coord2GlobalMapVector consistingTriangles(4);
    for (int i = 0; i < triangleEnumerator.size(); i++)
    {
        int x = triangleEnumerator[i][0];
        int y = triangleEnumerator[i][1];
        int z = nInterval - x - y;

        std::vector<int> faceInd_0 {x, y, 0};
        std::vector<int> faceInd_1 {x, 0, y};
        std::vector<int> faceInd_2 {0, x, y};
        std::vector<int> faceInd_3 {z, x, y};

        consistingTriangles[0][triangleEnumerator[i]] = parametricToIndex[faceInd_0];
        consistingTriangles[1][triangleEnumerator[i]] = parametricToIndex[faceInd_1];
        consistingTriangles[2][triangleEnumerator[i]] = parametricToIndex[faceInd_2];
        consistingTriangles[3][triangleEnumerator[i]] = parametricToIndex[faceInd_3];
    }

    // Discretizes resulting triangles using triangle discretization routine
    for (int iFace = 0; iFace < 4; iFace++)
    {
    	SubEntityIndexVector faceSub = refineEntitySubset<FACE_CODIM, FACE_CODIM> (triangleEnumeratorReduced, nInterval, consistingTriangles[iFace]);
    	rez.insert (rez.begin(), faceSub.begin(), faceSub.end());
    }

    return rez;
}

} // Namespace VTKEntitySubset










  template<class GridType>
  class CurvilinearVTKWriter
  {

/** \brief This class takes curved elements, samples them on a grid, connects points into mesh of
 * small straight-sided triangles and edges, ands writes them to the .vtk file
 *
 */
  public:

      static const int dimension = GridType::dimension;
      static const bool isCached = GridType::is_cached;

      // Codimensions of entity types for better code readability
      static const int   VERTEX_CODIM   = 3;
      static const int   EDGE_CODIM     = 2;
      static const int   FACE_CODIM     = 1;
      static const int   ELEMENT_CODIM  = 0;

      typedef FieldVector< double, dimension >            GlobalVector;
      typedef std::vector<int>                            IndexVector;
      typedef std::vector<int>                            TagVector;
      typedef std::vector<std::vector<int> >              SubEntityIndexVector;
      typedef std::map<std::vector<int>, int>             LocalCoordinate2GlobalIdMap;
      typedef std::vector< LocalCoordinate2GlobalIdMap >  Coord2GlobalMapVector;
      typedef typename LocalCoordinate2GlobalIdMap::iterator  LC2GIMapIter;

      typedef std::map<std::string, unsigned int>         FieldNameMap;
      typedef std::map<int, GlobalVector>                 FieldCoordMap;

      typedef typename FieldNameMap::iterator             FieldNameMapIter;
      typedef typename FieldCoordMap::iterator            FieldCoordMapIter;

      typedef std::vector< IndexVector >                  ElemGridEnumerate;

      // Logging Message Typedefs
      typedef typename GridType::LoggingMessage           LoggingMessage;
      static const unsigned int LOG_CATEGORY_DEBUG = LoggingMessage::Category::DEBUG;

      static const unsigned int BOUNDARY_SEGMENT_PARTITION_TYPE = GridType::BOUNDARY_SEGMENT_PARTITION_TYPE;

      // VTK constants
      const std::string VTK_XML_VERSION = "1.0";
      const std::string VTK_GRID_TYPE = "UnstructuredGrid";
      const std::string VTK_VTU_VERSION = "0.1";
      const std::string VTK_BYTE_ORDER = "LittleEndian";


  public:

    CurvilinearVTKWriter (MPIHelper &mpihelper)
    {
    	tmpNInterval_ = 0;         // Defined later
        rank_ = mpihelper.rank();
        size_ = mpihelper.size();
    }


    CurvilinearVTKWriter (int rank, int size)
    {
    	tmpNInterval_ = 0;         // Defined later
        rank_ = rank;
        size_ = size;
    }


    /** \brief Takes curvilinear element, discretizes it into linear element, adds linear elements to the writer
     *
     *  \param[in]  geomtype                    Geometry Type of the element
     *  \param[in]  nodeSet                 Coordinates of the element's interpolation points
     *  \param[in]  tagSet                  Set of 0 to 3 tags, in priority order. [Physical Tag, Structural Type, ProcessRank]
     *  \param[in]  elementOrder                   Curvilinear interpolation order of the element
     *  \param[in]  nDiscretizationPoint          Number of discretization points-per-edge
     *  \param[in]  interpolate                    If set to false, instead of interpolating, the provided interpolation points will be used as discretization points
     *  \param[in]  interpolate                    If set to true, there will be gaps between all elements, obtained by scaling them away from the center of mass
     *  \param[in]  writeEdgeData                  If true, discretization points are connected using linear edges forming triangular mesh, and added to writer
     *  \param[in]  writeTriangleData              If true, discretization points are connected using linear triangles forming triangular mesh, and added to writer
     *
     *    \note Discretization is performed using a regular grid over the reference element, for example,
     *    below the stars represent the triangular discretization points with nDiscretizationPoint = 5
     *
     *  *
     *  **
     *  ***
     *  ****
     *  *****
     */
    template<int mydim>
    void addCurvilinearElement(
            const Dune::GeometryType & geomtype,
            const std::vector<GlobalVector> & nodeSet,
            std::vector<int> & tagSet,
            int elementOrder,
            int nDiscretizationPoint,
            bool interpolate,
            bool explode,
            bool writeEdgeData,
            bool writeTriangleData

            )
    {
        int   thisElmPhysTag        = (tagSet.size() > 0 ) ? tagSet[0] : 0;
        int   thisElmPartitionType  = (tagSet.size() > 1 ) ? tagSet[1] : 0;
        int   thisElmProcessRank    = (tagSet.size() > 2 ) ? tagSet[2] : 0;

        // 0.0 - no shrinking, 0.99 - very small element (must be < 1)
        double shrinkMagnitude = explode ? 0.2 : 0.0;

        // Expand all domain boundary surfaces a little bit so that they do not interlay with element surfaces
        double boundaryMagnification;
        std::string pname;

        if (thisElmPartitionType == BOUNDARY_SEGMENT_PARTITION_TYPE)
        {
        	boundaryMagnification = 1.2;
        	pname = "BoundarySegment";
        } else
        {
        	boundaryMagnification = 1.0;
        	pname = Dune::PartitionName(static_cast<Dune::PartitionType> (thisElmPartitionType));
        }

        // It is not possible to write triangle data for edges
        if (mydim == 1)  { writeTriangleData = false; }

        std::stringstream log_message;
        log_message << "VTK_WRITER: Adding a curvilinear element Type=" << Dune::CurvilinearGeometryHelper::geometryName(geomtype);
        log_message << " Order="               << elementOrder;
        log_message << " PhysicalTag="         << thisElmPhysTag;
        log_message << " StructuralType="      << pname;
        log_message << " ProcessRank="         << thisElmProcessRank;
        log_message << " nDiscretization="     << nDiscretizationPoint;
        log_message << " useInterpolation="    << interpolate;
        log_message << " explodeElements="     << explode;
        log_message << " explosionMagnitude="  << shrinkMagnitude;
        log_message << " writeVTK_edges="      << writeEdgeData;
        log_message << " writeVTK_triangles="  << writeTriangleData;
        //log_message << " vertices=" << Dune::VectorHelper::vector2string(nodeSet);
        LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, log_message.str());

        addCurvilinearSimplex<mydim>(geomtype, nodeSet, tagSet, elementOrder, nDiscretizationPoint, shrinkMagnitude, boundaryMagnification, interpolate, writeEdgeData, writeTriangleData);
    }

    /** \brief After an element has been added, one has the option attach one or more vector fields to the discretization vertices */
    template <class VTKElementaryFunction>
    void addField(std::string fieldname, VTKElementaryFunction & vtkfunction)
    {
    	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Adding element field " + fieldname);

    	int fieldIndex;                                                 // The index of the field associated with this field name
    	FieldNameMapIter iter = fieldName2Index_.find(fieldname);

    	if (iter != fieldName2Index_.end() )  {       // If field exists, find its index
    		fieldIndex = (*iter).second;
    	} else {                                      // Otherwise, make new index and add field name to them map
    		fieldIndex = vtkFieldVector_.size();
    		fieldName2Index_.insert(std::pair<std::string, int>(fieldname, fieldIndex));
    		vtkFieldVector_.resize(fieldIndex + 1);
    	}

    	// Loop over all points
    	for (LC2GIMapIter iter = tmpParameter2Index_.begin(); iter != tmpParameter2Index_.end(); iter++)
    	{
    		std::vector<int> integerCoordinate = (*iter).first;
    		unsigned int     vertexIndex  = (*iter).second;

    		// Compute local coordinate
    		GlobalVector local;
    		for (int i = 0; i < dimension; i++)  { local[i] = (double((*iter).first[i])) / tmpNInterval_; }

    		// Evaluate field
    		GlobalVector field = vtkfunction.evaluate(local);

    		// Append coordinate index and field
    		vtkFieldVector_[fieldIndex].insert(std::pair<int, GlobalVector>(vertexIndex, field));
    	}
    }


    // Writes serial VTK file
    // Takes a vector of vertices, a vector of edges and a vector of triangles, writes them all to a .VTK file
    void writeVTK( std::string filename)
    {
        FILE* vtkFile = fopen(filename.c_str(), "w");

        int nElements = vtkEdgeVertexIndex_.size() + vtkTriangleVertexIndex_.size();
        int nCells = vtkEdgeVertexIndex_.size()*3 + vtkTriangleVertexIndex_.size()*4;

        fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
        fprintf(vtkFile, "CurvilinearGmshReader test output\n");
        fprintf(vtkFile, "ASCII\n");
        fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(vtkFile, "POINTS %d double\n", vtkPoint_.size() );

        // Write all points
        for (int i = 0; i < vtkPoint_.size(); i++ ) {
            for (int d = 0; d < dimension; d++)  { fprintf(vtkFile, "%lg ", vtkPoint_[i][d]); }
            fprintf(vtkFile, "\n");
        }

        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELLS %d %d\n", nElements, nCells );

        // Write all elements
        for (int i = 0; i < vtkEdgeVertexIndex_.size(); i++ )        { fprintf(vtkFile, "2 %d %d \n", vtkEdgeVertexIndex_[i][0], vtkEdgeVertexIndex_[i][1]); }
        for (int i = 0; i < vtkTriangleVertexIndex_.size(); i++ )    { fprintf(vtkFile, "3 %d %d %d\n", vtkTriangleVertexIndex_[i][0], vtkTriangleVertexIndex_[i][1], vtkTriangleVertexIndex_[i][2] ); }

        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELL_TYPES %d\n", nElements );

        // Write edge and triangle cell types
        for (int i = 0; i < vtkEdgeVertexIndex_.size(); i++ )        { fprintf(vtkFile, "3\n"); }
        for (int i = 0; i < vtkTriangleVertexIndex_.size(); i++ )    { fprintf(vtkFile, "5\n"); }


        // If are defined fields associated with vertices, then write them too
        if (fieldName2Index_.size() > 0)
        {
        	fprintf(vtkFile, "\n");
        	fprintf(vtkFile, "POINT_DATA %d\n", vtkPoint_.size() );

        	for (FieldNameMapIter iterName = fieldName2Index_.begin(); iterName != fieldName2Index_.end(); ++iterName)
        	{
        		std::string fieldName = (*iterName).first;
        		int fieldIndex  = (*iterName).second;

        		fprintf(vtkFile, "VECTORS %s FLOAT\n", fieldName.c_str());

                for (int i = 0; i < vtkPoint_.size(); i++ ) {
                	// If there is no field defined for this vertex, just print a zero vector for consistency
                	FieldCoordMapIter iterField = vtkFieldVector_[fieldIndex].find(i);
                	if (iterField == vtkFieldVector_[fieldIndex].end())  { fprintf(vtkFile, "0.0 0.0 0.0\n"); }
                	else {
                		GlobalVector v = (*iterField).second;
                		fprintf(vtkFile, "%lg %lg %lg\n", v[0], v[1], v[2]);
                	}
                }
                fprintf(vtkFile, "\n");
        	}
        }



        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELL_DATA %d\n", nElements );

        // Write edge and triangle Structural type
        fprintf(vtkFile, "SCALARS physicalTag FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgePhysicalTag_.size(); i++ )        { fprintf(vtkFile, "%d\n", vtkEdgePhysicalTag_[i]); }
        for (int i = 0; i < vtkTrianglePhysicalTag_.size(); i++ )    { fprintf(vtkFile, "%d\n", vtkTrianglePhysicalTag_[i]); }

        // Write edge and triangle physicalTags
        fprintf(vtkFile, "SCALARS structuralType FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgeStructuralType_.size(); i++ )        { fprintf(vtkFile, "%d\n", vtkEdgeStructuralType_[i]); }
        for (int i = 0; i < vtkTriangleStructuralType_.size(); i++ )    { fprintf(vtkFile, "%d\n", vtkTriangleStructuralType_[i]); }

        // Write edge and triangle provider process ranks
        fprintf(vtkFile, "SCALARS processRank FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgeProcessRank_.size(); i++ )        { fprintf(vtkFile, "%d\n", vtkEdgeProcessRank_[i]); }
        for (int i = 0; i < vtkTriangleProcessRank_.size(); i++ )    { fprintf(vtkFile, "%d\n", vtkTriangleProcessRank_[i]); }

        // Empty line at the end of file
        fprintf(vtkFile, "\n");
        fclose(vtkFile);
    }

    // Writes a PVTU parallel file (no data in this file)
    void writePVTU(std::string filenameBody, int size)
    {
    	std::string filename = filenameBody + ".pvtu";
        FILE* pvtuFile = fopen(filename.c_str(), "w");

        // Write header
        // *****************************************************
        fprintf(pvtuFile, "<?xml version=\"%s\"?>\n", VTK_XML_VERSION.c_str());
        fprintf(pvtuFile, "<VTKFile type=\"P%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(pvtuFile, "<P%s GhostLevel=\"0\">\n", VTK_GRID_TYPE.c_str());

        // PointData could be provided here
        // (For example, tags associated with vertices, which we do not have at the moment)
        // *****************************************************
        // <PPointData>...</PPointData>


        // Write edge and triangle physicalTags
        // *****************************************************
        fprintf(pvtuFile, "<PCellData Scalars=\"physicalTag\">\n");
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" Name=\"physicalTag\" NumberOfComponents=\"1\"/>\n");

        // Write edge and triangle structural types
        // *****************************************************
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\"/>\n");

        // Write edge and triangle process ranks
        // *****************************************************
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" Name=\"processRank\" NumberOfComponents=\"1\"/>\n");
        fprintf(pvtuFile, "</PCellData>\n");


        // Write coordinates of vertices
        // *****************************************************
        fprintf(pvtuFile, "<PPoints>\n");
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
        fprintf(pvtuFile, "</PPoints>\n");


        // Write element information
        // *****************************************************
        fprintf(pvtuFile, "<PCells>\n");
        fprintf(pvtuFile, "<PDataArray type=\"Int32\" Name=\"connectivity\"/>\n");   // Vertex indices
        fprintf(pvtuFile, "<PDataArray type=\"Int32\" Name=\"offsets\"/>\n");
        fprintf(pvtuFile, "<PDataArray type=\"UInt8\" Name=\"types\"/>\n");   // Element types
        fprintf(pvtuFile, "</PCells>\n");


        // Write all .vtu data files file
        // *****************************************************
        for (int iProc = 0; iProc < size; iProc++ )
        {
            std::string vtuFilename = filenameBody + "_process_" + std::to_string(iProc) + ".vtu";
            fprintf(pvtuFile, "<Piece  Source=\"%s\"/>\n", vtuFilename.c_str());
        }


        // Finish writing file
        // *****************************************************
        fprintf(pvtuFile, "</P%s>\n", VTK_GRID_TYPE.c_str());
        fprintf(pvtuFile, "</VTKFile>\n");
        fclose(pvtuFile);
    }

    // Writes serial VTU file
    void writeVTU( std::string filename)
    {
        FILE* vtuFile = fopen(filename.c_str(), "w");

        int nVertices = vtkPoint_.size();
        int nElements = vtkEdgeVertexIndex_.size() + vtkTriangleVertexIndex_.size();

        // Compute offsets which is a general way to determine the number of vertices per element
        std::vector<int> offsets;
        int tmp_offset = 0;
        for (int i = 0; i < vtkEdgeVertexIndex_.size(); i++ )        { tmp_offset += 2;  offsets.push_back(tmp_offset); }
        for (int i = 0; i < vtkTriangleVertexIndex_.size(); i++ )    { tmp_offset += 3;  offsets.push_back(tmp_offset); }

        // Write header
        // *****************************************************
        fprintf(vtuFile, "<?xml version=\"%s\"?>\n", VTK_XML_VERSION.c_str());
        fprintf(vtuFile, "<VTKFile type=\"%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(vtuFile, "<%s>\n", VTK_GRID_TYPE.c_str());
        fprintf(vtuFile, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nVertices, nElements);


        // PointData could be provided here
        // (For example, tags associated with vertices, which we do not have at the moment)
        // *****************************************************
        // <PointData>...</PointData>


        // Write edge and triangle physicalTags and structural types
        // *****************************************************
        fprintf(vtuFile, "<CellData Scalars=\"physicalTag\">\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"physicalTag\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle physicalTags
        for (int i = 0; i < vtkEdgePhysicalTag_.size(); i++ )        { fprintf(vtuFile, "%d ", vtkEdgePhysicalTag_[i]); }
        for (int i = 0; i < vtkTrianglePhysicalTag_.size(); i++ )    { fprintf(vtuFile, "%d ", vtkTrianglePhysicalTag_[i]); }

        fprintf(vtuFile, "\n</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle structural type
        for (int i = 0; i < vtkEdgeStructuralType_.size(); i++ )        { fprintf(vtuFile, "%d ", vtkEdgeStructuralType_[i]); }
        for (int i = 0; i < vtkTriangleStructuralType_.size(); i++ )    { fprintf(vtuFile, "%d ", vtkTriangleStructuralType_[i]); }

        fprintf(vtuFile, "\n</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"processRank\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle provider process ranks
        for (int i = 0; i < vtkEdgeProcessRank_.size(); i++ )        { fprintf(vtuFile, "%d ", vtkEdgeProcessRank_[i]); }
        for (int i = 0; i < vtkTriangleProcessRank_.size(); i++ )    { fprintf(vtuFile, "%d ", vtkTriangleProcessRank_[i]); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</CellData>\n");


        // Write coordinates of vertices
        // *****************************************************
        fprintf(vtuFile, "<Points>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

        // Write all points
        for (int i = 0; i < vtkPoint_.size(); i++ ) {
            for (int d = 0; d < dimension; d++)  { fprintf(vtuFile, "%lg ", vtkPoint_[i][d]); }
            fprintf(vtuFile, "\n");
        }

        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</Points>\n");


        // Write element information
        // *****************************************************
        fprintf(vtuFile, "<Cells>\n");
        fprintf(vtuFile, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");   // Vertex indices

        for (int i = 0; i < vtkEdgeVertexIndex_.size(); i++ )        { fprintf(vtuFile, "%d %d ", vtkEdgeVertexIndex_[i][0], vtkEdgeVertexIndex_[i][1]); }
        for (int i = 0; i < vtkTriangleVertexIndex_.size(); i++ )    { fprintf(vtuFile, "%d %d %d ", vtkTriangleVertexIndex_[i][0], vtkTriangleVertexIndex_[i][1], vtkTriangleVertexIndex_[i][2] ); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

        for (int i = 0; i < offsets.size(); i++ )        { fprintf(vtuFile, "%d ", offsets[i]); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");   // Element types

        for (int i = 0; i < vtkEdgeVertexIndex_.size(); i++ )        { fprintf(vtuFile, "3 "); }
        for (int i = 0; i < vtkTriangleVertexIndex_.size(); i++ )    { fprintf(vtuFile, "5 "); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</Cells>\n");


        // Finish writing file
        // *****************************************************
        fprintf(vtuFile, "</Piece>\n");
        fprintf(vtuFile, "</%s>\n", VTK_GRID_TYPE.c_str());
        fprintf(vtuFile, "</VTKFile>\n");
        fclose(vtuFile);
    }

    // Writes a VTU file on all processes and a PVTU on Master Process
    void writeParallelVTU(std::string filenameBody)
    {
        // Write a PVTU file on master process
        if (rank_ == 0) { writePVTU(filenameBody, size_); }

        // Write a VTU file on all processes
        writeVTU(filenameBody  + "_process_" + std::to_string(rank_) + ".vtu");
    }


  protected:

    // ***********************************************
    // Auxiliary Methods
    // ***********************************************

    /** \brief Adds an edge which will be explicitly written to the file
     *  Note: Only add tags if they are provided by the user
     * */
    void addDiscretizationEdge(IndexVector & indexSet, std::vector<int> & tagSet)
    {
        vtkEdgeVertexIndex_.push_back(indexSet);
        vtkEdgePhysicalTag_.push_back(    (tagSet.size() > 0) ? tagSet[0] : 0);
        vtkEdgeStructuralType_.push_back( (tagSet.size() > 1) ? tagSet[1] : 0);
        vtkEdgeProcessRank_.push_back(    (tagSet.size() > 2) ? tagSet[2] : 0);
    }

    /** \brief Adds a triangle which will be explicitly written to the file
     *  Note: Only add tags if they are provided by the user
     * */
    void addDiscretizationTriangle(IndexVector & indexSet, std::vector<int> & tagSet)
    {
        vtkTriangleVertexIndex_.push_back(indexSet);
        vtkTrianglePhysicalTag_.push_back(    (tagSet.size() > 0) ? tagSet[0] : 0);
        vtkTriangleStructuralType_.push_back( (tagSet.size() > 1) ? tagSet[1] : 0);
        vtkTriangleProcessRank_.push_back(    (tagSet.size() > 2) ? tagSet[2] : 0);
    }

    /** \brief Calculates the centre of mass of a vector of points (equal-weighted) */
    GlobalVector vectorCentreOfMass( const std::vector<GlobalVector> & cornerVector)
    {
        GlobalVector rez;
        for (int i = 0; i < cornerVector.size(); i++) { rez += cornerVector[i]; }
        rez /= cornerVector.size();
        return rez;
    }






    /** \brief Checks if the tetrahedral discretization point corresponds to tetrahedral boundary
     *
     *  \param[in]  p                               tetrahedral discretization point key
     *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
     *
     */
    bool onTetrahedronBoundary(const std::vector<int> & p, const int nInterval)
    {
          return ((p[0] == 0) || (p[1] == 0) || (p[2] == 0) || (p[0] + p[1] + p[2] == nInterval));
    }




    template <int mydim>
    void addCurvilinearSimplex(
            const Dune::GeometryType & geomtype,
            const std::vector<GlobalVector> & nodeSet,
            std::vector<int> & tagSet,
            int elementOrder,
            int nDiscretizationPoint,
            double shrinkMagnitude,
            double boundaryMagnification,
            bool interpolate,
            bool writeEdgeData,
            bool writeTriangleData
            )
    {
        typedef FieldVector< double, mydim >      LocalVector;

        CurvilinearElementInterpolator<double, mydim, dimension> elementInterpolator(geomtype, nodeSet, elementOrder);

        const int codim = dimension - mydim;


        // *******************************************************************************
        // Step 1. Find coordinates of the corners and the (linear) center of mass
        // *******************************************************************************
        int nDofThis    = nodeSet.size();
        int nCornerThis = elementInterpolator.nCorner();


        LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKWriter: Calculating CoM" );
        std::vector<GlobalVector> cornerVector;
        for (int i = 0; i < nCornerThis; i++) { cornerVector.push_back(elementInterpolator.corner(i)); }
        GlobalVector CoM = vectorCentreOfMass(cornerVector);


        // *******************************************************************************
        // Step 2: Construct sampling grid over the element. Sample the points, store them
        //   and map to their insertion coordinate
        // For a Generic triangle the sampling looks like this
        //
        //  (i,j)=(n,0) -> *
        //                 **
        //                 ***
        //                 ****
        //  (i,j)=(0,0) -> ****** <- (i,j) = (0,n)
        //
        // *******************************************************************************
        LocalCoordinate2GlobalIdMap parametricToIndex;
        int nInterval = interpolate ? nDiscretizationPoint - 1 : elementOrder;


        LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKWriter: Calculating Enumerators" );
        ElemGridEnumerate  simplexEnumerate        = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval);
        ElemGridEnumerate  simplexEnumerateReduced = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval-1);
        std::vector< LocalVector > simplexLocalGrid = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, mydim>(simplexEnumerate, nInterval);

        LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and inserting refinement vertices" );
        for (int i = 0; i < simplexEnumerate.size(); i++)
        {
            // Find if this vertex is internal or boundary
            bool isBoundaryPoint = (mydim == 3) ? onTetrahedronBoundary(simplexEnumerate[i], nInterval) : true;

            // Write this vertex only if we are going to write an element using it
            if (writeEdgeData || (writeTriangleData && isBoundaryPoint)) {
                // If we interpolate, then all points will be taken from new sample grid
                // Otherwise we take the intrinsic interpolation point grid which has the same shape
                GlobalVector tmpPoint = interpolate ? elementInterpolator.realCoordinate(simplexLocalGrid[i]) : nodeSet[i];
                for (int d = 0; d < dimension; d++)  {
                    tmpPoint[d] = (tmpPoint[d] + (CoM[d] - tmpPoint[d]) * shrinkMagnitude) * boundaryMagnification;
                }

                // Add coordinates to the coordinates array
                vtkPoint_.push_back(tmpPoint);

                // Add point to the point map
                parametricToIndex[simplexEnumerate[i]] = vtkPoint_.size() - 1;
            }
            //std::cout << "* coords " << parUV[0] << ", " << parUV[1] << ", in the map cooresponds to " << parametricToIndex[parUV] << std::endl;
        }

        // *******************************************************************************
        // Step 3: Store the index of vertices used for future use
        // *******************************************************************************
        tmpNInterval_ = nInterval;
        tmpParameter2Index_ = parametricToIndex;

        // *******************************************************************************
        // Step 4: Write edges that discretize this element to VTK
        // *******************************************************************************
        if (writeEdgeData)
        {
        	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and writing refinement-edges" );
        	SubEntityIndexVector thisEntitySubset = VTKEntitySubset::refineEntitySubset<codim, EDGE_CODIM>(simplexEnumerateReduced, nInterval, parametricToIndex);
        	for (int i = 0; i < thisEntitySubset.size(); i++) { addDiscretizationEdge(thisEntitySubset[i], tagSet); }
        }

        // *******************************************************************************
        // Step 5: Split this face into tiny triangles and add them to the triangle array
        // *******************************************************************************
        if (writeTriangleData)
        {
        	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and writing refinement-triangles" );
        	SubEntityIndexVector thisEntitySubset = VTKEntitySubset::refineEntitySubset<codim, FACE_CODIM>(simplexEnumerateReduced, nInterval, parametricToIndex);

        	for (int i = 0; i < thisEntitySubset.size(); i++) {
        		//std::cout << "subentity " << i << " consists of " << Dune::VectorHelper::vector2string(thisEntitySubset[i]) << std::endl;
        		addDiscretizationTriangle(thisEntitySubset[i], tagSet);
        	}
        }
    }


  private:
    // MPI
    int rank_;
    int size_;

    int tmpNInterval_;
    LocalCoordinate2GlobalIdMap tmpParameter2Index_;   // Stores the set of vertex indices used when inserting the last entity

    std::vector<GlobalVector> vtkPoint_;
    std::vector<IndexVector> vtkEdgeVertexIndex_;
    std::vector<IndexVector> vtkTriangleVertexIndex_;

    TagVector vtkEdgePhysicalTag_;
    TagVector vtkEdgeStructuralType_;
    TagVector vtkEdgeProcessRank_;

    TagVector vtkTrianglePhysicalTag_;
    TagVector vtkTriangleStructuralType_;
    TagVector vtkTriangleProcessRank_;

    // Field storage
    FieldNameMap                 fieldName2Index_;
    std::vector<FieldCoordMap>   vtkFieldVector_;   // For each field maps vertex index to field

  };



}  // namespace Dune

#endif /** DUNE_CURVILINEARVTKWRITER_HH **/
