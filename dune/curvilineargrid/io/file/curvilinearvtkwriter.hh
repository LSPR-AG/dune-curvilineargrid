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

#include <dune/curvilineargeometry/interpolation/lagrangeinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>



namespace Dune
{

namespace CurvGrid {

namespace VTKEntitySubset
{


typedef int                                         IndexType;
typedef std::vector<IndexType>                      IndexVector;
typedef std::vector< IndexVector >                  ElemGridEnumerate;

typedef std::map<std::vector<int>, int>             LocalCoordinate2GlobalIdMap;
typedef std::vector< LocalCoordinate2GlobalIdMap >  Coord2GlobalMapVector;

typedef std::vector<std::vector<int> >              SubEntityIndexVector;

// Codimensions of entity types for better code readability
const int   VERTEX_CODIM   = 3;
const int   EDGE_CODIM     = 2;
const int   FACE_CODIM     = 1;
const int   ELEMENT_CODIM  = 0;



/** \brief Splits provided triangular prism into 2 or 3 tetrahedrons, depending on whether the prism is reduced or not  */
SubEntityIndexVector prism2tetrahedron(const IndexVector & prismIndex) {
	assert((prismIndex.size() == 5) ||(prismIndex.size() == 6));  // Allow only triangular and reduced-triangular prisms

	SubEntityIndexVector rez;

	rez.push_back(IndexVector {prismIndex[0], prismIndex[1], prismIndex[2], prismIndex[4]}  );
	rez.push_back(IndexVector {prismIndex[0], prismIndex[2], prismIndex[3], prismIndex[4]}  );

	if (prismIndex.size() == 6)  // Add the 3rd tetrahedron only if the prism is not reduced
	{
		rez.push_back(IndexVector {prismIndex[2], prismIndex[3], prismIndex[4], prismIndex[5]}  );
	}

	return rez;
}


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
SubEntityIndexVector refineEntitySubset(const ElemGridEnumerate & entityEnumeratorReduced, int nInterval, const LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	std::cout << "ERROR: called refineErrorSubset generic (non-specialized) form that does nothing!" << std::endl;
	assert(0);
    return SubEntityIndexVector();
}


/** \brief Connects discretization points of a curvilinear edge using linear edges, adds these edges for writing to VTK */
template <>
SubEntityIndexVector refineEntitySubset<EDGE_CODIM, EDGE_CODIM>(
		const ElemGridEnumerate & edgeEnumeratorReduced,
        int nInterval,
        const LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	SubEntityIndexVector rez;

    // Construct all edges and add them to the edge array
    for (unsigned int i = 0; i < edgeEnumeratorReduced.size(); i++)
    {
    	IndexType x = edgeEnumeratorReduced[i][0];

        // Construct triangle (123)
        IndexVector parUV_0 {x};
        IndexVector parUV_1 {x + 1};

        IndexVector edge01 {parametricToIndex.at(parUV_0), parametricToIndex.at(parUV_1)};

        rez.push_back(edge01);
    }
    return rez;
}


/** \brief Connects discretization points of a curvilinear triangle using linear edges, adds these edges for writing to VTK */
template <>
SubEntityIndexVector refineEntitySubset<FACE_CODIM, EDGE_CODIM> (
		const ElemGridEnumerate & triangleEnumeratorReduced,
        int nInterval,
        const LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	SubEntityIndexVector rez;

    // Construct all edges and add them to the edge array
    for (unsigned int i = 0; i < triangleEnumeratorReduced.size(); i++)
    {
    	IndexType x = triangleEnumeratorReduced[i][0];
    	IndexType y = triangleEnumeratorReduced[i][1];

        // Construct triangle (123)
        IndexVector parUV_0 {x,      y,   };
        IndexVector parUV_1 {x + 1,  y,   };
        IndexVector parUV_2 {x,      y + 1};

        IndexVector edge01 { parametricToIndex.at(parUV_0), parametricToIndex.at(parUV_1) };
        IndexVector edge12 { parametricToIndex.at(parUV_1), parametricToIndex.at(parUV_2) };
        IndexVector edge20 { parametricToIndex.at(parUV_2), parametricToIndex.at(parUV_0) };

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
		const ElemGridEnumerate & tetrahedronEnumeratorReduced,
        int nInterval,
        const LocalCoordinate2GlobalIdMap & parametricToIndex)
{
	SubEntityIndexVector rez;

    // Construct all edges and add them to the edge array
    for (unsigned int i = 0; i < tetrahedronEnumeratorReduced.size(); i++)
    {
    	IndexType x = tetrahedronEnumeratorReduced[i][0];
    	IndexType y = tetrahedronEnumeratorReduced[i][1];
    	IndexType z = tetrahedronEnumeratorReduced[i][2];

        // Construct triangle (123)
        IndexVector parUVW_0 {x,      y,      z    };
        IndexVector parUVW_1 {x + 1,  y,      z    };
        IndexVector parUVW_2 {x,      y + 1,  z    };
        IndexVector parUVW_3 {x,      y,      z + 1};

        IndexVector edge01 { parametricToIndex.at(parUVW_0), parametricToIndex.at(parUVW_1) };
        IndexVector edge12 { parametricToIndex.at(parUVW_1), parametricToIndex.at(parUVW_2) };
        IndexVector edge20 { parametricToIndex.at(parUVW_2), parametricToIndex.at(parUVW_0) };
        IndexVector edge30 { parametricToIndex.at(parUVW_3), parametricToIndex.at(parUVW_0) };
        IndexVector edge31 { parametricToIndex.at(parUVW_3), parametricToIndex.at(parUVW_1) };
        IndexVector edge32 { parametricToIndex.at(parUVW_3), parametricToIndex.at(parUVW_2) };

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
		const ElemGridEnumerate & triangleEnumeratorReduced,
        int nInterval,
        const LocalCoordinate2GlobalIdMap & parametricToIndex
        )
{
	SubEntityIndexVector rez;

    for (unsigned int i = 0; i < triangleEnumeratorReduced.size(); i++) {
        // Construct two triangles (123) and (234) from 4 points, where 1 is the index point
        //   3--4
        //   |\ |
        //   | \|
        //   1--2
        // First triangle always exists, because we iterate over (i,j) in such a way that there is 1 free point at the edge
        // Second triangle we construct only if point 4 is still in the triangle

    	IndexType x = triangleEnumeratorReduced[i][0];
    	IndexType y = triangleEnumeratorReduced[i][1];

        // Construct triangle (123)
        IndexVector parUV_1 {x    , y    };
        IndexVector parUV_2 {x + 1, y    };
        IndexVector parUV_3 {x    , y + 1};

        IndexVector triangle123 { parametricToIndex.at(parUV_1), parametricToIndex.at(parUV_2), parametricToIndex.at(parUV_3) };

        // Add triangle (123)
        rez.push_back(triangle123);

        // Check if point 4 is still within the original triangle
        // If yes then add its symmetric triangle too
        if (x + y + 1 < nInterval)
        {
            IndexVector parUV_4 {x + 1, y + 1};
            IndexVector triangle234 { parametricToIndex.at(parUV_2), parametricToIndex.at(parUV_3), parametricToIndex.at(parUV_4) };

            // Add triangle (234)
            rez.push_back(triangle234);
        }
    }
    return rez;
}


/** \brief Connects discretization points of a curvilinear tetrahedron using linear triangles, adds these edges for writing to VTK */
template<>
SubEntityIndexVector refineEntitySubset<ELEMENT_CODIM, FACE_CODIM> (
		const ElemGridEnumerate & tetrahedralEnumeratorReduced,
        int nInterval,
        const LocalCoordinate2GlobalIdMap & parametricToIndex
        )
{
	SubEntityIndexVector rez;

    ElemGridEnumerate triangleEnumerator = CurvilinearGeometryHelper::simplexGridEnumerate<2>(nInterval);
    ElemGridEnumerate triangleEnumeratorReduced = CurvilinearGeometryHelper::simplexGridEnumerate<2>(nInterval - 1);

    // Takes tetrahedral discretization point key to globalID 3D map, and splits it into 4 triangle 2D maps, corresponding to the faces of the tetrahedron
    Coord2GlobalMapVector consistingTriangles(4);
    for (unsigned int i = 0; i < triangleEnumerator.size(); i++)
    {
    	IndexType x = triangleEnumerator[i][0];
    	IndexType y = triangleEnumerator[i][1];
    	IndexType z = nInterval - x - y;

        IndexVector faceInd_0 {x, y, 0};
        IndexVector faceInd_1 {x, 0, y};
        IndexVector faceInd_2 {0, x, y};
        IndexVector faceInd_3 {z, x, y};

        consistingTriangles[0][triangleEnumerator[i]] = parametricToIndex.at(faceInd_0);
        consistingTriangles[1][triangleEnumerator[i]] = parametricToIndex.at(faceInd_1);
        consistingTriangles[2][triangleEnumerator[i]] = parametricToIndex.at(faceInd_2);
        consistingTriangles[3][triangleEnumerator[i]] = parametricToIndex.at(faceInd_3);
    }

    // Discretizes resulting triangles using triangle discretization routine
    // [FIXME] Replace Nr 4 by subentity number
    for (unsigned int iFace = 0; iFace < 4; iFace++)
    {
    	SubEntityIndexVector faceSub = refineEntitySubset<FACE_CODIM, FACE_CODIM> (triangleEnumeratorReduced, nInterval, consistingTriangles[iFace]);
    	rez.insert (rez.end(), faceSub.begin(), faceSub.end());
    }

    return rez;
}


/** \brief Connects discretization points of a curvilinear tetrahedron using linear tetrahedrons  */
template<>
SubEntityIndexVector refineEntitySubset<ELEMENT_CODIM, ELEMENT_CODIM> (
		const ElemGridEnumerate & tetrahedralEnumeratorReduced,
        int nInterval,
        const LocalCoordinate2GlobalIdMap & parametricToIndex
        )
{
	SubEntityIndexVector rez;

    // For each discretization point, takes associated cube, splits it into 2 prisms, then each prism into 3 tets
    for (unsigned int i = 0; i < tetrahedralEnumeratorReduced.size(); i++)
    {
    	IndexType x = tetrahedralEnumeratorReduced[i][0];
    	IndexType y = tetrahedralEnumeratorReduced[i][1];
    	IndexType z = tetrahedralEnumeratorReduced[i][2];

        // If this is the boundary tetrahedron, it will not have any opposite tetrahedrons
        if (x + y + z + 1 == nInterval)
        {
        	rez.push_back(
                IndexVector
                {
        			parametricToIndex.at( IndexVector{x,     y,     z}     ),
        			parametricToIndex.at( IndexVector{x + 1, y,     z}     ),
        			parametricToIndex.at( IndexVector{x,     y + 1, z}     ),
        			parametricToIndex.at( IndexVector{x,     y,     z + 1} )
                }
        	);
        } else {  // Otherwise, all opposite tetrahedrons need to be enumerated as well

            std::vector<IndexVector > cubeIndex   // Enumerate all vertices of a cube associated with this discretization point
            {
            	{x,     y,     z},
            	{x + 1, y,     z},
            	{x,     y + 1, z},
            	{x + 1, y + 1, z},
            	{x,     y,     z + 1},
            	{x + 1, y,     z + 1},
            	{x,     y + 1, z + 1},
            	{x + 1, y + 1, z + 1}
            };

            IndexVector prismIndex1  {            // Split the cube into 2 prisms
            	parametricToIndex.at(cubeIndex[0]),
                parametricToIndex.at(cubeIndex[1]),
                parametricToIndex.at(cubeIndex[2]),
                parametricToIndex.at(cubeIndex[4]),
                parametricToIndex.at(cubeIndex[5]),
                parametricToIndex.at(cubeIndex[6])
            };

            IndexVector prismIndex2  {            // This prism might only have 5 points if the last corner of the cube is not inside original tetrahedron
            	parametricToIndex.at(cubeIndex[1]),
                parametricToIndex.at(cubeIndex[2]),
                parametricToIndex.at(cubeIndex[3]),
                parametricToIndex.at(cubeIndex[5]),
                parametricToIndex.at(cubeIndex[6])
            };

            // If the top-most corner of the cube exists, then add it as well to the 2nd prism
            // Otherwise the 2nd prism is a reduced prism with only 5 corners
            if (x + y + z + 2 < nInterval) { prismIndex2.push_back(parametricToIndex.at(cubeIndex[7])); }

            SubEntityIndexVector prism2tetIndex1 = prism2tetrahedron(prismIndex1);
            SubEntityIndexVector prism2tetIndex2 = prism2tetrahedron(prismIndex2);

            rez.insert (rez.end(), prism2tetIndex1.begin(), prism2tetIndex1.end());
            rez.insert (rez.end(), prism2tetIndex2.begin(), prism2tetIndex2.end());
        }
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

	  typedef typename GridType::ctype  ctype;
      static const int dimension = GridType::dimension;
      static const bool isCached = GridType::is_cached;

      // Codimensions of entity types for better code readability
      static const int   VERTEX_CODIM   = 3;
      static const int   EDGE_CODIM     = 2;
      static const int   FACE_CODIM     = 1;
      static const int   ELEMENT_CODIM  = 0;

      typedef FieldVector< ctype, dimension >             GlobalCoordinate;
      typedef std::vector<GlobalCoordinate>               GlobalCoordinateVec;
      typedef std::vector<int>                            IndexVector;
      typedef std::vector<int>                            TagVector;
      typedef std::vector<IndexVector >                   SubEntityIndexVector;
      typedef std::map<std::vector<int>, int>             LocalCoordinate2GlobalIdMap;
      typedef std::vector< LocalCoordinate2GlobalIdMap >  Coord2GlobalMapVector;
      typedef typename LocalCoordinate2GlobalIdMap::iterator  LC2GIMapIter;

      typedef std::map<std::string, unsigned int>         FieldNameMap;
      typedef std::pair<std::string, unsigned int>        FieldNamePair;
      typedef std::map<int, ctype>                        FieldScalarMap;
      typedef std::pair<int, ctype>                       FieldScalarPair;
      typedef std::map<int, GlobalCoordinate>             FieldCoordMap;
      typedef std::pair<int, GlobalCoordinate>            FieldCoordPair;

      typedef typename FieldNameMap::iterator             FieldNameMapIter;
      typedef typename FieldNameMap::const_iterator       FieldNameMapConstIter;
      typedef typename FieldScalarMap::iterator           FieldScalarMapIter;
      typedef typename FieldScalarMap::const_iterator     FieldScalarMapConstIter;
      typedef typename FieldCoordMap::iterator            FieldCoordMapIter;
      typedef typename FieldCoordMap::const_iterator      FieldCoordMapConstIter;

      typedef std::vector< IndexVector >                  ElemGridEnumerate;

      // VTK constants
      const std::string VTK_XML_VERSION = "1.0";
      const std::string VTK_GRID_TYPE = "UnstructuredGrid";
      const std::string VTK_VTU_VERSION = "0.1";
      const std::string VTK_BYTE_ORDER = "LittleEndian";


  public:

    CurvilinearVTKWriter (MPIHelper &mpihelper)
    {
    	nInterval_ = 0;         // Defined later
        rank_ = mpihelper.rank();
        size_ = mpihelper.size();
    }


    CurvilinearVTKWriter (int rank, int size)
    {
    	nInterval_ = 0;         // Defined later
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
     *  \param[in]  interpolate                    If set to false, instead of interpolating, the provided interpolation points will be used as discretization points. Otherwise, a sub-sampling grid over the element will be interpolated
     *  \param[in]  explode	                    there will be gaps between all elements, obtained by scaling them away from the center of mass
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
            const GlobalCoordinateVec & nodeSet,
            const TagVector & tagSet,
            int elementOrder,
            int nDiscretizationPoint,
            bool interpolate,
            bool explode,
            std::vector<bool> writeCodim)
    {
    	assert((mydim > 0) && (mydim <= 3));  // Forbid writing vertices and hypergeometric entities

        int   thisElmPhysTag        = (tagSet.size() > 0 ) ? tagSet[0] : 0;
        int   thisElmPartitionType  = (tagSet.size() > 1 ) ? tagSet[1] : 0;
        int   thisElmProcessRank    = (tagSet.size() > 2 ) ? tagSet[2] : 0;

        // Treat boundary segments in a special way
        std::string pname;
        bool magnify;
        if (thisElmPartitionType == BOUNDARY_SEGMENT_PARTITION_TYPE) {
        	magnify = true;
        	pname = "BoundarySegment";
        } else if (thisElmPartitionType == INTERIOR_BOUNDARY_SEGMENT_PARTITION_TYPE) {
            	magnify = false;
            	pname = "InteriorBoundarySegment";
        } else {
        	magnify = false;
        	pname = Dune::PartitionName(static_cast<Dune::PartitionType> (thisElmPartitionType));
        }

        // It is not possible to write entities of dimension higher than that of the discretized entity
        assert(writeCodim.size() >= dimension);  // Ensure that there enough info provided
        if (mydim == 2)  { writeCodim[ELEMENT_CODIM] = false; }
        if (mydim == 1)  { writeCodim[ELEMENT_CODIM] = false;  writeCodim[FACE_CODIM] = false; }

        std::stringstream log_message;
        log_message << "VTK_WRITER: Adding a curvilinear element Type=" << CurvilinearGeometryHelper::geometryName(geomtype);
        log_message << " Order="               << elementOrder;
        log_message << " PhysicalTag="         << thisElmPhysTag;
        log_message << " StructuralType="      << pname;
        log_message << " ProcessRank="         << thisElmProcessRank;
        log_message << " nDiscretization="     << nDiscretizationPoint;
        log_message << " useInterpolation="    << interpolate;
        log_message << " explodeElements="     << explode;
        log_message << " writeVTK_edges="      << writeCodim[EDGE_CODIM];
        log_message << " writeVTK_triangles="  << writeCodim[FACE_CODIM];
        log_message << " writeVTK_triangles="  << writeCodim[ELEMENT_CODIM];
        //log_message << " vertices=" << VectorHelper::vector2string(nodeSet);
        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, log_message.str());

        addCurvilinearSimplex<mydim>(geomtype, nodeSet, tagSet, elementOrder, nDiscretizationPoint, explode, magnify, interpolate, writeCodim);
    }


    // Initialize field name. Useful, if master process does not have any entities of a given type, but still needs to write a correct .pvtu file
    template <class VTKVectorFunction>
    void initScalarField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Initializing scalar field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(scalarFieldName2Index_, vtkFieldScalar_, vtkfunction.name());
    }

    // Initialize field name. Useful, if master process does not have any entities of a given type, but still needs to write a correct .pvtu file
    template <class VTKVectorFunction>
    void initVectorField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Initializing vector field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(vectorFieldName2Index_, vtkFieldVector_, vtkfunction.name());
    }



    template <class VTKVectorFunction>
    void addScalarField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Adding element field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(scalarFieldName2Index_, vtkFieldScalar_, vtkfunction.name());

    	// Loop over all points
    	for (LC2GIMapIter iter = parameter2Index_.begin(); iter != parameter2Index_.end(); iter++)
    	{
    		std::vector<int> integerCoordinate = (*iter).first;
    		unsigned int     vertexIndex       = (*iter).second;

    		// Compute local coordinate
    		Domain local;
    		for (int i = 0; i < mydimension; i++)  { local[i] = (ctype((*iter).first[i])) / nInterval_; }

    		// Evaluate field
    		Range field = vtkfunction.evaluate(local);

    		//std::cout << "VTKWriter sampling local coordinate " << local << " coordinate index " << vertexIndex << " field=" << field << std::endl;

    		// Append coordinate index and field
    		vtkFieldScalar_[fieldIndex].insert(FieldScalarPair(vertexIndex, field));
    	}
    }

    template <class VTKVectorFunction>
    void addVectorField(const VTKVectorFunction & vtkfunction) {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKGridWriter: Adding element field " + vtkfunction.name());

    	typedef typename VTKVectorFunction::Domain  Domain;
    	typedef typename VTKVectorFunction::Range   Range;
    	const int mydimension = VTKVectorFunction::mydimension;

        // The index of the field associated with this field name
    	int fieldIndex = updateFieldIndex(vectorFieldName2Index_, vtkFieldVector_, vtkfunction.name());

    	// Loop over all points
    	for (LC2GIMapIter iter = parameter2Index_.begin(); iter != parameter2Index_.end(); iter++)
    	{
    		std::vector<int> integerCoordinate = (*iter).first;
    		unsigned int     vertexIndex       = (*iter).second;

    		// Compute local coordinate
    		Domain local;
    		for (int i = 0; i < mydimension; i++)  { local[i] = (ctype((*iter).first[i])) / nInterval_; }

    		// Evaluate field
    		Range field = vtkfunction.evaluate(local);

    		//std::cout << "VTKWriter sampling local coordinate " << local << " coordinate index " << vertexIndex << " field=" << field << std::endl;

    		// Append coordinate index and field
    		vtkFieldVector_[fieldIndex].insert(FieldCoordPair(vertexIndex, field));
    	}
    }


    // Writes serial VTK file
    // Takes a vector of vertices, a vector of edges and a vector of triangles, writes them all to a .VTK file
    void writeVTK(std::string filename) const
    {
        FILE* vtkFile = fopen(filename.c_str(), "w");

        std::vector<std::size_t> nEntity
        {
        	vtkCodimVertexIndex_[ELEMENT_CODIM].size(),
        	vtkCodimVertexIndex_[FACE_CODIM].size(),
        	vtkCodimVertexIndex_[EDGE_CODIM].size(),
        	vtkPoint_.size()
        };

        unsigned int nEntityTot =   nEntity[EDGE_CODIM] +   nEntity[FACE_CODIM] +   nEntity[ELEMENT_CODIM];
        unsigned int nCellTot   = 3*nEntity[EDGE_CODIM] + 4*nEntity[FACE_CODIM] + 5*nEntity[ELEMENT_CODIM];

        // Check if the datasize is consistent for writing
        writerSelfCheck(nEntity);

        fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
        fprintf(vtkFile, "CurvilinearGmshReader test output\n");
        fprintf(vtkFile, "ASCII\n");
        fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");


        // Write all points
        fprintf(vtkFile, "POINTS %d double\n", vtkPoint_.size() );
        for (unsigned int i = 0; i < nEntity[VERTEX_CODIM]; i++ ) {
            for (int d = 0; d < dimension; d++)  { fprintf(vtkFile, "%lg ", vtkPoint_[i][d]); }
            fprintf(vtkFile, "\n");
        }


        // Write all elements
        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELLS %d %d\n", nEntityTot, nCellTot );
        for (int iCodim = EDGE_CODIM; iCodim >= ELEMENT_CODIM; iCodim--)  // Write edges first, elements last
        {
            for (unsigned int i = 0; i < nEntity[iCodim]; i++)     {
            	unsigned int nSubVertex = vtkCodimVertexIndex_[iCodim][i].size();
				fprintf(vtkFile, "%d ", nSubVertex);

            	for (unsigned int iV = 0; iV < nSubVertex; iV++)  { fprintf(vtkFile, "%d ", vtkCodimVertexIndex_[iCodim][i][iV]); }

            	fprintf(vtkFile, "\n");
            }
        }


        // Write edge and triangle cell types
        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELL_TYPES %d\n", nEntityTot );
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++ )     { fprintf(vtkFile, "3\n"); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++ )     { fprintf(vtkFile, "5\n"); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++ )  { fprintf(vtkFile, "10\n"); }


        // If are defined fields associated with vertices, then write them too
        if ((vectorFieldName2Index_.size() > 0) || (scalarFieldName2Index_.size() > 0))
        {
        	fprintf(vtkFile, "\n");
        	fprintf(vtkFile, "POINT_DATA %d\n", vtkPoint_.size() );

        	// Write all scalar fields
        	for (FieldNameMapConstIter iterName = scalarFieldName2Index_.begin(); iterName != scalarFieldName2Index_.end(); ++iterName)
        	{
        		std::string fieldName = (*iterName).first;
        		int fieldIndex  = (*iterName).second;

        		fprintf(vtkFile, "SCALARS %s FLOAT\n", fieldName.c_str());

                for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                	// If there is no field defined for this vertex, just print a zero vector for consistency
                	FieldScalarMapIter iterField = vtkFieldScalar_[fieldIndex].find(i);
                	if (iterField == vtkFieldScalar_[fieldIndex].end())  { fprintf(vtkFile, "0.0\n"); }
                	else {                                                 fprintf(vtkFile, "%lg\n", (*iterField).second); }
                }
                fprintf(vtkFile, "\n");
        	}

        	// Write all vector fields
        	for (FieldNameMapConstIter iterName = vectorFieldName2Index_.begin(); iterName != vectorFieldName2Index_.end(); ++iterName)
        	{
        		std::string fieldName = (*iterName).first;
        		int fieldIndex  = (*iterName).second;

        		fprintf(vtkFile, "VECTORS %s FLOAT\n", fieldName.c_str());

                for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                	// If there is no field defined for this vertex, just print a zero vector for consistency
                	FieldCoordMapIter iterField = vtkFieldVector_[fieldIndex].find(i);
                	if (iterField == vtkFieldVector_[fieldIndex].end())  { fprintf(vtkFile, "0.0 0.0 0.0\n"); }
                	else {
                		GlobalCoordinate v = (*iterField).second;
                		fprintf(vtkFile, "%lg %lg %lg\n", v[0], v[1], v[2]);
                	}
                }
                fprintf(vtkFile, "\n");
        	}
        }

        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELL_DATA %d\n", nEntityTot );


        // Write edge and triangle Structural type
        fprintf(vtkFile, "SCALARS physicalTag FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { fprintf(vtkFile, "%d\n", vtkCodimPhysicalTag_[EDGE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { fprintf(vtkFile, "%d\n", vtkCodimPhysicalTag_[FACE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { fprintf(vtkFile, "%d\n", vtkCodimPhysicalTag_[ELEMENT_CODIM][i]); }


        // Write edge and triangle physicalTags
        fprintf(vtkFile, "SCALARS structuralType FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { fprintf(vtkFile, "%d\n", vtkCodimStructuralType_[EDGE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { fprintf(vtkFile, "%d\n", vtkCodimStructuralType_[FACE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { fprintf(vtkFile, "%d\n", vtkCodimStructuralType_[ELEMENT_CODIM][i]); }


        // Write edge and triangle provider process ranks
        fprintf(vtkFile, "SCALARS processRank FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { fprintf(vtkFile, "%d\n", vtkCodimProcessRank_[EDGE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { fprintf(vtkFile, "%d\n", vtkCodimProcessRank_[FACE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { fprintf(vtkFile, "%d\n", vtkCodimProcessRank_[ELEMENT_CODIM][i]); }


        // Empty line at the end of file
        fprintf(vtkFile, "\n");
        fclose(vtkFile);
    }


    // Writes a PVTU parallel file (no data in this file)
    // Requite path and filename as separate strings, because
    void writePVTU(std::string path, std::string filenameBody, int size) const
    {
    	std::string filename = path + filenameBody + ".pvtu";
        FILE* pvtuFile = fopen(filename.c_str(), "w");


        // Write header
        // *****************************************************
        fprintf(pvtuFile, "<?xml version=\"%s\"?>\n", VTK_XML_VERSION.c_str());
        fprintf(pvtuFile, "<VTKFile type=\"P%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(pvtuFile, "<P%s GhostLevel=\"0\">\n", VTK_GRID_TYPE.c_str());


        // Write PointData - scalar and vector fields sampled over all vertices
        // *****************************************************
        bool haveScalar = scalarFieldName2Index_.size() > 0;
        bool haveVector = vectorFieldName2Index_.size() > 0;
        if (haveScalar || haveVector)
        {
        	// Generate and write PointData preamble
            std::stringstream pointDataNames;
            if (haveScalar)  { pointDataNames << "Scalars=\"" << namesPile(scalarFieldName2Index_) << "\"";  if (haveVector) { pointDataNames << " "; } }
            if (haveVector)  { pointDataNames << "Vectors=\"" << namesPile(vectorFieldName2Index_) << "\"";                                             }
        	fprintf(pvtuFile, "<PPointData %s>\n", pointDataNames.str().c_str());

        	// Write point data arrays
        	for (FieldNameMapConstIter iterName = scalarFieldName2Index_.begin(); iterName != scalarFieldName2Index_.end(); ++iterName)
        	{
        		fprintf(pvtuFile, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"ascii\"/>\n", (*iterName).first.c_str());
        	}

        	for (FieldNameMapConstIter iterName = vectorFieldName2Index_.begin(); iterName != vectorFieldName2Index_.end(); ++iterName)
        	{
        		fprintf(pvtuFile, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"/>\n", (*iterName).first.c_str());
        	}

            fprintf(pvtuFile, "</PPointData>\n");
        }


        // Write scalars associated with entities
        //*****************************************************
        fprintf(pvtuFile, "<PCellData Scalars=\"physicalTag\">\n");
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" Name=\"physicalTag\" NumberOfComponents=\"1\"/>\n");      // Write physicalTags
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\"/>\n");   // Write structural types
        fprintf(pvtuFile, "<PDataArray type=\"Float32\" Name=\"processRank\" NumberOfComponents=\"1\"/>\n");      // Write process ranks
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
    void writeVTU(std::string filename) const
    {
        FILE* vtuFile = fopen(filename.c_str(), "w");

        std::vector<std::size_t> nEntity
        {
        	vtkCodimVertexIndex_[ELEMENT_CODIM].size(),
        	vtkCodimVertexIndex_[FACE_CODIM].size(),
        	vtkCodimVertexIndex_[EDGE_CODIM].size(),
        	vtkPoint_.size()
        };

        // Check if the datasize is consistent for writing
        writerSelfCheck(nEntity);

        // Compute offsets which is a general way to determine the number of vertices per element
        std::vector<unsigned int> offsets;
        int tmp_offset = 0;
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { tmp_offset += 2;  offsets.push_back(tmp_offset); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { tmp_offset += 3;  offsets.push_back(tmp_offset); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { tmp_offset += 4;  offsets.push_back(tmp_offset); }

        // Write header
        // *****************************************************
        fprintf(vtuFile, "<?xml version=\"%s\"?>\n", VTK_XML_VERSION.c_str());
        fprintf(vtuFile, "<VTKFile type=\"%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(vtuFile, "<%s>\n", VTK_GRID_TYPE.c_str());
        fprintf(vtuFile, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nEntity[VERTEX_CODIM], nEntity[EDGE_CODIM] + nEntity[FACE_CODIM] + nEntity[ELEMENT_CODIM]);


        // Write PointData - scalar and vector fields sampled over all vertices
        // *****************************************************
        bool haveScalar = scalarFieldName2Index_.size() > 0;
        bool haveVector = vectorFieldName2Index_.size() > 0;
        if (haveScalar || haveVector)
        {
        	// Generate and write PointData preamble
            std::stringstream pointDataNames;
            if (haveScalar)  { pointDataNames << "Scalars=\"" << namesPile(scalarFieldName2Index_) << "\"";  if (haveVector) { pointDataNames << " "; } }
            if (haveVector)  { pointDataNames << "Vectors=\"" << namesPile(vectorFieldName2Index_) << "\"";                                             }
        	fprintf(vtuFile, "<PointData %s>\n", pointDataNames.str().c_str());

        	// Write point data arrays
        	for (FieldNameMapConstIter iterName = scalarFieldName2Index_.begin(); iterName != scalarFieldName2Index_.end(); ++iterName)
        	{
        		std::string fieldName = (*iterName).first;
        		int fieldIndex  = (*iterName).second;

        		fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"ascii\">\n", fieldName.c_str());

                for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                	// If there is no field defined for this vertex, just print a zero vector for consistency
                	FieldScalarMapConstIter iterField = vtkFieldScalar_[fieldIndex].find(i);
                	if (iterField == vtkFieldScalar_[fieldIndex].end())  { fprintf(vtuFile, "0.0\n"); }
                	else {                                                 fprintf(vtuFile, "%lg\n", (*iterField).second); }
                }

        		fprintf(vtuFile, "</DataArray>\n");
        	}

        	for (FieldNameMapConstIter iterName = vectorFieldName2Index_.begin(); iterName != vectorFieldName2Index_.end(); ++iterName)
        	{
        		std::string fieldName = (*iterName).first;
        		int fieldIndex  = (*iterName).second;

        		fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", fieldName.c_str());

                for (unsigned int i = 0; i < vtkPoint_.size(); i++ ) {
                	// If there is no field defined for this vertex, just print a zero vector for consistency
                	FieldCoordMapConstIter iterField = vtkFieldVector_[fieldIndex].find(i);
                	if (iterField == vtkFieldVector_[fieldIndex].end())  { fprintf(vtuFile, "0.0 0.0 0.0\n"); }
                	else {
                		GlobalCoordinate v = (*iterField).second;
                		fprintf(vtuFile, "%lg %lg %lg\n", v[0], v[1], v[2]);
                	}
                }

        		fprintf(vtuFile, "</DataArray>\n");
        	}

            fprintf(vtuFile, "</PointData>\n");
        }


        // Write edge and triangle physicalTags and structural types
        // *****************************************************
        fprintf(vtuFile, "<CellData Scalars=\"physicalTag\">\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"physicalTag\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle physicalTags
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { fprintf(vtuFile, "%d ", vtkCodimPhysicalTag_[EDGE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { fprintf(vtuFile, "%d ", vtkCodimPhysicalTag_[FACE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { fprintf(vtuFile, "%d ", vtkCodimPhysicalTag_[ELEMENT_CODIM][i]); }

        fprintf(vtuFile, "\n</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle structural type
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { fprintf(vtuFile, "%d ", vtkCodimStructuralType_[EDGE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { fprintf(vtuFile, "%d ", vtkCodimStructuralType_[FACE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { fprintf(vtuFile, "%d ", vtkCodimStructuralType_[ELEMENT_CODIM][i]); }

        fprintf(vtuFile, "\n</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"processRank\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle provider process ranks
        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++)     { fprintf(vtuFile, "%d ", vtkCodimProcessRank_[EDGE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++)     { fprintf(vtuFile, "%d ", vtkCodimProcessRank_[FACE_CODIM][i]); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++)  { fprintf(vtuFile, "%d ", vtkCodimProcessRank_[ELEMENT_CODIM][i]); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</CellData>\n");


        // Write coordinates of vertices
        // *****************************************************
        fprintf(vtuFile, "<Points>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

        // Write all points
        for (unsigned int i = 0; i < nEntity[VERTEX_CODIM]; i++ ) {
            for (int d = 0; d < dimension; d++)  { fprintf(vtuFile, "%lg ", vtkPoint_[i][d]); }
            fprintf(vtuFile, "\n");
        }

        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</Points>\n");


        // Write element information
        // *****************************************************
        fprintf(vtuFile, "<Cells>\n");
        fprintf(vtuFile, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");   // Vertex indices

        for (int iCodim = EDGE_CODIM; iCodim >= ELEMENT_CODIM; iCodim--)  // Write edges first, elements last
        {
            for (unsigned int i = 0; i < nEntity[iCodim]; i++)     {
            	for (unsigned int iV = 0; iV < vtkCodimVertexIndex_[iCodim][i].size(); iV++)  { fprintf(vtuFile, "%d ", vtkCodimVertexIndex_[iCodim][i][iV]); }
            	fprintf(vtuFile, "\n");
            }
        }


        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

        for (unsigned int i = 0; i < offsets.size(); i++ )        { fprintf(vtuFile, "%d ", offsets[i]); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");   // Element types

        for (unsigned int i = 0; i < nEntity[EDGE_CODIM]; i++ )     { fprintf(vtuFile, "3 "); }
        for (unsigned int i = 0; i < nEntity[FACE_CODIM]; i++ )     { fprintf(vtuFile, "5 "); }
        for (unsigned int i = 0; i < nEntity[ELEMENT_CODIM]; i++ )  { fprintf(vtuFile, "10 "); }

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
    void writeParallelVTU(std::string path, std::string filenameBody) const
    {
        // Write a PVTU file on master process
        if (rank_ == 0) { writePVTU(path, filenameBody, size_); }

        // Write a VTU file on all processes
        writeVTU(path + filenameBody  + "_process_" + std::to_string(rank_) + ".vtu");
    }


    // Checks if the storage arrays are consistent before writing to file
    void writerSelfCheck(std::vector<std::size_t> nEntity) const
    {
        for (int iCodim = ELEMENT_CODIM; iCodim <= EDGE_CODIM; iCodim++)
        {
            assert(vtkCodimPhysicalTag_[iCodim].size()    == nEntity[iCodim]);
            assert(vtkCodimStructuralType_[iCodim].size() == nEntity[iCodim]);
            assert(vtkCodimProcessRank_[iCodim].size()    == nEntity[iCodim]);
        }
    }

  protected:

    // ***********************************************
    // Auxiliary Methods
    // ***********************************************

    /** \brief Adds a set of linear discretization entities which will be explicitly written to the file
     *  Note: Only add tags if they are provided by the user
     * */
    template<int codim, int subcodim>
    void refineEntity(const ElemGridEnumerate & simplexEnumerateReduced, int nInterval, const LocalCoordinate2GlobalIdMap & parametricToIndex, const std::vector<int> & tagSet)
    {
    	LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and writing refinement-edges" );
    	SubEntityIndexVector thisEntitySubset = VTKEntitySubset::refineEntitySubset<codim, subcodim>(simplexEnumerateReduced, nInterval, parametricToIndex);

    	for (unsigned int i = 0; i < thisEntitySubset.size(); i++) {
            vtkCodimVertexIndex_[subcodim].push_back(thisEntitySubset[i]);
            vtkCodimPhysicalTag_[subcodim].push_back(    (tagSet.size() > 0) ? tagSet[0] : 0);
            vtkCodimStructuralType_[subcodim].push_back( (tagSet.size() > 1) ? tagSet[1] : 0);
            vtkCodimProcessRank_[subcodim].push_back(    (tagSet.size() > 2) ? tagSet[2] : 0);
    	}
    }

    /** \brief Calculates the centre of mass of a vector of points (equal-weighted) */
    GlobalCoordinate vectorCentreOfMass( const GlobalCoordinateVec & cornerVector) const
    {
        GlobalCoordinate rez;
        for (unsigned int i = 0; i < cornerVector.size(); i++) { rez += cornerVector[i]; }
        rez /= cornerVector.size();
        return rez;
    }


    /** \brief Checks if the tetrahedral discretization point corresponds to tetrahedral boundary
     *
     *  \param[in]  p                               tetrahedral discretization point key
     *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
     *
     */
    bool onTetrahedronBoundary(const std::vector<int> & p, const int nInterval) const
    {
          return ((p[0] == 0) || (p[1] == 0) || (p[2] == 0) || (p[0] + p[1] + p[2] == nInterval));
    }


    // Pass CoM, because it is more optimal to only calculate it once for all samples over the element
    void shrinkMagnify(GlobalCoordinate & point, const GlobalCoordinate & CoM, bool explode, bool magnify) const
    {
        ctype shrinkMagnitude       = explode ? 0.2 : 0.0;  // 0.0 - no shrinking, 0.99 - very small element (must be < 1)
        ctype boundaryMagnification = magnify ? 1.2 : 1.0;  // Expand all domain boundary surfaces a little bit so that they do not interlay with element surfaces

        for (int d = 0; d < dimension; d++)  {
        	point[d] = (point[d] + (CoM[d] - point[d]) * shrinkMagnitude) * boundaryMagnification;
        }
    }


    template <int mydim>
    void addCurvilinearSimplex(
            const Dune::GeometryType & geomtype,
            const GlobalCoordinateVec & nodeSet,
            const std::vector<int> & tagSet,
            int elementOrder,
            int nDiscretizationPoint,
            bool shrink,
            bool magnify,
            bool interpolate,
            std::vector<bool> writeCodim)
    {

        typedef FieldVector< ctype, mydim >      LocalVector;

        LagrangeInterpolator<ctype, mydim, dimension> elementInterpolator(geomtype, nodeSet, elementOrder);

        const int codim = dimension - mydim;


        // *******************************************************************************
        // Step 1. Find coordinates of the corners and the (linear) center of mass
        // *******************************************************************************
        int nDofThis    = nodeSet.size();
        int nCornerThis = elementInterpolator.nCorner();


        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Calculating CoM" );
        GlobalCoordinateVec cornerVector;
        for (int i = 0; i < nCornerThis; i++) { cornerVector.push_back(elementInterpolator.corner(i)); }
        GlobalCoordinate CoM = vectorCentreOfMass(cornerVector);


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
        // Step 3: Store the index of vertices used for future use
        // *******************************************************************************
        nInterval_ = interpolate ? nDiscretizationPoint - 1 : elementOrder;
        parameter2Index_.clear();


        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Calculating Enumerators" );
        ElemGridEnumerate  simplexEnumerate        = CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval_);
        ElemGridEnumerate  simplexEnumerateReduced = CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval_ - 1);
        std::vector< LocalVector > simplexLocalGrid = CurvilinearGeometryHelper::simplexGridCoordinateSet<ctype, mydim>(simplexEnumerate, nInterval_);

        LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, "CurvilinearVTKWriter: Computing and inserting refinement vertices" );
        for (unsigned int i = 0; i < simplexEnumerate.size(); i++)
        {
            // Find if this vertex is internal or boundary
            bool isBoundaryPoint = (mydim == 3) ? onTetrahedronBoundary(simplexEnumerate[i], nInterval_) : true;

            // If we interpolate, then all points will be taken from new sample grid
            // Otherwise we take the intrinsic interpolation point grid which has the same shape
            GlobalCoordinate tmpPoint = interpolate ? elementInterpolator.global(simplexLocalGrid[i]) : nodeSet[i];

            //Optionally, transform the element. Do not try to multipy element by 1.0 if it is not going to be transformed
            if (shrink || magnify) { shrinkMagnify(tmpPoint, CoM, shrink, magnify); }

            // Add coordinates to the coordinates array
            vtkPoint_.push_back(tmpPoint);

            // Add point to the point map
            parameter2Index_[simplexEnumerate[i]] = vtkPoint_.size() - 1;
        }

        // *******************************************************************************
        // Step 4: Write entities that discretize this element to VTK
        // Current paradigm is to either discretize all inserted entities with edges, making sort of a net
        // Or to discretize each inserted entity with smaller entities of the same codim
        // *******************************************************************************
        if      (writeCodim[EDGE_CODIM])  { refineEntity<codim, EDGE_CODIM>(simplexEnumerateReduced, nInterval_, parameter2Index_, tagSet);  }
        else if (writeCodim[codim])       { refineEntity<codim, codim>(simplexEnumerateReduced, nInterval_, parameter2Index_, tagSet);  }
        //if (writeCodim[FACE_CODIM])     { refineEntity<codim, FACE_CODIM>(simplexEnumerateReduced, nInterval, parametricToIndex, tagSet);  }
        //if (writeCodim[ELEMENT_CODIM])  { refineEntity<codim, ELEMENT_CODIM>(simplexEnumerateReduced, nInterval, parametricToIndex, tagSet);  }
    }


    template <class Map, class Storage>
    unsigned int updateFieldIndex(Map & namemap, Storage & storage, std::string name) {
    	FieldNameMapIter iter = namemap.find(name);

    	if (iter != namemap.end() )  {       // If field exists, find its index
    		return (*iter).second;
    	} else {                                      // Otherwise, make new index and add field name to them map
    		unsigned int rez = storage.size();
    		namemap.insert(FieldNamePair(name, rez));
    		storage.push_back(typename Storage::value_type());
    		return rez;
    	}
    }


    // Stacks up all names from the map in a single string, separated by commas, with no extra spaces
    std::string namesPile(const FieldNameMap & map) const {
    	std::stringstream rez;

    	int i = 0;
    	for (FieldNameMapConstIter iterName = map.begin(); iterName != map.end(); ++iterName) {
    		if (i > 0) { rez << ","; }
    		rez << (*iterName).first;
    		i++;
     	}

    	return rez.str();
    }




  private:
    // MPI
    int rank_;
    int size_;

    int nInterval_;
    LocalCoordinate2GlobalIdMap parameter2Index_;   // Stores the set of vertex indices used when inserting the last entity

    GlobalCoordinateVec vtkPoint_;

    // Discretization entity storage
    std::vector<IndexVector> vtkCodimVertexIndex_[4];   // Vertices that discretize the entity of a given codimension
    TagVector vtkCodimPhysicalTag_[4];                  // Physical tag of the entity
    TagVector vtkCodimStructuralType_[4];               // Structural type of the entity
    TagVector vtkCodimProcessRank_[4];                  // Process rank of the entity

    // Field storage
    FieldNameMap                 scalarFieldName2Index_;
    FieldNameMap                 vectorFieldName2Index_;
    std::vector<FieldScalarMap>  vtkFieldScalar_;   // For each field maps vertex index to field
    std::vector<FieldCoordMap>   vtkFieldVector_;   // For each field maps vertex index to field



  };


} // namespace CurvGrid

}  // namespace Dune

#endif /** DUNE_CURVILINEARVTKWRITER_HH **/
