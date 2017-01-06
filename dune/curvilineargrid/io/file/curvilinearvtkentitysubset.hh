#ifndef CURV_VTK_ENTITY_SUBSET
#define CURV_VTK_ENTITY_SUBSET

namespace Dune
{

namespace CurvGrid
{

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

} // Namespace CurvGrid

} // Namespace Dune
#endif // CURV_VTK_ENTITY_SUBSET
