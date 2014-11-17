// curvilinear gmsh reader functionality - in this implementation of the curvilinear gmsh reader we ONLY read higher order triangles and tetrahedra


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

#include <dune/geometry/polynomialinterpolation/curvilinearelementinterpolator.hh>

namespace Dune
{

  template<int cdim>
  class CurvilinearVTKWriter
  {

/** \brief This class takes curved elements, samples them on a grid, connects points into mesh of
 * small straight-sided triangles and edges, ands writes them to the .vtk file
 *
 */
  public:

	  const std::string VTK_XML_VERSION = "1.0";
	  const std::string VTK_GRID_TYPE = "UnstructuredGrid";
	  const std::string VTK_VTU_VERSION = "0.1";
	  const std::string VTK_BYTE_ORDER = "LittleEndian";

	  enum VTK_ENTITY_STRUCTURAL_TYPE {
	   	ENTITY_INTERNAL = 0,			// Elements internal to this process
	   	ENTITY_DOMAIN_BOUNDARY = 1,	// Faces on the boundary of the computational domain
	   	ENTITY_PROCESS_BOUNDARY = 2,   // Faces on the interprocessor boundary (not including overlap)
	   	ENTITY_GHOST = 3,              // Elements borrowed from the neighboring process
	   	ENTITY_OVERLAP = 4,			// Elements overlapping between processes
	   	ENTITY_FRONT = 5,				// Faces on the overlap boundary
	   	ENTITY_PERIODIC_BOUNDARY = 6   // Boundary faces which are periodic
	  };


  protected:
      typedef FieldVector< double, cdim >      GlobalVector;
      typedef std::vector<int>                 IndexVector;
      typedef std::vector<int>                 TagVector;
      typedef std::map<std::vector<int>, int>  LocalCoordinate2GlobalIdMap;
      typedef std::vector< LocalCoordinate2GlobalIdMap >  Coord2GlobalMapVector;

      typedef std::vector< IndexVector >       ElemGridEnumerate;

      std::vector<GlobalVector> vtkPoint;
      std::vector<IndexVector> vtkEdge;
      std::vector<IndexVector> vtkTriangle;

      TagVector vtkEdgeTag;
      TagVector vtkTriangleTag;

      TagVector vtkEdgeStructuralType;
      TagVector vtkTriangleStructuralType;

      /** \brief Adds an edge which will be explicitly written to the file */
      void addDiscretizationEdge(IndexVector indices, int tag, int strType)
      {
    	  vtkEdge.push_back(indices);
    	  vtkEdgeTag.push_back(tag);
    	  vtkEdgeStructuralType.push_back(strType);
      }

      /** \brief Adds a triangle which will be explicitly written to the file */
      void addDiscretizationTriangle(IndexVector indices, int tag, int strType)
      {
    	  vtkTriangle.push_back(indices);
    	  vtkTriangleTag.push_back(tag);
    	  vtkTriangleStructuralType.push_back(strType);
      }

      /** \brief Calculates the centre of mass of a vector of points (equal-weighted) */
      GlobalVector vectorCentreOfMass( const std::vector<GlobalVector> & corners)
      {
    	  GlobalVector rez;
    	  for (int i = 0; i < corners.size(); i++) { rez += corners[i]; }
    	  rez /= corners.size();
    	  return rez;
      }

      /** \brief Connects discretization points of a curvilinear triangle using linear edges, adds these edges for writing to VTK
       *
       *  \param[in]  triangleEnumeratorReduced      A vector of keys that enumerate triangular discretization points.
       *                                             Reduced means that enumerator's points-per-edge is one less than the number used to discretize the element.
       *  \param[in]  nIntervals                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *  \param[in]  thisElmPhysTag                 Physical Tag of the element being discretized
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *
       */
      void addTriangularInterpolationEdgeSet(
    		  ElemGridEnumerate & triangleEnumeratorReduced,
    		  int nIntervals,
    		  int thisElmPhysTag,
    		  int thisElmStructuralType,
    		  LocalCoordinate2GlobalIdMap & parametricToIndex)
      {
          // Construct all edges and add them to the edge array
    	  for (int i = 0; i < triangleEnumeratorReduced.size(); i++)
    	  {
    		  int x = triangleEnumeratorReduced[i][0];
    		  int y = triangleEnumeratorReduced[i][1];

              // Construct triangle (123)
              std::vector<int> parUV_0 {x,      y,   };
              std::vector<int> parUV_1 {x + 1,  y,   };
              std::vector<int> parUV_2 {x,      y + 1};

              std::vector<int> edge01;    edge01.push_back(parametricToIndex[parUV_0]);    edge01.push_back(parametricToIndex[parUV_1]);
              std::vector<int> edge12;    edge12.push_back(parametricToIndex[parUV_1]);    edge12.push_back(parametricToIndex[parUV_2]);
              std::vector<int> edge20;    edge20.push_back(parametricToIndex[parUV_2]);    edge20.push_back(parametricToIndex[parUV_0]);

              // Add all edges
              addDiscretizationEdge(edge01, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge12, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge20, thisElmPhysTag, thisElmStructuralType);
          }
      }

      /** \brief Connects discretization points of a curvilinear tetrahedron using linear edges, adds these edges for writing to VTK
       *
       *  \param[in]  tetrahedronEnumeratorReduced   A vector of keys that enumerate tetrahedral discretization points.
       *                                             Reduced means that enumerator's points-per-edge is one less than the number used to discretize the element.
       *  \param[in]  nIntervals                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *  \param[in]  thisElmPhysTag                 Physical Tag of the element being discretized
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *
       */
      void addTetrahedralInterpolationEdgeSet(
    		  ElemGridEnumerate & tetrahedronEnumeratorReduced,
    		  int nIntervals,
    		  int thisElmPhysTag,
    		  int thisElmStructuralType,
    		  LocalCoordinate2GlobalIdMap & parametricToIndex)
      {
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

              std::vector<int> edge01;    edge01.push_back(parametricToIndex[parUVW_0]);    edge01.push_back(parametricToIndex[parUVW_1]);
              std::vector<int> edge12;    edge12.push_back(parametricToIndex[parUVW_1]);    edge12.push_back(parametricToIndex[parUVW_2]);
              std::vector<int> edge20;    edge20.push_back(parametricToIndex[parUVW_2]);    edge20.push_back(parametricToIndex[parUVW_0]);
              std::vector<int> edge30;    edge30.push_back(parametricToIndex[parUVW_3]);    edge30.push_back(parametricToIndex[parUVW_0]);
              std::vector<int> edge31;    edge31.push_back(parametricToIndex[parUVW_3]);    edge31.push_back(parametricToIndex[parUVW_1]);
              std::vector<int> edge32;    edge32.push_back(parametricToIndex[parUVW_3]);    edge32.push_back(parametricToIndex[parUVW_2]);

              // Add all edges
              addDiscretizationEdge(edge01, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge12, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge20, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge30, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge31, thisElmPhysTag, thisElmStructuralType);
              addDiscretizationEdge(edge32, thisElmPhysTag, thisElmStructuralType);
          }
      }

      /** \brief Connects discretization points of a curvilinear triangle using linear triangles, adds these edges for writing to VTK
       *
       *  \param[in]  triangleEnumeratorReduced      A vector of keys that enumerate triangular discretization points.
       *                                             Reduced means that enumerator's points-per-edge is one less than the number used to discretize the element.
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *  \param[in]  thisElmPhysTag                 Physical Tag of the element being discretized
       *  \param[in]  nDiscretizationPoints          Number of discretization points-per-edge
       *
       */
      void addTriangularInterpolationTriangleSet(
    		  ElemGridEnumerate & triangleEnumeratorReduced,
    		  LocalCoordinate2GlobalIdMap & parametricToIndex,
    		  int thisElmPhysTag,
    		  int thisElmStructuralType,
    		  int nDiscretizationPoints
    		  )
      {

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

              std::vector<int> triangle123;
              triangle123.push_back(parametricToIndex[parUV_1]);
              triangle123.push_back(parametricToIndex[parUV_2]);
              triangle123.push_back(parametricToIndex[parUV_3]);

              //std::cout << "- coords " << parUV_1[0] << ", " << parUV_1[1] << ", in the map cooresponds to " << parametricToIndex[parUV_1] << std::endl;
              //std::cout << "- coords " << parUV_2[0] << ", " << parUV_2[1] << ", in the map cooresponds to " << parametricToIndex[parUV_2] << std::endl;
              //std::cout << "- coords " << parUV_3[0] << ", " << parUV_3[1] << ", in the map cooresponds to " << parametricToIndex[parUV_3] << std::endl;

              // Add triangle (123)
              addDiscretizationTriangle(triangle123, thisElmPhysTag, thisElmStructuralType);

              // Check if point 4 is still in the triangle
              if (x + y + 2 < nDiscretizationPoints)
              {
            	  std::vector<int> parUV_4 {x + 1, y + 1};

                  std::vector<int> triangle234;
                  triangle234.push_back(parametricToIndex[parUV_2]);
                  triangle234.push_back(parametricToIndex[parUV_3]);
                  triangle234.push_back(parametricToIndex[parUV_4]);

                  // Add triangle (234)
                  addDiscretizationTriangle(triangle234, thisElmPhysTag, thisElmStructuralType);
              }
          }

      }

      /** \brief Checks if the tetrahedral discretization point corresponds to tetrahedral boundary
       *
       *  \param[in]  p      						 tetrahedral discretization point key
       *  \param[in]  nIntervals                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *
       */
      bool onTetrahedronBoundary(const std::vector<int> & p, const int nIntervals)
      {
  		  return ((p[0] == 0) || (p[1] == 0) || (p[2] == 0) || (p[0] + p[1] + p[2] == nIntervals));
      }

      /** \brief Takes tetrahedral discretization point key to globalID 3D map, and splits it into 4 triangle 2D maps, corresponding to the faces of the tetrahedron
       *
       *  \param[in]  triangleEnumerator             A vector of keys that enumerate tetrahedral discretization points.
       *  \param[in]  nIntervals                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *
       *  \return vector of 4 face triangle maps from 2D discretization point enumerator key to globalId of the point
       *
       */
      Coord2GlobalMapVector tetrahedralInterpolationFaceSet(
    		  ElemGridEnumerate & triangleEnumerator,
    		  int nIntervals,
    		  LocalCoordinate2GlobalIdMap & parametricToIndex)
      {
    	  Coord2GlobalMapVector rez(4);

    	  for (int i = 0; i < triangleEnumerator.size(); i++)
    	  {
    		  int x = triangleEnumerator[i][0];
    		  int y = triangleEnumerator[i][1];
    		  int z = nIntervals - x - y;

              std::vector<int> faceInd_0 {x, y, 0};
              std::vector<int> faceInd_1 {x, 0, y};
              std::vector<int> faceInd_2 {0, x, y};
              std::vector<int> faceInd_3 {z, x, y};

              rez[0][triangleEnumerator[i]] = parametricToIndex[faceInd_0];
              rez[1][triangleEnumerator[i]] = parametricToIndex[faceInd_1];
              rez[2][triangleEnumerator[i]] = parametricToIndex[faceInd_2];
              rez[3][triangleEnumerator[i]] = parametricToIndex[faceInd_3];
    	  }
    	  return rez;
      }


      bool isBoundary(int thisElmStructuralType)
      {
    	  return thisElmStructuralType == ENTITY_DOMAIN_BOUNDARY;
      }

      // Test Function - Takes a triangle with all its interpolatory points, splits it into visualisation triangles
      // and adds all of them to vtk_triangle
      template <int mydim>
      void addCurvilinearSimplex(
              const Dune::GeometryType & thisElmType,
              const std::vector<GlobalVector> & elementNodeSet,
              int thisElmOrder,
              int thisElmPhysTag,
              int nDiscretizationPoints,
              double shrinkMagnitude,
              bool interpolate,
              bool writeEdgeData,
              bool writeTriangleData,
              int thisElmStructuralType
              )
      {
    	  typedef FieldVector< double, mydim >      LocalVector;

    	  CurvilinearElementInterpolator<double, mydim, cdim> thisElementInt(thisElmType, elementNodeSet, thisElmOrder);




    	  // *******************************************************************************
          // Step 1. Find coordinates of the corners and the (linear) center of mass
    	  // *******************************************************************************
    	  int thisElmDofNo = elementNodeSet.size();
    	  int thisElmCornerNo = thisElementInt.corners();

    	  std::vector<GlobalVector> corners;
    	  for (int i = 0; i < thisElementInt.corners(); i++) { corners.push_back(thisElementInt.corner(i)); }
          GlobalVector CoM = vectorCentreOfMass(corners);


          //std::cout << "Element with init dof " << elementDofs[0] << std::endl;



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
          int nIntervals = interpolate ? nDiscretizationPoints - 1 : thisElmOrder;

          // Expand all boundary surfaces a little bit so that they do not interlay with element surfaces
          double boundary_magnification = isBoundary(thisElmStructuralType) ? 1.2 : 1.0;

          ElemGridEnumerate  simplexEnumerate        = thisElementInt.simplexGridEnumerate(nIntervals);
          ElemGridEnumerate  simplexEnumerateReduced = thisElementInt.simplexGridEnumerate(nIntervals-1);
          std::vector< LocalVector > simplexLocalGrid = thisElementInt.simplexGridCoordinates(simplexEnumerate, nIntervals);

          for (int i = 0; i < simplexEnumerate.size(); i++)
          {
        	  // Find if this vertex is internal or boundary
        	  bool boundary_point = (mydim == 3) ? onTetrahedronBoundary(simplexEnumerate[i], nIntervals) : true;

        	  // Write this vertex only if we are going to write an element using it
        	  if (writeEdgeData || (writeTriangleData && boundary_point)) {
            	  // If we interpolate, then all points will be taken from new sample grid
            	  // Otherwise we take the intrinsic interpolation point grid which has the same shape
            	  GlobalVector tmpPoint = interpolate ? thisElementInt.realCoordinate(simplexLocalGrid[i]) : elementNodeSet[i];
            	  for (int d = 0; d < cdim; d++)  {
            		  tmpPoint[d] = (tmpPoint[d] + (CoM[d] - tmpPoint[d]) * shrinkMagnitude) * boundary_magnification;
            	  }

            	  // Add coordinates to the coordinates array
            	  vtkPoint.push_back(tmpPoint);

            	  // Add point to the point map
            	  parametricToIndex[simplexEnumerate[i]] = vtkPoint.size() - 1;
              }
              //std::cout << "* coords " << parUV[0] << ", " << parUV[1] << ", in the map cooresponds to " << parametricToIndex[parUV] << std::endl;
          }


          // *******************************************************************************
          // Step 3: Write edges discretizing this element to VTK
          // *******************************************************************************
          if (writeEdgeData)
          {
        	  switch (mydim)
        	  {
        	  case 2:  addTriangularInterpolationEdgeSet(simplexEnumerateReduced, nIntervals, thisElmPhysTag, thisElmStructuralType, parametricToIndex);  break;
        	  case 3:  addTetrahedralInterpolationEdgeSet(simplexEnumerateReduced, nIntervals, thisElmPhysTag, thisElmStructuralType, parametricToIndex);  break;
        	  }
          }

          // *******************************************************************************
          // Step 4: Split this face into tiny triangles and add them to the triangle array
          // *******************************************************************************

          if (writeTriangleData)
          {
              switch (mydim)
              {
        	  	  case 2:  addTriangularInterpolationTriangleSet(simplexEnumerateReduced, parametricToIndex, thisElmPhysTag, thisElmStructuralType, nDiscretizationPoints);  break;
        	  	  case 3:
        	  	  {
        	  		  ElemGridEnumerate triangleEnumerate = CurvilinearElementInterpolator<double, 2, 3>::simplexGridEnumerate(nIntervals);
        	  		  ElemGridEnumerate triangleEnumerateReduced = CurvilinearElementInterpolator<double, 2, 3>::simplexGridEnumerate(nIntervals - 1);

        	          Coord2GlobalMapVector consistingTriangles = tetrahedralInterpolationFaceSet(triangleEnumerate, nIntervals, parametricToIndex);

        	          for (int iFace = 0; iFace < 4; iFace++)
        	          {
        	        	  addTriangularInterpolationTriangleSet(triangleEnumerateReduced, consistingTriangles[iFace], thisElmPhysTag, thisElmStructuralType, nDiscretizationPoints);
        	          }
        	  	  }  break;
              }
          }
      }


  public:

    CurvilinearVTKWriter () { }

    /** \brief Takes curvilinear element, discretizes it into linear element, adds linear elements to the writer
     *
     *  \param[in]  thisElmType             	   Geometry Type of the element
     *  \param[in]  elementNodeSet                   Coordinates of the element's interpolation points
     *  \param[in]  thisElmOrder                   Curvilinear interpolation order of the element
     *  \param[in]  thisElmPhysTag                 Physical Tag of the element
     *  \param[in]  nDiscretizationPoints          Number of discretization points-per-edge
     *  \param[in]  interpolate			           If set to false, instead of interpolating, the provided interpolation points will be used as discretization points
     *  \param[in]  interpolate			           If set to true, there will be gaps between all elements, obtained by scaling them away from the center of mass
     *  \param[in]  writeEdgeData          		   If true, discretization points are connected using linear edges forming triangular mesh, and added to writer
     *  \param[in]  writeTriangleData       		   If true, discretization points are connected using linear triangles forming triangular mesh, and added to writer
     *  \param[in]  thisElmStructuralType  		   Defines the structural type of the element (internal, boundary, ghost, etc )
     *
     *	\note Discretization is performed using a regular grid over the reference element, for example,
     *	below the stars represent the triangular discretization points with nDiscretizationPoints = 5
     *
     *  *
     *  **
     *  ***
     *  ****
     *  *****
     */
    void addCurvilinearElement(
    		const Dune::GeometryType & thisElmType,
    		const std::vector<GlobalVector> & elementNodeSet,
            int thisElmOrder,
            int thisElmPhysTag,
            int nDiscretizationPoints,
            bool interpolate,
            bool explode,
            bool writeEdgeData,
            bool writeTriangleData,
            int thisElmStructuralType
    		)
    {
    	// 0.0 - no shrinking, 0.99 - very small element (must be < 1)
    	double shrinkMagnitude = explode ? 0.2 : 0.0;

    	if (thisElmType.isTriangle())
    	{
    		addCurvilinearSimplex<2>(thisElmType, elementNodeSet, thisElmOrder, thisElmPhysTag, nDiscretizationPoints, shrinkMagnitude, interpolate, writeEdgeData, writeTriangleData, thisElmStructuralType);
    	}
    	else if (thisElmType.isTetrahedron())
    	{
    		addCurvilinearSimplex<3>(thisElmType, elementNodeSet, thisElmOrder, thisElmPhysTag, nDiscretizationPoints, shrinkMagnitude, interpolate, writeEdgeData, writeTriangleData, thisElmStructuralType);
    	}
    	else
    	{
    		DUNE_THROW(Dune::IOError, "CURVILINEAR_VTK_WRITER: only implemented elements are triangles and tetrahedra at the moment");
    	}
    }


    // Writes serial VTK file
    // Takes a vector of vertices, a vector of edges and a vector of triangles, writes them all to a .VTK file
    void VTK_Write( std::string file_name)
    {
        FILE* vtk_file = fopen(file_name.c_str(), "w");

        int total_elements = vtkEdge.size() + vtkTriangle.size();
        int total_cell_entries = vtkEdge.size()*3 + vtkTriangle.size()*4;

        fprintf(vtk_file, "# vtk DataFile Version 2.0\n");
        fprintf(vtk_file, "CurvilinearGmshReader test output\n");
        fprintf(vtk_file, "ASCII\n");
        fprintf(vtk_file, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(vtk_file, "POINTS %d double\n", vtkPoint.size() );

        // Write all points
        for (int i = 0; i < vtkPoint.size(); i++ ) {
        	for (int d = 0; d < cdim; d++)  { fprintf(vtk_file, "%lg ", vtkPoint[i][d]); }
        	fprintf(vtk_file, "\n");
        }

        fprintf(vtk_file, "\n");
        fprintf(vtk_file, "CELLS %d %d\n", total_elements, total_cell_entries );

        // Write all elements
        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtk_file, "2 %d %d \n", vtkEdge[i][0], vtkEdge[i][1]); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtk_file, "3 %d %d %d\n", vtkTriangle[i][0], vtkTriangle[i][1], vtkTriangle[i][2] ); }

        fprintf(vtk_file, "\n");
        fprintf(vtk_file, "CELL_TYPES %d\n", total_elements );

        // Write edge and triangle cell types
        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtk_file, "3\n"); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtk_file, "5\n"); }


        fprintf(vtk_file, "\n");
        fprintf(vtk_file, "CELL_DATA %d\n", total_elements );

        // Write edge and triangle Structural type
        fprintf(vtk_file, "SCALARS physicalTag FLOAT\n");
        fprintf(vtk_file, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgeTag.size(); i++ )        { fprintf(vtk_file, "%d\n", vtkEdgeTag[i]); }
        for (int i = 0; i < vtkTriangleTag.size(); i++ )    { fprintf(vtk_file, "%d\n", vtkTriangleTag[i]); }

        // Write edge and triangle physicalTags
        fprintf(vtk_file, "SCALARS structuralType FLOAT\n");
        fprintf(vtk_file, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgeStructuralType.size(); i++ )        { fprintf(vtk_file, "%d\n", vtkEdgeStructuralType[i]); }
        for (int i = 0; i < vtkTriangleStructuralType.size(); i++ )    { fprintf(vtk_file, "%d\n", vtkTriangleStructuralType[i]); }

        // Empty line at the end of file
        fprintf(vtk_file, "\n");
        fclose(vtk_file);
    }

    // Writes a PVTU parallel file (no data in this file)
    void PVTU_Write(std::string file_name_body, int size)
    {
        std::string file_name = file_name_body + ".pvtu";
    	FILE* pvtu_file = fopen(file_name.c_str(), "w");

        // Write header
        // *****************************************************
        fprintf(pvtu_file, "<?xml version=\"%s\"?>\n", VTK_XML_VERSION.c_str());
        fprintf(pvtu_file, "<VTKFile type=\"%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(pvtu_file, "<%s GhostLevel=\"0\">\n", VTK_GRID_TYPE.c_str());

        // PointData could be provided here
        // (For example, tags associated with vertices, which we do not have at the moment)
        // *****************************************************
        // <PPointData>...</PPointData>


        // Write edge and triangle physicalTags
        // *****************************************************
        fprintf(pvtu_file, "<PCellData Scalars=\"physicalTag\">\n");
        fprintf(pvtu_file, "<PDataArray type=\"Float32\" Name=\"physicalTag\" NumberOfComponents=\"1\"/>\n");

        // Write edge and triangle structural types
        // *****************************************************
        fprintf(pvtu_file, "<PDataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\"/>\n");
        fprintf(pvtu_file, "</PCellData>\n");


        // Write coordinates of vertices
        // *****************************************************
        fprintf(pvtu_file, "<PPoints>\n");
        fprintf(pvtu_file, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
        fprintf(pvtu_file, "</PPoints>\n");


        // Write element information
        // *****************************************************
        fprintf(pvtu_file, "<PCells>\n");
        fprintf(pvtu_file, "<PDataArray type=\"Int32\" Name=\"connectivity\"/>\n");   // Vertex indices
        fprintf(pvtu_file, "<PDataArray type=\"Int32\" Name=\"offsets\"/>\n");
        fprintf(pvtu_file, "<PDataArray type=\"UInt8\" Name=\"types\"/>\n");   // Element types
        fprintf(pvtu_file, "</PCells>\n");


        // Write all .vtu data files file
        // *****************************************************
        for (int iProc = 0; iProc < size; iProc++ )
        {
        	std::string vtu_file_name = file_name_body + "_process_" + std::to_string(iProc) + ".vtu";
        	fprintf(pvtu_file, "<Piece  Source=\"%s\"/>\n", vtu_file_name.c_str());
        }


        // Finish writing file
        // *****************************************************
        fprintf(pvtu_file, "</%s>\n", VTK_GRID_TYPE.c_str());
        fprintf(pvtu_file, "</VTKFile>\n");
        fclose(pvtu_file);
    }

    // Writes serial VTU file
    void VTU_Write( std::string file_name)
    {
        FILE* vtu_file = fopen(file_name.c_str(), "w");

        int total_vertices = vtkPoint.size();
        int total_elements = vtkEdge.size() + vtkTriangle.size();

        // Compute offsets which is a general way to determine the number of vertices per element
        std::vector<int> offsets;
        int tmp_offset = 0;
        for (int i = 0; i < vtkEdge.size(); i++ )        { tmp_offset += 2;  offsets.push_back(tmp_offset); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { tmp_offset += 3;  offsets.push_back(tmp_offset); }

        // Write header
        // *****************************************************
        fprintf(vtu_file, "<?xml version=\"%s\"?>\n", VTK_XML_VERSION.c_str());
        fprintf(vtu_file, "<VTKFile type=\"%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(vtu_file, "<%s>\n", VTK_GRID_TYPE);
        fprintf(vtu_file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", total_vertices, total_elements);


        // PointData could be provided here
        // (For example, tags associated with vertices, which we do not have at the moment)
        // *****************************************************
        // <PointData>...</PointData>


        // Write edge and triangle physicalTags and structural types
        // *****************************************************
        fprintf(vtu_file, "<CellData Scalars=\"physicalTag\">\n");
        fprintf(vtu_file, "<DataArray type=\"Float32\" Name=\"physicalTag\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle physicalTags
        for (int i = 0; i < vtkEdgeTag.size(); i++ )        { fprintf(vtu_file, "%d ", vtkEdgeTag[i]); }
        for (int i = 0; i < vtkTriangleTag.size(); i++ )    { fprintf(vtu_file, "%d ", vtkTriangleTag[i]); }

        fprintf(vtu_file, "\n</DataArray>\n");
        fprintf(vtu_file, "<DataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle structural type
        for (int i = 0; i < vtkEdgeStructuralType.size(); i++ )        { fprintf(vtu_file, "%d ", vtkEdgeStructuralType[i]); }
        for (int i = 0; i < vtkTriangleStructuralType.size(); i++ )    { fprintf(vtu_file, "%d ", vtkTriangleStructuralType[i]); }


        fprintf(vtu_file, "\n");
        fprintf(vtu_file, "</DataArray>\n");
        fprintf(vtu_file, "</CellData>\n");


        // Write coordinates of vertices
        // *****************************************************
        fprintf(vtu_file, "<Points>\n");
        fprintf(vtu_file, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

        // Write all points
        for (int i = 0; i < vtkPoint.size(); i++ ) {
        	for (int d = 0; d < cdim; d++)  { fprintf(vtu_file, "%lg ", vtkPoint[i][d]); }
        	fprintf(vtu_file, "\n");
        }

        fprintf(vtu_file, "</DataArray>\n");
        fprintf(vtu_file, "</Points>\n");


        // Write element information
        // *****************************************************
        fprintf(vtu_file, "<Cells>\n");
        fprintf(vtu_file, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");   // Vertex indices

        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtu_file, "%d %d ", vtkEdge[i][0], vtkEdge[i][1]); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtu_file, "%d %d %d ", vtkTriangle[i][0], vtkTriangle[i][1], vtkTriangle[i][2] ); }

        fprintf(vtu_file, "\n");
        fprintf(vtu_file, "</DataArray>\n");
        fprintf(vtu_file, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

        for (int i = 0; i < offsets.size(); i++ )        { fprintf(vtu_file, "%d ", offsets[i]); }

        fprintf(vtu_file, "\n");
        fprintf(vtu_file, "</DataArray>\n");
        fprintf(vtu_file, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");   // Element types

        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtu_file, "3 "); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtu_file, "5 "); }

        fprintf(vtu_file, "\n");
        fprintf(vtu_file, "</DataArray>\n");
        fprintf(vtu_file, "</Cells>\n");


        // Finish writing file
        // *****************************************************
        fprintf(vtu_file, "</Piece>\n");
        fprintf(vtu_file, "</%s>\n", VTK_GRID_TYPE);
        fprintf(vtu_file, "</VTKFile>\n");
        fclose(vtu_file);
    }

    // Writes a VTU file on all processes and a PVTU on Master Process
    void VTU_ParallelWrite(std::string file_name_body, int rank, int size)
    {
    	// Write a PVTU file on master process
    	if (rank == 0) { PVTU_Write(file_name_body, size); }

    	// Write a VTU file on all processes
    	VTU_Write(file_name_body  + "_process_" + std::to_string(rank) + ".vtu");
    }

  };



}  // namespace Dune

#endif /** DUNE_CURVILINEARVTKWRITER_HH **/
