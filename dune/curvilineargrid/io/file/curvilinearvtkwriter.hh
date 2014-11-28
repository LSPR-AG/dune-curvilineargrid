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

namespace Dune
{

  
  struct VtkEntityStructuralType
  {
    enum {
        Internal = 0,          // Elements internal to this process
        DomainBoundary = 1,    // Faces on the boundary of the computational domain
        ProcessBoundary = 2,   // Faces on the interprocessor boundary (not including overlap)
        Ghost = 3,             // Elements borrowed from the neighboring process
        Overlap = 4,           // Elements overlapping between processes
        Front = 5,             // Faces on the overlap boundary
        Periodic = 6           // Boundary faces which are periodic
    };
  };
  
  
  
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
      GlobalVector vectorCentreOfMass( const std::vector<GlobalVector> & cornerVector)
      {
    	  GlobalVector rez;
    	  for (int i = 0; i < cornerVector.size(); i++) { rez += cornerVector[i]; }
    	  rez /= cornerVector.size();
    	  return rez;
      }

      /** \brief Connects discretization points of a curvilinear triangle using linear edges, adds these edges for writing to VTK
       *
       *  \param[in]  triangleEnumeratorReduced      A vector of keys that enumerate triangular discretization points.
       *                                             Reduced means that enumerator's points-per-edge is one less than the number used to discretize the element.
       *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *  \param[in]  thisElmPhysTag                 Physical Tag of the element being discretized
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *
       */
      void addTriangularInterpolationEdgeSet(
    		  ElemGridEnumerate & triangleEnumeratorReduced,
    		  int nInterval,
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
       *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *  \param[in]  thisElmPhysTag                 Physical Tag of the element being discretized
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *
       */
      void addTetrahedralInterpolationEdgeSet(
    		  ElemGridEnumerate & tetrahedronEnumeratorReduced,
    		  int nInterval,
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
       *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *
       */
      bool onTetrahedronBoundary(const std::vector<int> & p, const int nInterval)
      {
  		  return ((p[0] == 0) || (p[1] == 0) || (p[2] == 0) || (p[0] + p[1] + p[2] == nInterval));
      }

      /** \brief Takes tetrahedral discretization point key to globalID 3D map, and splits it into 4 triangle 2D maps, corresponding to the faces of the tetrahedron
       *
       *  \param[in]  triangleEnumerator             A vector of keys that enumerate tetrahedral discretization points.
       *  \param[in]  nInterval                     Number of discretization intervals per edge (number of discretization points-per-edge minus 1)
       *  \param[in]  parametricToIndex              Map from discretization point enumerator key to that point's globalId (within this writer)
       *
       *  \return vector of 4 face triangle maps from 2D discretization point enumerator key to globalId of the point
       *
       */
      Coord2GlobalMapVector tetrahedralInterpolationFaceSet(
    		  ElemGridEnumerate & triangleEnumerator,
    		  int nInterval,
    		  LocalCoordinate2GlobalIdMap & parametricToIndex)
      {
    	  Coord2GlobalMapVector rez(4);

    	  for (int i = 0; i < triangleEnumerator.size(); i++)
    	  {
    		  int x = triangleEnumerator[i][0];
    		  int y = triangleEnumerator[i][1];
    		  int z = nInterval - x - y;

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
    	  return thisElmStructuralType == VtkEntityStructuralType::DomainBoundary;
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
    	  int thisElmCornerNo = thisElementInt.nCorner();

    	  std::vector<GlobalVector> cornerVector;
    	  for (int i = 0; i < thisElmCornerNo; i++) { cornerVector.push_back(thisElementInt.corner(i)); }
          GlobalVector CoM = vectorCentreOfMass(cornerVector);


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
          int nInterval = interpolate ? nDiscretizationPoints - 1 : thisElmOrder;

          // Expand all boundary surfaces a little bit so that they do not interlay with element surfaces
          double BOUNDARY_MAGNIFICATION = isBoundary(thisElmStructuralType) ? 1.2 : 1.0;


          ElemGridEnumerate  simplexEnumerate        = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval);
          ElemGridEnumerate  simplexEnumerateReduced = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<mydim>(nInterval-1);
          std::vector< LocalVector > simplexLocalGrid = Dune::CurvilinearGeometryHelper::simplexGridCoordinateSet<double, mydim>(simplexEnumerate, nInterval);

          for (int i = 0; i < simplexEnumerate.size(); i++)
          {
        	  // Find if this vertex is internal or boundary
        	  bool isBoundaryPoint = (mydim == 3) ? onTetrahedronBoundary(simplexEnumerate[i], nInterval) : true;

        	  // Write this vertex only if we are going to write an element using it
        	  if (writeEdgeData || (writeTriangleData && isBoundaryPoint)) {
            	  // If we interpolate, then all points will be taken from new sample grid
            	  // Otherwise we take the intrinsic interpolation point grid which has the same shape
            	  GlobalVector tmpPoint = interpolate ? thisElementInt.realCoordinate(simplexLocalGrid[i]) : elementNodeSet[i];
            	  for (int d = 0; d < cdim; d++)  {
            		  tmpPoint[d] = (tmpPoint[d] + (CoM[d] - tmpPoint[d]) * shrinkMagnitude) * BOUNDARY_MAGNIFICATION;
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
        	  case 2:  addTriangularInterpolationEdgeSet(simplexEnumerateReduced, nInterval, thisElmPhysTag, thisElmStructuralType, parametricToIndex);  break;
        	  case 3:  addTetrahedralInterpolationEdgeSet(simplexEnumerateReduced, nInterval, thisElmPhysTag, thisElmStructuralType, parametricToIndex);  break;
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
        	  		  ElemGridEnumerate triangleEnumerate = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<2>(nInterval);
        	  		  ElemGridEnumerate triangleEnumerateReduced = Dune::CurvilinearGeometryHelper::simplexGridEnumerate<2>(nInterval - 1);

        	          Coord2GlobalMapVector consistingTriangles = tetrahedralInterpolationFaceSet(triangleEnumerate, nInterval, parametricToIndex);

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
    void writeVTK( std::string filename)
    {
        FILE* vtkFile = fopen(filename.c_str(), "w");

        int nElements = vtkEdge.size() + vtkTriangle.size();
        int nCells = vtkEdge.size()*3 + vtkTriangle.size()*4;

        fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
        fprintf(vtkFile, "CurvilinearGmshReader test output\n");
        fprintf(vtkFile, "ASCII\n");
        fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(vtkFile, "POINTS %d double\n", vtkPoint.size() );

        // Write all points
        for (int i = 0; i < vtkPoint.size(); i++ ) {
        	for (int d = 0; d < cdim; d++)  { fprintf(vtkFile, "%lg ", vtkPoint[i][d]); }
        	fprintf(vtkFile, "\n");
        }

        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELLS %d %d\n", nElements, nCells );

        // Write all elements
        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtkFile, "2 %d %d \n", vtkEdge[i][0], vtkEdge[i][1]); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtkFile, "3 %d %d %d\n", vtkTriangle[i][0], vtkTriangle[i][1], vtkTriangle[i][2] ); }

        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELL_TYPES %d\n", nElements );

        // Write edge and triangle cell types
        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtkFile, "3\n"); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtkFile, "5\n"); }


        fprintf(vtkFile, "\n");
        fprintf(vtkFile, "CELL_DATA %d\n", nElements );

        // Write edge and triangle Structural type
        fprintf(vtkFile, "SCALARS physicalTag FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgeTag.size(); i++ )        { fprintf(vtkFile, "%d\n", vtkEdgeTag[i]); }
        for (int i = 0; i < vtkTriangleTag.size(); i++ )    { fprintf(vtkFile, "%d\n", vtkTriangleTag[i]); }

        // Write edge and triangle physicalTags
        fprintf(vtkFile, "SCALARS structuralType FLOAT\n");
        fprintf(vtkFile, "LOOKUP_TABLE default\n");
        for (int i = 0; i < vtkEdgeStructuralType.size(); i++ )        { fprintf(vtkFile, "%d\n", vtkEdgeStructuralType[i]); }
        for (int i = 0; i < vtkTriangleStructuralType.size(); i++ )    { fprintf(vtkFile, "%d\n", vtkTriangleStructuralType[i]); }

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
        fprintf(pvtuFile, "<VTKFile type=\"%s\" version=\"%s\" byte_order=\"%s\">\n", VTK_GRID_TYPE.c_str(), VTK_VTU_VERSION.c_str(), VTK_BYTE_ORDER.c_str());
        fprintf(pvtuFile, "<%s GhostLevel=\"0\">\n", VTK_GRID_TYPE.c_str());

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
        fprintf(pvtuFile, "</%s>\n", VTK_GRID_TYPE.c_str());
        fprintf(pvtuFile, "</VTKFile>\n");
        fclose(pvtuFile);
    }

    // Writes serial VTU file
    void writeVTU( std::string filename)
    {
        FILE* vtuFile = fopen(filename.c_str(), "w");

        int nVertices = vtkPoint.size();
        int nElements = vtkEdge.size() + vtkTriangle.size();

        // Compute offsets which is a general way to determine the number of vertices per element
        std::vector<int> offsets;
        int tmp_offset = 0;
        for (int i = 0; i < vtkEdge.size(); i++ )        { tmp_offset += 2;  offsets.push_back(tmp_offset); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { tmp_offset += 3;  offsets.push_back(tmp_offset); }

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
        for (int i = 0; i < vtkEdgeTag.size(); i++ )        { fprintf(vtuFile, "%d ", vtkEdgeTag[i]); }
        for (int i = 0; i < vtkTriangleTag.size(); i++ )    { fprintf(vtuFile, "%d ", vtkTriangleTag[i]); }

        fprintf(vtuFile, "\n</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" Name=\"structuralType\" NumberOfComponents=\"1\" format=\"ascii\">\n");

        // Write edge and triangle structural type
        for (int i = 0; i < vtkEdgeStructuralType.size(); i++ )        { fprintf(vtuFile, "%d ", vtkEdgeStructuralType[i]); }
        for (int i = 0; i < vtkTriangleStructuralType.size(); i++ )    { fprintf(vtuFile, "%d ", vtkTriangleStructuralType[i]); }


        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</CellData>\n");


        // Write coordinates of vertices
        // *****************************************************
        fprintf(vtuFile, "<Points>\n");
        fprintf(vtuFile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

        // Write all points
        for (int i = 0; i < vtkPoint.size(); i++ ) {
        	for (int d = 0; d < cdim; d++)  { fprintf(vtuFile, "%lg ", vtkPoint[i][d]); }
        	fprintf(vtuFile, "\n");
        }

        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "</Points>\n");


        // Write element information
        // *****************************************************
        fprintf(vtuFile, "<Cells>\n");
        fprintf(vtuFile, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");   // Vertex indices

        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtuFile, "%d %d ", vtkEdge[i][0], vtkEdge[i][1]); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtuFile, "%d %d %d ", vtkTriangle[i][0], vtkTriangle[i][1], vtkTriangle[i][2] ); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");

        for (int i = 0; i < offsets.size(); i++ )        { fprintf(vtuFile, "%d ", offsets[i]); }

        fprintf(vtuFile, "\n");
        fprintf(vtuFile, "</DataArray>\n");
        fprintf(vtuFile, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");   // Element types

        for (int i = 0; i < vtkEdge.size(); i++ )        { fprintf(vtuFile, "3 "); }
        for (int i = 0; i < vtkTriangle.size(); i++ )    { fprintf(vtuFile, "5 "); }

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
    void writeParallelVTU(std::string filenameBody, int rank, int size)
    {
    	// Write a PVTU file on master process
    	if (rank == 0) { writePVTU(filenameBody, size); }

    	// Write a VTU file on all processes
    	writeVTU(filenameBody  + "_process_" + std::to_string(rank) + ".vtu");
    }

  };



}  // namespace Dune

#endif /** DUNE_CURVILINEARVTKWRITER_HH **/
