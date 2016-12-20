/***************************************************************************
                          tetmesh.h  -  description
                             -------------------
    begin                : Mon Dec 15 2003
    copyright            : (C) 2003 by Roman Geus
    email                : roman.geus@psi.ch
    edited by            : Hua Guo, Sep 2010
***************************************************************************/

/***************************************************************************
                          curvilineargridbase.hh
                             -------------------
    begin                : Tue Nov 25 2014
    copyright            : (C) 2014 by Aleksejs Fomins, LSPR AG
    description          : Upgraded the mesh to curvilinear grid
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DUNE_CURVILINEAROCTREENODE_HH
#define DUNE_CURVILINEAROCTREENODE_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>







namespace Dune {

namespace CurvGrid {


// Forward-declaration of GridBase because the modules include each other
template <class ct, int cdim, bool isCached>
class CurvilinearGridBase;


/** \brief Wraps CurvilinearGridBase element with functions necessary for CurvilinearOctree    */

template <class ct, int cdim, bool isCached>
class CurvilinearOctreeNode {
public:

    /* public types */
    typedef Dune::FieldVector<ct, cdim>      GlobalCoordinate;
    typedef std::vector<GlobalCoordinate>              VertexVector;

    typedef CurvilinearGridBase<ct, cdim, isCached>  GridBaseType;

    typedef typename GridBaseType::EntityStorage                        EntityStorage;
    typedef typename GridBaseType::template Codim<0>::EntityGeometry    ElementGeometry;



public: /* public methods */

    CurvilinearOctreeNode(const GridBaseType & grid, int elementIndex) :
    	gridbase_(grid),
    	elementIndex_(elementIndex)
	{
    	// Calculate it once since it will be asked for multiple times
    	calculateBoundingBox();
	}

    // Note: This operation is expensive - do not use too frequently
    ElementGeometry elementGeometry() { return gridbase_.template entityGeometry<0>(elementIndex_); }

    // Gets a box in which this Tetrahedron fits
    void elementBoundingBox(GlobalCoordinate & center, GlobalCoordinate & extent) const {
    	center = boundingBoxCenter_;
    	extent = boundingBoxExtent_;
    }

    int elementIndex() { return elementIndex_; }

private:

    /** \brief Calculates a box which bounds the curvilinear element.
     *
     * Naive method:
     * 1) Find bounding box which bounds interpolatory vertex set
     * 2) Enlarge the box to account for excess curvature
     *
     * Problem:
     * 1) The enlargement factor is completely empirical - need some reasonable estimates
     * 2) Wasteful to use the same enlargement factor in all directions
     *
     * TODO: Exact method:
     * 1) For each coordinate (x,y,z) select candidate faces for minimizing and maximizing.
     * 1.1) For example, 3 faces that share the vertex with largest X, 3 faces that share the lowest
     * 2) For each face solve all optimization problems, then maximize/minimize over solutions
     * 2.2) Sample problem: max { x(u,v) } over face
     *
     * Advantage:
     * 1) Assuming excess factor of 1.1, there will be approx 1.83 times less elements per octant,
     * hence that much faster per one octree request
     *
     * Problem:
     * 1) Optimization method tricky since not convex
     * 2) Increased precomputing time
     * 2.1) 18 2D non-convex optimization problems per element, vs 1.83 times more 3D single-extrema Newton Methods.
     *      It is not obvious that this effort is worthwhile
     *
     *
     * */

    void calculateBoundingBox()
    {

    	// Might make sense that it get smaller with order
    	const int EXCESS_CURVATURE_FACTOR = 1.1;

    	ElementGeometry thisGeometry = elementGeometry();

    	GlobalCoordinate maxBoxCorner = thisGeometry.vertex(0);
    	GlobalCoordinate minBoxCorner = maxBoxCorner;


    	// 1) Find corners of surrounding box for the interpolatory vertex set
    	for (int i = 1; i < thisGeometry.nVertex(); i++)  {
    		GlobalCoordinate tmpVertex = thisGeometry.vertex(i);

    		if (tmpVertex[0] > maxBoxCorner[0])  { maxBoxCorner[0] = tmpVertex[0]; }
    		if (tmpVertex[1] > maxBoxCorner[1])  { maxBoxCorner[1] = tmpVertex[1]; }
    		if (tmpVertex[2] > maxBoxCorner[2])  { maxBoxCorner[2] = tmpVertex[2]; }
    		if (tmpVertex[0] < minBoxCorner[0])  { minBoxCorner[0] = tmpVertex[0]; }
    		if (tmpVertex[1] < minBoxCorner[1])  { minBoxCorner[1] = tmpVertex[1]; }
    		if (tmpVertex[2] < minBoxCorner[2])  { minBoxCorner[2] = tmpVertex[2]; }
    	}

    	// Calculate center and extent from two corners
    	boundingBoxCenter_ = maxBoxCorner + minBoxCorner;
    	boundingBoxCenter_ *= 0.5;

    	boundingBoxExtent_ = maxBoxCorner - minBoxCorner;
    	boundingBoxExtent_ *= 0.5 * EXCESS_CURVATURE_FACTOR;
    }





private: // Private members

    const GridBaseType & gridbase_;

    int elementIndex_;

    GlobalCoordinate boundingBoxCenter_;
    GlobalCoordinate boundingBoxExtent_;

};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEAROCTREENODE_HH
