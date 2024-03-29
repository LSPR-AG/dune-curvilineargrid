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

#ifndef DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
#define DUNE_CURVILINEARGRIDDIAGNOSTIC_HH

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
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargrid/grid.hh>

#include <dune/curvilineargrid/utility/diagnosticshelper.hh>



namespace Dune {

namespace CurvGrid {


template <class GridType>
class CurvilinearGridDiagnostic
{
private:

	typedef typename  GridType::ctype            ctype;

	typedef typename  GridType::LeafGridView LeafGridView;
	typedef typename  LeafGridView::IntersectionIterator IntersectionIterator;

	static const int  cdim     = GridType::dimension;
	static const bool isCached = GridType::is_cached;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridType::GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridType::GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridType::GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridType::GridStorageType::ELEMENT_CODIM;

    typedef typename GridType::GridStorageType::StructuralType          StructuralType;
    typedef typename GridType::GridStorageType::PhysicalTagType         PhysicalTagType;
    typedef typename GridType::GridStorageType::InterpolatoryOrderType  InterpolatoryOrderType;

    typedef typename GridType::GridStorageType  GridStorageType;
    static const int NO_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::None;
    static const int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;

	typedef typename LeafGridView::template Codim<ELEMENT_CODIM >::Entity   ElementType;
	typedef typename LeafGridView::template Codim<FACE_CODIM >::Entity      FaceType;

	typedef typename LeafGridView::template Codim<ELEMENT_CODIM>::Iterator  ElementIterator;

	typedef typename GridType::template Codim<ELEMENT_CODIM>::EntityGeometryMappingImpl  ElementGeometryImpl;
	typedef typename GridType::template Codim<ELEMENT_CODIM>::Geometry                   ElementGeometry;
	typedef typename GridType::template Codim<FACE_CODIM>::Geometry                      FaceGeometry;


    typedef std::vector<std::vector<double> >  MeshStatType;

public:

	CurvilinearGridDiagnostic(
		MPIHelper & mpihelper,
		GridType & grid) :
			mpihelper_(mpihelper),
			grid_(grid)
	{
		rank_ = mpihelper.rank();
		size_ = mpihelper.size();
	}


	// Runs analytic tests that collect statistics on the elements and the faces of the mesh
	// Statistics is communicated to MASTER_PROCESS, which writes it to a file
	// [TODO] When calculating total volume, calculate it separately for all tags and provide tags
	// [TODO] add test to check number of ghost elements
	void runAnalyticTest(std::string filename)
	{
		MeshStatType meshStatistics(14, std::vector<double>());
		meshStatistics[0].push_back(grid_.numInternal(ELEMENT_CODIM));
		meshStatistics[1].push_back(grid_.numBoundarySegments());
		meshStatistics[2].push_back(grid_.numProcessBoundaries());
		meshStatistics[12].push_back(0.0);  // processBoundarySurfaceArea
		meshStatistics[13].push_back(0.0);  // domainBoundarySurfaceArea

		LeafGridView leafView = grid_.leafGridView();
		int elemIterCount = 0;         // Count elements
		for (auto&& e : elements(leafView, Dune::Partitions::interior))
		{
	        LoggingMessage::writePatience("CurvilinearDiagnostics: Collecting mesh statistics", elemIterCount++, grid_.size(ELEMENT_CODIM));

	        //std::cout << "element index " << grid_.leafIndexSet().index(e) << std::endl;

			ElementGeometry thisGeometry = e.geometry();
			ElementGeometryImpl thisGeometryImpl = grid_.template entityBaseGeometry<ELEMENT_CODIM>(e);

			DiagnosticsHelper<GridType>::consistencyTests(thisGeometryImpl);                  // Test if element is well-formed
			DiagnosticsHelper<GridType>::elementStatistics(thisGeometry, meshStatistics);     // Collect statistics on the element

	        const IntersectionIterator nend = leafView.iend(e);
			for( IntersectionIterator nit = leafView.ibegin(e); nit != nend; ++nit )
			{
			  if ((nit->boundary()) && (!nit->neighbor()) )   {
				  // Domain Boundaries
				  FaceGeometry thisGeometry = nit->geometry();
				  DiagnosticsHelper<GridType>::domainBoundaryStatistics(thisGeometry, meshStatistics);
			  }
			  else if ((nit->neighbor())&&(nit->outside().partitionType() == Dune::PartitionType::GhostEntity))
			  {
				  // Process Boundaries
				  FaceGeometry thisGeometry = nit->geometry();
				  DiagnosticsHelper<GridType>::processBoundaryStatistics(thisGeometry, meshStatistics);
			  }
			}
		}

		// Write diagnostics result to a file
		DiagnosticsHelper<GridType>::writeAnalyticTestResult(filename, meshStatistics, mpihelper_);
	}


	// Writes OCTree to VTK
	// [TODO] Move this functionality here from OCTree impl, edit it
	void vtkWriteOctree() {}



private:

	MPIHelper & mpihelper_;
	int rank_;
	int size_;

	const GridType & grid_;


};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
