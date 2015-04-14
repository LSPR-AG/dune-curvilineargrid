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
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/curvilineargeometry/interpolation/curvilinearelementinterpolator.hh>
#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>
#include <dune/curvilineargeometry/curvilineargeometry.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>
#include <dune/curvilineargrid/common/vectorhelper.hh>

#include <dune/curvilineargrid/curvilineargrid/grid.hh>

#include <dune/curvilineargrid/utility/diagnosticshelper.hh>



namespace Dune {

template <class GridType>
class CurvilinearGridDiagnostic
{
private:

	typedef typename  GridType::ctype            ctype;
	typedef typename  GridType::LoggingMessage   LoggingMessage;

	typedef typename  GridType::LeafGridView LeafGridView;
	typedef typename  LeafGridView::IntersectionIterator IntersectionIterator;

	static const int  cdim     = GridType::dimension;
	static const bool isCached = GridType::is_cached;

    static const unsigned int LOG_CATEGORY_DEBUG  = LoggingMessage::Category::DEBUG;
    static const unsigned int LOG_CATEGORY_ERROR  = LoggingMessage::Category::ERROR;

    const int MASTER_PROCESS = 0;

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
    static const int BOUNDARY_SEGMENT_PARTITION_TYPE = GridStorageType::BOUNDARY_SEGMENT_PARTITION_TYPE;


	typedef typename LeafGridView::template Codim<ELEMENT_CODIM >::Entity   ElementType;
	typedef typename LeafGridView::template Codim<FACE_CODIM >::Entity      FaceType;

	typedef typename LeafGridView::template Codim<ELEMENT_CODIM>::Iterator  ElementIterator;

	typedef typename LeafGridView::template Codim<ELEMENT_CODIM>::Geometry  ElementGeometry;
	typedef typename LeafGridView::template Codim<FACE_CODIM>::Geometry     FaceGeometry;


    typedef std::vector<std::vector<double> >  MeshStatType;

public:

	CurvilinearGridDiagnostic(
		MPIHelper & mpihelper,
		GridType & grid) :
			mpihelper_(mpihelper),
			grid_(grid),
			loggingmessage_(LoggingMessage::getInstance())
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
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Started collecting mesh statistics");

		MeshStatType meshStatistics;
		meshStatistics[0].push_back(grid_.numInternal(ELEMENT_CODIM));
		meshStatistics[1].push_back(grid_.numBoundarySegments());
		meshStatistics[2].push_back(grid_.numProcessBoundaries());
		meshStatistics[12].push_back(0.0);  // processBoundarySurfaceArea
		meshStatistics[13].push_back(0.0);  // domainBoundarySurfaceArea

		LeafGridView leafView = grid_.leafGridView();

		// Iterate over entities of this codimension
		ElementIterator iElemB = leafView.template begin<ELEMENT_CODIM>();
		ElementIterator iElemE   = leafView.template end<ELEMENT_CODIM>();
		for (ElementIterator it = iElemB; it != iElemE; ++it)
		{
			const ElementType & e = *it;

			if (e.partitionType() != Dune::PartitionType::GhostEntity)
			{
				// Perform Volume tests
				ElementGeometry thisGeometry = e.geometry ();
				DiagnosticsHelper<GridType>::volumeTests(thisGeometry, meshStatistics);


				// Perform Face Tests
		        const IntersectionIterator nend = leafView.iend(e);
				for( IntersectionIterator nit = leafView.ibegin(e); nit != nend; ++nit )
				{
				  if (!nit->boundary())   {
					  // Domain Boundaries
					  FaceGeometry thisGeometry = nit->geometry();
					  DiagnosticsHelper<GridType>::domainBoundaryTests(thisGeometry, meshStatistics);
				  }
				  else if (nit->outside().partitionType() == Dune::PartitionType::GhostEntity)
				  {
					  // Process Boundaries
					  FaceGeometry thisGeometry = nit->geometry();
					  DiagnosticsHelper<GridType>::processBoundaryTests(thisGeometry, meshStatistics);
				  }
				}
			}
		}

		// Write diagnostics result to a file
		Dune::DiagnosticsHelper<GridType>::writeAnalyticTestResult(filename, meshStatistics, mpihelper_);
	}



	// VTK Tests
	// *************************************************************

	// Writes the mesh to VTK, including the additional constructions made by the mesh generator
	// [TODO] Implement LoggingMessage
	// [TODO] Use physicalTag=-1 for entities that do not have physicalTag
	// [TODO] Use another scalar = "OwnerRank" for entities
	void vtkWriteMesh (
		std::vector<bool> withElements,
		std::vector<bool> withFaces,
		std::vector<bool> withEdges,
		int nDiscretizationPoints,
		bool interpolate,
		bool explode,
		bool writeVtkEdges,
		bool writeVtkTriangles
	)
	{
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Started Writing Grid to VTK");

    	Dune::CurvilinearVTKWriter<GridType> vtkCurvWriter(mpihelper_);

    	std::vector<PartitionType>  structTypeSet
    	{
    		PartitionType::InteriorEntity,
    		PartitionType::GhostEntity,
    		PartitionType::BorderEntity
    	};

    	for (int iType = 0; iType < structTypeSet.size(); iType++)
    	{
    		if (withEdges[iType])     { addVTKentitySet<EDGE_CODIM>    (vtkCurvWriter, structTypeSet[iType], nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles); }
    		if (withFaces[iType])     { addVTKentitySet<FACE_CODIM>    (vtkCurvWriter, structTypeSet[iType], nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles); }

    		// For elements there is no Domain and Process boundaries, so only Internal and Ghost requests are processed
    		if ((iType < 2) && (withElements[iType]))
    		                          { addVTKentitySet<ELEMENT_CODIM> (vtkCurvWriter, structTypeSet[iType], nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles); }
    	}

    	// Write boundarySegments if requested
    	if (withFaces[3])  { addVTKentitySet<FACE_CODIM> (vtkCurvWriter, PartitionType::InteriorEntity, nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles, DOMAIN_BOUNDARY_TYPE); }


		// Writing Mesh
    	// *************************************************************************
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Writing VTK File");
    	//vtkCurvWriter.writeVTK("./curvreader_output_process_" + std::to_string(rank_) + ".vtk");
    	vtkCurvWriter.writeParallelVTU("./curvreader_output");
    	loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Finished writing");
	}

	// Writes OCTree to VTK
	// [TODO] Move this functionality here from OCTree impl, edit it
	void vtkWriteOctree() {}






protected:



	template<class Entity, int codim>
	void addVTKentity(
			Dune::CurvilinearVTKWriter<GridType> & vtkCurvWriter,
			const Entity & e,
			StructuralType typeTag,
			int nDiscretizationPoints,
			bool interpolate,
			bool explode,
			bool VTK_WRITE_EDGES,
			bool VTK_WRITE_TRIANGLES
	)
	{
		typedef typename GridType::GridBaseType::template Codim<codim>::EntityGeometry   EntityGeometry;
		typedef typename EntityGeometry::GlobalCoordinate  GlobalCoordinate;

		Dune::GeometryType gt              = e.type();

		// Constructing a geometry is quite expensive, do it only once
		//EntityGeometry geom = it->geometry();
		EntityGeometry geom = grid_.template entityBaseGeometry<codim>(e);
		std::vector<GlobalCoordinate>  interpVertices = geom.vertexSet();

		PhysicalTagType physicalTag        = grid_.template entityPhysicalTag<codim>(e);
		InterpolatoryOrderType interpOrder = grid_.template entityInterpolationOrder<codim>(e);
		std::vector<int>         tags  { physicalTag, typeTag, mpihelper_.rank() };

		vtkCurvWriter.template addCurvilinearElement<cdim - codim>(
	    			gt,
	    			interpVertices,
	    			tags,
	    			interpOrder,
	    			nDiscretizationPoints,
	    			interpolate,
	    			explode,
	    			VTK_WRITE_EDGES,
	    			VTK_WRITE_TRIANGLES);
	}



	template <int codim>
	void addVTKentitySet(
		Dune::CurvilinearVTKWriter<GridType> & vtkCurvWriter,
		PartitionType ptype,
		int nDiscretizationPoints,
		bool interpolate,
		bool explode,
		bool VTK_WRITE_EDGES,
		bool VTK_WRITE_TRIANGLES,
		StructuralType boundaryType = NO_BOUNDARY_TYPE
	)
	{


		LeafGridView leafView = grid_.leafGridView();

		// To insert boundary segments into writer
		if (boundaryType == DOMAIN_BOUNDARY_TYPE)
		{
			assert(codim == FACE_CODIM);

			for( auto && e : elements( leafView, Partitions::interior ) )
			{
		        const IntersectionIterator nend = leafView.iend(e);
				for( IntersectionIterator nit = leafView.ibegin(e); nit != nend; ++nit )
				{
				  // Checks if intersection is border
				  if(nit->boundary())
				  {
					// Gets indexInInside from intersection
					const FaceType & face = e.template subEntity<FACE_CODIM>(nit->indexInInside());
					addVTKentity<FaceType, FACE_CODIM>(vtkCurvWriter, face, BOUNDARY_SEGMENT_PARTITION_TYPE, nDiscretizationPoints, interpolate, explode, VTK_WRITE_EDGES, VTK_WRITE_TRIANGLES);
				  }
				}
			}
		} else // To insert entities
		{
			for( auto && e : entities( leafView, Dune::Dim<cdim - codim>()) )
			{
			  Dune::PartitionType thisPType      = e.partitionType();

			  // If we requested to output entities of this type, we will write them to VTK
			  if (thisPType == ptype)
			  {
				  typedef typename LeafGridView::template Codim<codim >::Entity   EntityType;
				  addVTKentity<EntityType, codim>(vtkCurvWriter, e, ptype, nDiscretizationPoints, interpolate, explode, VTK_WRITE_EDGES, VTK_WRITE_TRIANGLES);
			  }
			}
		}
	}






private:

	LoggingMessage & loggingmessage_;

	MPIHelper & mpihelper_;
	int rank_;
	int size_;

	GridType & grid_;


};


} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
