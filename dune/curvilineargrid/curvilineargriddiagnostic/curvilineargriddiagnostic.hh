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

#include <dune/curvilineargrid/feedback/loggingmessage.hh>

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>



namespace Dune {

template <class ct>
class CurvilinearGridDiagnostic
{
private:
	typedef typename Dune::CurvilinearGridBase<ct>::Vertex                GridVertex;
	typedef typename Dune::CurvilinearGridBase<ct>::Index2IndexMap        GridIndexMap;
	typedef typename Dune::CurvilinearGridBase<ct>::IndexMapIterator      GridIndexMapIterator;

	typedef typename Dune::CurvilinearGridBase<ct>::FaceGeometry          GridFaceGeometry;
	typedef typename Dune::CurvilinearGridBase<ct>::ElementGeometry       GridElementGeometry;

    // Logging Message Typedefs
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;
    static const unsigned int LOG_CATEGORY_ERROR = Dune::LoggingMessage::Category::ERROR;

	const int VTK_INTERNAL          = Dune::VtkEntityStructuralType::Internal;
	const int VTK_GHOST             = Dune::VtkEntityStructuralType::Ghost;
	const int VTK_DOMAIN_BOUNDARY   = Dune::VtkEntityStructuralType::DomainBoundary;
	const int VTK_PROCESS_BOUNDARY  = Dune::VtkEntityStructuralType::ProcessBoundary;

	const int MASTER_PROCESS = 0;

public:

	CurvilinearGridDiagnostic(
		bool verbose,
		bool processVerbose,
		MPIHelper & mpihelper,
		Dune::CurvilinearGridBase<ct> & gridbase) :
			verbose_(verbose),
			processVerbose_(processVerbose),
			mpihelper_(mpihelper),
			gridbase_(gridbase)
	{
		rank_ = mpihelper.rank();
		size_ = mpihelper.size();
	}


	// Runs analytic tests that collect statistics on the elements and the faces of the mesh
	// Statistics is communicated to MASTER_PROCESS, which writes it to a file
	void runAnalyticTest(std::string filename)
	{
		 std::ofstream diagnosticFile;

		 if (rank_ == MASTER_PROCESS)  { diagnosticFile.open(filename.c_str()); }

		 std::vector<std::vector<double> > meshStatistics (14, std::vector<double>());
		 analyticTests(meshStatistics);

		 writeCommunicateVector(diagnosticFile, "NELEMENT",                       meshStatistics[0]);
		 writeCommunicateVector(diagnosticFile, "NDBFACE",                        meshStatistics[1]);
		 writeCommunicateVector(diagnosticFile, "NPBFACE",                        meshStatistics[2]);
		 writeCommunicateVector(diagnosticFile, "ElementShortestEdge",            meshStatistics[3]);
		 writeCommunicateVector(diagnosticFile, "ElementLongestEdge",             meshStatistics[4]);
		 writeCommunicateVector(diagnosticFile, "ElementLinearVolume",            meshStatistics[5]);
		 writeCommunicateVector(diagnosticFile, "ElementCurvilinearVolume",       meshStatistics[6]);
		 writeCommunicateVector(diagnosticFile, "ElementEnclosedSphereRadius",    meshStatistics[7]);
		 writeCommunicateVector(diagnosticFile, "ElementSurroundingSphereRadius", meshStatistics[8]);
		 writeCommunicateVector(diagnosticFile, "ElementLinearQuality1",          meshStatistics[9]);
		 writeCommunicateVector(diagnosticFile, "ElementLinearQuality2",          meshStatistics[10]);
		 writeCommunicateVector(diagnosticFile, "ElementCurvilinearQuality1",     meshStatistics[11]);
		 writeCommunicateSum(diagnosticFile, "ProcssBoundarySurfaceArea",         meshStatistics[12][0]);
		 writeCommunicateSum(diagnosticFile, "DomainBoundarySurfaceArea",         meshStatistics[13][0]);

		 if (rank_ == MASTER_PROCESS)  { diagnosticFile.close(); }
	}



	// VTK Tests
	// *************************************************************

	// Writes the mesh to VTK, including the additional constructions made by the mesh generator
	// [TODO] Implement LoggingMessage
	// [TODO] Use physicalTag=-1 for entities that do not have physicalTag
	// [TODO] Use another scalar = "OwnerRank" for entities
	void vtkWriteMesh (
		bool withElements,
		bool withGhostElements,
		bool withDomainBoundaries,
		bool withProcessBoundaries,
		int nDiscretizationPoints,
		bool interpolate,
		bool explode
	)
	{
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Started Writing Grid to VTK");

    	Dune::CurvilinearVTKWriter<3> vtkCurvWriter(verbose_, processVerbose_, mpihelper_);
    	bool VTK_WRITE_EDGES = false;         // Whether to write volume discretization edges to file
    	bool VTK_WRITE_TRIANGLES = true;      // Whether to write surface discretization triangles to file


    	// Writing Elements
    	// *************************************************************************
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Started Writing Elements");
		if (withElements)
		{
			GridIndexMapIterator elemIterB =  gridbase_.elementIndexBegin();
			GridIndexMapIterator elemIterE =  gridbase_.elementIndexEnd();

			for (GridIndexMapIterator elemIter = elemIterB;  elemIter != elemIterE;  elemIter++)
			{
				int                      elementLocalIndex = (*elemIter).second;
				int                      physicalTag       = gridbase_.elementPhysicalTag(elementLocalIndex);
				GridElementGeometry      thisGeometry      = gridbase_.elementGeometry(elementLocalIndex);
				Dune::GeometryType       gt                = thisGeometry.type();
				int                      order             = thisGeometry.order();
				std::vector<GridVertex>  point             = thisGeometry.vertexSet();

				Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "writing element");

		    	vtkCurvWriter.addCurvilinearElement(
		    			gt,
		    			point,
		    			order,
		    			physicalTag,
		    			nDiscretizationPoints,
		    			interpolate,
		    			explode,
		    			VTK_WRITE_EDGES,
		    			VTK_WRITE_TRIANGLES,
		    			VTK_INTERNAL);
			}
		}


    	// Writing Ghost Elements
    	// *************************************************************************
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Started Writing Ghost Elements");
		if (withGhostElements)
		{
			GridIndexMapIterator ghostIterB =  gridbase_.ghostIndexBegin();
			GridIndexMapIterator ghostIterE =  gridbase_.ghostIndexEnd();

			for (GridIndexMapIterator ghostIter = ghostIterB;  ghostIter != ghostIterE;  ghostIter++)
			{
				int                      ghostLocalIndex   = (*ghostIter).second;
				int                      physicalTag       = gridbase_.ghostElementPhysicalTag(ghostLocalIndex);
				GridElementGeometry      thisGeometry      = gridbase_.ghostGeometry(ghostLocalIndex);
				Dune::GeometryType       gt                = thisGeometry.type();
				int                      order             = thisGeometry.order();
				std::vector<GridVertex>  point             = thisGeometry.vertexSet();

		    	vtkCurvWriter.addCurvilinearElement(
		    			gt,
		    			point,
		    			order,
		    			physicalTag,
		    			nDiscretizationPoints,
		    			interpolate,
		    			explode,
		    			VTK_WRITE_EDGES,
		    			VTK_WRITE_TRIANGLES,
		    			VTK_GHOST);
			}
		}


    	// Writing Domain Boundary Faces
    	// *************************************************************************
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Started Writing Domain Boundaries");
		if (withDomainBoundaries)
		{
			// Add DomainBoundaries to VTK
			GridIndexMapIterator faceIterB =  gridbase_.faceDomainBoundaryIndexBegin();
			GridIndexMapIterator faceIterE =  gridbase_.faceDomainBoundaryIndexEnd();

			for (GridIndexMapIterator faceIter = faceIterB;  faceIter != faceIterE;  faceIter++)
			{
				int                      faceLocalIndex    = (*faceIter).second;
				int                      physicalTag       = gridbase_.facePhysicalTag(faceLocalIndex);
				GridFaceGeometry         thisGeometry      = gridbase_.faceGeometry(faceLocalIndex);
				Dune::GeometryType       gt                = thisGeometry.type();
				int                      order             = thisGeometry.order();
				std::vector<GridVertex>  point             = thisGeometry.vertexSet();

		    	vtkCurvWriter.addCurvilinearElement(
		    			gt,
		    			point,
		    			order,
		    			physicalTag,
		    			nDiscretizationPoints,
		    			interpolate,
		    			explode,
		    			VTK_WRITE_EDGES,
		    			VTK_WRITE_TRIANGLES,
		    			VTK_DOMAIN_BOUNDARY);
			}
		}


		// Writing Process Boundary Faces
    	// *************************************************************************
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Started Writing Process Boundaries");
		if (withProcessBoundaries)
		{
			// Add DomainBoundaries to VTK
			GridIndexMapIterator faceIterB =  gridbase_.faceProcessBoundaryIndexBegin();
			GridIndexMapIterator faceIterE =  gridbase_.faceProcessBoundaryIndexEnd();

			for (GridIndexMapIterator faceIter = faceIterB;  faceIter != faceIterE;  faceIter++)
			{
				int                      faceLocalIndex    = (*faceIter).second;
				int                      physicalTag       = gridbase_.facePhysicalTag(faceLocalIndex);
				GridFaceGeometry         thisGeometry      = gridbase_.faceGeometry(faceLocalIndex);
				Dune::GeometryType       gt                = thisGeometry.type();
				int                      order             = thisGeometry.order();
				std::vector<GridVertex>  point             = thisGeometry.vertexSet();

		    	vtkCurvWriter.addCurvilinearElement(
		    			gt,
		    			point,
		    			order,
		    			physicalTag,
		    			nDiscretizationPoints,
		    			interpolate,
		    			explode,
		    			VTK_WRITE_EDGES,
		    			VTK_WRITE_TRIANGLES,
		    			VTK_PROCESS_BOUNDARY);
			}
		}


		// Writing Mesh
    	// *************************************************************************
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Writing VTK File");
    	vtkCurvWriter.writeVTK("./curvreader_output_process_" + std::to_string(rank_) + ".vtk");
    	vtkCurvWriter.writeParallelVTU("./curvreader_output");
    	Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Finished writing");
	}

	// Writes OCTree to VTK
	// [TODO] Move this functionality here from OCTree impl, edit it
	void vtkWriteOctree() {}






protected:

	// Primal Analytic Tests
	// *************************************************************



	void analyticTests(std::vector<std::vector<double> > & rez)
	{
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Started collecting mesh statistics");


		rez[0].push_back(gridbase_.nElement());
		rez[1].push_back(gridbase_.nFaceDomainBoundary());
		rez[2].push_back(gridbase_.nFaceProcessBoundary());

		// 1) Collect statistics related to the elements of the mesh
		// ***********************************************************************8
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Collecting element statistics");
		GridIndexMapIterator elemIterB = gridbase_.elementIndexBegin();
		GridIndexMapIterator elemIterE = gridbase_.elementIndexEnd();



		for (GridIndexMapIterator elemIter = elemIterB;  elemIter != elemIterE;  elemIter++)
		{
			GridElementGeometry thisGeometry = gridbase_.elementGeometry((*elemIter).second);
			std::vector<GridVertex> cr       = thisGeometry.cornerSet();

			GridVertex CoM = cr[0] + cr[1] + cr[2] + cr[3];
			CoM /= 4;

			std::vector<double>     comCornerRadius { (CoM-cr[0]).two_norm(), (CoM-cr[1]).two_norm(), (CoM-cr[2]).two_norm(), (CoM-cr[3]).two_norm() };
			std::vector<GridVertex> linearEdges   { cr[3]-cr[0], cr[3]-cr[1], cr[3]-cr[2], cr[2]-cr[0], cr[2]-cr[1], cr[1]-cr[0] };
			std::vector<double>     linearEdgeLen;
			for (int i = 0; i < linearEdges.size(); i++)  { linearEdgeLen.push_back(linearEdges[i].two_norm()); }

			double elementShortestEdge = vectorMin(linearEdgeLen);
			double elementLongestEdge = vectorMax(linearEdgeLen);

			double surfaceAreaLinear = (
				GridVectorTimes(linearEdges[0], linearEdges[1]).two_norm() +
				GridVectorTimes(linearEdges[0], linearEdges[2]).two_norm() +
				GridVectorTimes(linearEdges[1], linearEdges[2]).two_norm() +
				GridVectorTimes(linearEdges[3], linearEdges[4]).two_norm() ) / 2.0;

			double volumeLinear             = fabs( GridVectorDot(linearEdges[0], GridVectorTimes(linearEdges[1], linearEdges[2])) / 6.0 );
			double volumeCurvilinear        = thisGeometry.volume(1.0e-5);
			double radiusSphereEnclosed     = 6 * volumeLinear / surfaceAreaLinear;
			double radiusSphereSurround     = vectorMax(comCornerRadius);
			double qualityFactorLinear1     = elementShortestEdge / elementLongestEdge;
			double qualityFactorLinear2     = radiusSphereEnclosed / radiusSphereSurround;
			double qualityFactorCurvilinear = volumeLinear / volumeCurvilinear;

			rez[3].push_back(elementShortestEdge);
			rez[4].push_back(elementLongestEdge);
			rez[5].push_back(volumeLinear);
			rez[6].push_back(volumeCurvilinear);
			rez[7].push_back(radiusSphereEnclosed);
			rez[8].push_back(radiusSphereSurround);
			rez[9].push_back(qualityFactorLinear1);
			rez[10].push_back(qualityFactorLinear2);
			rez[11].push_back(qualityFactorCurvilinear);
		}


		// 2) Collect statistics related to the process boundary of the mesh
		// ***********************************************************************
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Collecting Process Boundary statistics");
		rez[12].push_back(0.0);  // processBoundarySurfaceArea
		GridIndexMapIterator pbIterB = gridbase_.faceProcessBoundaryIndexBegin();
		GridIndexMapIterator pbIterE = gridbase_.faceProcessBoundaryIndexEnd();

		for (GridIndexMapIterator pbIter = pbIterB;  pbIter != pbIterE;  pbIter++)
		{
			GridFaceGeometry faceGeom = gridbase_.faceGeometry((*pbIter).second);
			double faceCurvilinearArea = faceGeom.volume(1.0e-5);
			rez[12][0] += faceCurvilinearArea;

			std::cout << "Area: " << faceCurvilinearArea << std::endl;
		}


		// 3) Collect statistics related to the domain boundary of the mesh
		// ***********************************************************************
		Dune::LoggingMessage::write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(mpihelper_, verbose_, processVerbose_, __FILE__, __LINE__, "CurvilinearDiagnostics: Collecting Domain Boundary Statistics");

		rez[13].push_back(0.0);  // domainBoundarySurfaceArea
		GridIndexMapIterator dbIterB = gridbase_.faceDomainBoundaryIndexBegin();
		GridIndexMapIterator dbIterE = gridbase_.faceDomainBoundaryIndexEnd();

		for (GridIndexMapIterator dbIter = dbIterB;  dbIter != dbIterE;  dbIter++)
		{
			GridFaceGeometry faceGeom = gridbase_.faceGeometry((*dbIter).second);
			double faceCurvilinearArea = faceGeom.volume(1.0e-5);
			rez[13][0] += faceCurvilinearArea;

			std::cout << "size=" << faceGeom.vertexSet().size() <<" vertices = " << vector2string(faceGeom.vertexSet()) << std::endl;

			Dune::CurvilinearGeometry<double,2,3>::PolynomialVector polyMap = faceGeom.interpolatoryVectorAnalytical();

			std::cout << "Polynomial[0] = " << polyMap[0].to_string() << std::endl;
			std::cout << "Polynomial[1] = " << polyMap[0].to_string() << std::endl;
			std::cout << "Polynomial[2] = " << polyMap[0].to_string() << std::endl;

			std::cout << "Area: " << faceCurvilinearArea << std::endl;
		}
	}


	// Auxiliary Methods
	// *************************************************************

    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
    template <class T>
    std::string vector2string(const T & V)
    {
        std::stringstream tmp_stream;

        int nEntry = V.size();
        if (nEntry == 0)  { tmp_stream << "Null"; }
        for (int i = 0; i < nEntry; i++) {
        	tmp_stream << V[i];
        	if (i != nEntry - 1) { tmp_stream << " "; }
        }
        return tmp_stream.str();
    }

	// Initializes a FieldVector in 1 line
	GridVertex initVector(ct a, ct b, ct c)
	{
		GridVertex rez;
		rez[0] = a;
		rez[1] = b;
		rez[2] = c;
		return rez;
	}

	// Dot product between FieldVectors
	// [TODO] Replace with existing Dune functionality if found
	ct GridVectorDot(GridVertex a, GridVertex b)  { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

	// Cross product between FieldVectors
	// [TODO] Replace with existing Dune functionality if found
	GridVertex GridVectorTimes(GridVertex a, GridVertex b)  { return initVector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]); }

	// Returns smallest entry of a vector of comparable objects
	template <class T>
	T vectorMin (std::vector<T> data)
	{
		T rez = data[0];
		for (int i = 0; i < data.size(); i++)  { rez = std::min(rez, data[i]); }
		return rez;
	}

	// Returns largest entry of a vector of comparable objects
	template <class T>
	T vectorMax (std::vector<T> data)
	{
		T rez = data[0];
		for (int i = 0; i < data.size(); i++)  { rez = std::max(rez, data[i]); }
		return rez;
	}


	// Collects data from all processes and writes it to a file on the master process
	template<class T>
	void writeCommunicateVector(std::ofstream & filestr, std::string title, std::vector<T> & data)
	{
		Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

		// 1) Communicate size of communication to root
		// *************************************************************
		int dataSize = data.size();
		std::vector<int> dataSizeArray(size_);
		collective_comm.gather (&dataSize, reinterpret_cast<int*> (dataSizeArray.data()), 1, MASTER_PROCESS);


		// 2) Communicate actual data to root
		// *************************************************************
		int dataSizeTotal = 0;
		std::vector<int> displ;
		for (int i = 0; i < size_; i++)
		{
			dataSizeTotal += dataSizeArray[i];
			displ.push_back((i == 0) ? 0  :  displ[i-1] + dataSizeArray[i-1]);
		}
		std::vector<T> dataTotal(dataSizeTotal);
		collective_comm.gatherv (data.data(), dataSize, reinterpret_cast<T*> (dataTotal.data()), dataSizeArray.data(), displ.data(), MASTER_PROCESS);


		// 3) Write output to file on Master process
		// *************************************************************
		if (rank_ == MASTER_PROCESS)
		{
			filestr << "<" << title << ">";
			for (int i = 0; i < dataTotal.size(); i++) { filestr << dataTotal[i] << " ";  }
			filestr << "</" << title << ">" << std::endl;
		}
	}


	// Sums data over all processes and writes result
	template <class T>
	void writeCommunicateSum(std::ofstream & filestr, std::string title, T data)
	{
		Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();
		T sum = collective_comm.sum(data);
		if (rank_ == MASTER_PROCESS)  { filestr << "<" << title << ">" << sum << "</" << title << ">" << std::endl; }
	}




private:

	bool verbose_;
	bool processVerbose_;

	MPIHelper & mpihelper_;
	int rank_;
	int size_;

	Dune::CurvilinearGridBase<ct> & gridbase_;


};


} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
