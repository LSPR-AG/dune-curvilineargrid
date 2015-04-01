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

#include <dune/curvilineargrid/curvilineargridbase/curvilineargridstorage.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilineargridbase.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearoctreenode.hh>
#include <dune/curvilineargrid/curvilineargridbase/curvilinearlooseoctree.hh>



namespace Dune {

template <class GridType>
class CurvilinearGridDiagnostic
{
private:

	typedef typename  GridType::ctype            ctype;
	typedef typename  GridType::GridStorageType  GridStorageType;
	typedef typename  GridType::GridBaseType     GridBaseType;
	typedef typename  GridType::LoggingMessage   LoggingMessage;

	static const int  cdim     = GridType::dimension;
	static const bool isCached = GridType::is_cached;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;
	typedef typename GridStorageType::InterpolatoryOrderType    InterpolatoryOrderType;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

	typedef typename GridStorageType::Vertex                Vertex;
	typedef typename GridStorageType::LocalIndexSet         LocalIndexSet;
	typedef typename GridStorageType::IndexSetIterator      IndexSetIterator;

	typedef typename GridBaseType::template Codim<FACE_CODIM>::EntityGeometry        GridFaceGeometry;
	typedef typename GridBaseType::template Codim<ELEMENT_CODIM>::EntityGeometry     GridElementGeometry;

    // Logging Message Typedefs
    static const unsigned int LOG_CATEGORY_DEBUG  = LoggingMessage::Category::DEBUG;
    static const unsigned int LOG_CATEGORY_ERROR  = LoggingMessage::Category::ERROR;

    static const unsigned int DomainBoundaryType   = GridStorageType::PartitionType::DomainBoundary;
    static const unsigned int ProcessBoundaryType  = GridStorageType::PartitionType::ProcessBoundary;
    static const unsigned int InternalType         = GridStorageType::PartitionType::Internal;
    static const unsigned int GhostType            = GridStorageType::PartitionType::Ghost;

	const int MASTER_PROCESS = 0;

public:

	CurvilinearGridDiagnostic(
		MPIHelper & mpihelper,
		GridBaseType & gridbase) :
			mpihelper_(mpihelper),
			gridbase_(gridbase),
			loggingmessage_(LoggingMessage::getInstance())
	{
		rank_ = mpihelper.rank();
		size_ = mpihelper.size();
	}


	// Runs analytic tests that collect statistics on the elements and the faces of the mesh
	// Statistics is communicated to MASTER_PROCESS, which writes it to a file
	// [TODO] When calculating total volume, calculate it separately for all tags and provide tags
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


    	typedef typename GridStorageType::PartitionType  CurvPT;
    	std::vector<StructuralType>  structTypeSet {InternalType, GhostType, DomainBoundaryType, ProcessBoundaryType};

    	for (int iType = 0; iType < structTypeSet.size(); iType++)
    	{
    		if (withEdges[iType])     { addVTKentity<EDGE_CODIM>    (vtkCurvWriter, structTypeSet[iType], nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles); }
    		if (withFaces[iType])     { addVTKentity<FACE_CODIM>    (vtkCurvWriter, structTypeSet[iType], nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles); }

    		// For elements there is no Domain and Process boundaries, so only Internal and Ghost requests are processed
    		if ((iType < 2) && (withElements[iType]))
    		                          { addVTKentity<ELEMENT_CODIM> (vtkCurvWriter, structTypeSet[iType], nDiscretizationPoints, interpolate, explode, writeVtkEdges, writeVtkTriangles); }
    	}


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

	// Primal Analytic Tests
	// *************************************************************

	void analyticTests(std::vector<std::vector<double> > & rez)
	{
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Started collecting mesh statistics");


		rez[0].push_back(gridbase_.template nEntity(ELEMENT_CODIM, InternalType));
		rez[1].push_back(gridbase_.template nEntity(FACE_CODIM, DomainBoundaryType));
		rez[2].push_back(gridbase_.template nEntity(FACE_CODIM, ProcessBoundaryType));

		// 1) Collect statistics related to the elements of the mesh
		// ***********************************************************************8
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Collecting element statistics");
		IndexSetIterator elemIterB = gridbase_.template entityIndexBegin(ELEMENT_CODIM);
		IndexSetIterator elemIterE = gridbase_.template entityIndexEnd(ELEMENT_CODIM);



		for (IndexSetIterator elemIter = elemIterB;  elemIter != elemIterE;  elemIter++)
		{
			GridElementGeometry thisGeometry = gridbase_.template entityGeometry<ELEMENT_CODIM>(*elemIter);
			std::vector<Vertex> cr       = thisGeometry.cornerSet();

			Vertex CoM = cr[0] + cr[1] + cr[2] + cr[3];
			CoM /= 4;

			std::vector<double>     comCornerRadius { (CoM-cr[0]).two_norm(), (CoM-cr[1]).two_norm(), (CoM-cr[2]).two_norm(), (CoM-cr[3]).two_norm() };
			std::vector<Vertex> linearEdges   { cr[3]-cr[0], cr[3]-cr[1], cr[3]-cr[2], cr[2]-cr[0], cr[2]-cr[1], cr[1]-cr[0] };
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
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Collecting Process Boundary statistics");
		rez[12].push_back(0.0);  // processBoundarySurfaceArea
		IndexSetIterator pbIterB = gridbase_.template entityIndexBegin(FACE_CODIM, ProcessBoundaryType);
		IndexSetIterator pbIterE = gridbase_.template entityIndexEnd(FACE_CODIM, ProcessBoundaryType);

		for (IndexSetIterator pbIter = pbIterB;  pbIter != pbIterE;  pbIter++)
		{
			GridFaceGeometry faceGeom = gridbase_.template entityGeometry<FACE_CODIM>(*pbIter);
			double faceCurvilinearArea = faceGeom.volume(1.0e-5);
			rez[12][0] += faceCurvilinearArea;

			std::cout << "Area: " << faceCurvilinearArea << std::endl;
		}


		// 3) Collect statistics related to the domain boundary of the mesh
		// ***********************************************************************
		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Collecting Domain Boundary Statistics");

		rez[13].push_back(0.0);  // domainBoundarySurfaceArea
		IndexSetIterator dbIterB = gridbase_.template entityIndexBegin(FACE_CODIM, DomainBoundaryType);
		IndexSetIterator dbIterE = gridbase_.template entityIndexEnd(FACE_CODIM, DomainBoundaryType);

		for (IndexSetIterator dbIter = dbIterB;  dbIter != dbIterE;  dbIter++)
		{
			GridFaceGeometry faceGeom = gridbase_.template entityGeometry<FACE_CODIM>(*dbIter);
			double faceCurvilinearArea = faceGeom.volume(1.0e-5);
			rez[13][0] += faceCurvilinearArea;

			std::cout << "size=" << faceGeom.vertexSet().size() <<" vertices = " << Dune::VectorHelper::vector2string(faceGeom.vertexSet()) << std::endl;

			Dune::CurvilinearGeometry<double,2,3>::PolynomialVector polyMap = faceGeom.interpolatoryVectorAnalytical();

			std::cout << "Polynomial[0] = " << polyMap[0].to_string() << std::endl;
			std::cout << "Polynomial[1] = " << polyMap[0].to_string() << std::endl;
			std::cout << "Polynomial[2] = " << polyMap[0].to_string() << std::endl;

			std::cout << "Area: " << faceCurvilinearArea << std::endl;
		}
	}


	template <int codim>
	void addVTKentity(
		Dune::CurvilinearVTKWriter<GridType> & vtkCurvWriter,
		StructuralType structtype,
		int nDiscretizationPoints,
		bool interpolate,
		bool explode,
		bool VTK_WRITE_EDGES,
		bool VTK_WRITE_TRIANGLES
	)
	{
		typedef typename GridBaseType::template Codim<codim>::EntityGeometry     EntityGeometry;

		loggingmessage_.template write<LOG_CATEGORY_DEBUG>( __FILE__, __LINE__, "CurvilinearDiagnostics: Started writing entities codim=" + std::to_string(codim) + " type=" + gridbase_.PartitonTypeName(structtype));

		IndexSetIterator elemIterB =  gridbase_.template entityIndexBegin(codim, structtype);
		IndexSetIterator elemIterE =  gridbase_.template entityIndexEnd(codim, structtype);

		for (IndexSetIterator elemIter = elemIterB;  elemIter != elemIterE;  elemIter++)
		{
			LocalIndexType           thisLocalIndex = *elemIter;
			StructuralType           thisStructType = gridbase_.entityStructuralType(codim, thisLocalIndex);

			assert(thisStructType == structtype);  // Checking grid self-consistency

			PhysicalTagType          physicalTag       = gridbase_.physicalTag(codim, thisLocalIndex);
			EntityGeometry           thisGeometry      = gridbase_.template entityGeometry<codim>(thisLocalIndex);
			Dune::GeometryType       gt                = thisGeometry.type();
			InterpolatoryOrderType   order             = thisGeometry.order();
			std::vector<Vertex>      point             = thisGeometry.vertexSet();
			std::vector<int>         tags  { physicalTag, structtype, rank_ };

	    	vtkCurvWriter.template addCurvilinearElement<cdim - codim>(
	    			gt,
	    			point,
	    			tags,
	    			order,
	    			nDiscretizationPoints,
	    			interpolate,
	    			explode,
	    			VTK_WRITE_EDGES,
	    			VTK_WRITE_TRIANGLES);
		}
	}

	// Auxiliary Methods
	// *************************************************************


	// Initializes a FieldVector in 1 line
	Vertex initVector(ctype a, ctype b, ctype c)
	{
		Vertex rez;
		rez[0] = a;
		rez[1] = b;
		rez[2] = c;
		return rez;
	}

	// Dot product between FieldVectors
	// [TODO] Replace with existing Dune functionality if found
	ctype GridVectorDot(Vertex a, Vertex b)  { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

	// Cross product between FieldVectors
	// [TODO] Replace with existing Dune functionality if found
	Vertex GridVectorTimes(Vertex a, Vertex b)  { return initVector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]); }

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

	LoggingMessage & loggingmessage_;

	MPIHelper & mpihelper_;
	int rank_;
	int size_;

	GridBaseType & gridbase_;


};


} // namespace Dune

#endif  // DUNE_CURVILINEARGRIDDIAGNOSTIC_HH
