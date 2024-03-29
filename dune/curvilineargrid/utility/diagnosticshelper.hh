#ifndef DUNE_CURVILINEARDIAGNOSTICHELPER_HH
#define DUNE_CURVILINEARDIAGNOSTICHELPER_HH

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargeometry/utility/curvilinearvaliditycheck.hh>


namespace Dune
{

namespace CurvGrid {

//using namespace CurvGrid;

template <class GridType>
struct DiagnosticsHelper
{
	typedef typename  GridType::ctype            ctype;
	typedef typename  GridType::GridStorageType  GridStorageType;
	typedef typename  GridType::GridBaseType     GridBaseType;

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

	typedef typename GridStorageType::GlobalCoordinate                GlobalCoordinate;



	typedef std::vector<std::vector<double> >  MeshStatType;



	// Checks if the entities of the grid are well-defined
	template <class CurvGeometry>
	static void consistencyTests(CurvGeometry & elemGeom)
	{
		// Validate that the element is not self-intersecting
		Dune::CurvilinearValidityCheck::SelfIntersection(elemGeom);
	}


	// Collect statistics about the volume element
	template <class Geometry>
	static void elementStatistics(Geometry & elemGeom, MeshStatType & meshStatistics)
	{
		std::vector<GlobalCoordinate> cr;
		for (int i = 0; i < elemGeom.corners(); i++)  { cr.push_back(elemGeom.corner(i)); }

		GlobalCoordinate CoM = cr[0] + cr[1] + cr[2] + cr[3];
		CoM /= 4;

		std::vector<double> comCornerRadius { (CoM-cr[0]).two_norm(), (CoM-cr[1]).two_norm(), (CoM-cr[2]).two_norm(), (CoM-cr[3]).two_norm() };
		std::vector<GlobalCoordinate> linearEdges     { cr[3]-cr[0], cr[3]-cr[1], cr[3]-cr[2], cr[2]-cr[0], cr[2]-cr[1], cr[1]-cr[0] };
		std::vector<double> linearEdgeLen;
		for (int i = 0; i < linearEdges.size(); i++)  { linearEdgeLen.push_back(linearEdges[i].two_norm()); }

		double elementShortestEdge = vectorMin(linearEdgeLen);
		double elementLongestEdge = vectorMax(linearEdgeLen);

		double surfaceAreaLinear = (
			GridVectorTimes(linearEdges[0], linearEdges[1]).two_norm() +
			GridVectorTimes(linearEdges[0], linearEdges[2]).two_norm() +
			GridVectorTimes(linearEdges[1], linearEdges[2]).two_norm() +
			GridVectorTimes(linearEdges[3], linearEdges[4]).two_norm() ) / 2.0;

		double volumeLinear             = fabs( GridVectorDot(linearEdges[0], GridVectorTimes(linearEdges[1], linearEdges[2])) / 6.0 );
		double volumeCurvilinear        = elemGeom.volume();
		double radiusSphereEnclosed     = 6 * volumeLinear / surfaceAreaLinear;
		double radiusSphereSurround     = vectorMax(comCornerRadius);
		double qualityFactorLinear1     = elementShortestEdge / elementLongestEdge;
		double qualityFactorLinear2     = radiusSphereEnclosed / radiusSphereSurround;
		double qualityFactorCurvilinear = volumeLinear / volumeCurvilinear;

		meshStatistics[3].push_back(elementShortestEdge);
		meshStatistics[4].push_back(elementLongestEdge);
		meshStatistics[5].push_back(volumeLinear);
		meshStatistics[6].push_back(volumeCurvilinear);
		meshStatistics[7].push_back(radiusSphereEnclosed);
		meshStatistics[8].push_back(radiusSphereSurround);
		meshStatistics[9].push_back(qualityFactorLinear1);
		meshStatistics[10].push_back(qualityFactorLinear2);
		meshStatistics[11].push_back(qualityFactorCurvilinear);
	}


	// Collect statistics about a process boundary
	template <class Geometry>
	static void processBoundaryStatistics(Geometry & faceGeom, MeshStatType & meshStatistics)
	{
		double faceCurvilinearArea = faceGeom.volume();
		meshStatistics[12][0] += faceCurvilinearArea;
		LoggingMessage::template write<LOG_MSG_DVERB> (__FILE__, __LINE__, "Area: " + std::to_string(faceCurvilinearArea) );
	}


	// Collect statistics about a domain boundary
	template <class Geometry>
	static void domainBoundaryStatistics(Geometry & faceGeom, MeshStatType & meshStatistics)
	{
		double faceCurvilinearArea = faceGeom.volume();
		meshStatistics[13][0] += faceCurvilinearArea;
		LoggingMessage::template write<LOG_MSG_DVERB> (__FILE__, __LINE__, "Area: " + std::to_string(faceCurvilinearArea) );
	}








	// Runs analytic tests that collect statistics on the elements and the faces of the mesh
	// Statistics is communicated to MPI_MASTER_RANK, which writes it to a file
	// [TODO] When calculating total volume, calculate it separately for all tags and provide tags
	static void writeAnalyticTestResult(std::string & filename, MeshStatType & meshStatistics, MPIHelper & mpihelper)
	{
		 std::ofstream diagnosticFile;
		 if (mpihelper.rank() == MPI_MASTER_RANK)  { diagnosticFile.open(filename.c_str()); }

		 writeCommunicateVector(mpihelper, diagnosticFile, "NELEMENT",                       meshStatistics[0]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "NDBFACE",                        meshStatistics[1]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "NPBFACE",                        meshStatistics[2]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementShortestEdge",            meshStatistics[3]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementLongestEdge",             meshStatistics[4]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementLinearVolume",            meshStatistics[5]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementCurvilinearVolume",       meshStatistics[6]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementEnclosedSphereRadius",    meshStatistics[7]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementSurroundingSphereRadius", meshStatistics[8]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementLinearQuality1",          meshStatistics[9]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementLinearQuality2",          meshStatistics[10]);
		 writeCommunicateVector(mpihelper, diagnosticFile, "ElementCurvilinearQuality1",     meshStatistics[11]);
		 writeCommunicateSum(   mpihelper, diagnosticFile, "ProcssBoundarySurfaceArea",      meshStatistics[12][0]);
		 writeCommunicateSum(   mpihelper, diagnosticFile, "DomainBoundarySurfaceArea",      meshStatistics[13][0]);

		 if (mpihelper.rank() == MPI_MASTER_RANK)  { diagnosticFile.close(); }
	}


	// *****************************************************************
	// Auxiliary members
	// *****************************************************************

	// Initializes a FieldVector in 1 line
	static GlobalCoordinate initVector(ctype a, ctype b, ctype c)
	{
		GlobalCoordinate rez;
		rez[0] = a;
		rez[1] = b;
		rez[2] = c;
		return rez;
	}

	// Dot product between FieldVectors
	// [TODO] Replace with existing Dune functionality if found
	static ctype GridVectorDot(GlobalCoordinate a, GlobalCoordinate b)  { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

	// Cross product between FieldVectors
	// [TODO] Replace with existing Dune functionality if found
	static GlobalCoordinate GridVectorTimes(GlobalCoordinate a, GlobalCoordinate b)  { return initVector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]); }

	// Returns smallest entry of a vector of comparable objects
	template <class T>
	static T vectorMin (std::vector<T> data)
	{
		T rez = data[0];
		for (int i = 0; i < data.size(); i++)  { rez = std::min(rez, data[i]); }
		return rez;
	}

	// Returns largest entry of a vector of comparable objects
	template <class T>
	static T vectorMax (std::vector<T> data)
	{
		T rez = data[0];
		for (int i = 0; i < data.size(); i++)  { rez = std::max(rez, data[i]); }
		return rez;
	}


	// Collects data from all processes and writes it to a file on the master process
	template<class T>
	static void writeCommunicateVector(MPIHelper & mpihelper, std::ofstream & filestr, std::string title, std::vector<T> & data)
	{
		int rank = mpihelper.rank();
		int size = mpihelper.size();

		Dune::Communication<MPI_Comm> collective_comm = mpihelper.getCommunication();

		// 1) Communicate size of communication to root
		// *************************************************************
		int dataSize = data.size();
		std::vector<int> dataSizeArray(size);
		collective_comm.gather (&dataSize, reinterpret_cast<int*> (dataSizeArray.data()), 1, MPI_MASTER_RANK);


		// 2) Communicate actual data to root
		// *************************************************************
		int dataSizeTotal = 0;
		std::vector<int> displ;
		for (int i = 0; i < size; i++)
		{
			dataSizeTotal += dataSizeArray[i];
			displ.push_back((i == 0) ? 0  :  displ[i-1] + dataSizeArray[i-1]);
		}
		std::vector<T> dataTotal(dataSizeTotal);
		collective_comm.gatherv (data.data(), dataSize, reinterpret_cast<T*> (dataTotal.data()), dataSizeArray.data(), displ.data(), MPI_MASTER_RANK);


		// 3) Write output to file on Master process
		// *************************************************************
		if (rank == MPI_MASTER_RANK)
		{
			filestr << "<" << title << ">";
			for (int i = 0; i < dataTotal.size(); i++) { filestr << dataTotal[i] << " ";  }
			filestr << "</" << title << ">" << std::endl;
		}
	}


	// Sums data over all processes and writes result
	template <class T>
	static void writeCommunicateSum(MPIHelper & mpihelper, std::ofstream & filestr, std::string title, T data)
	{
		int rank = mpihelper.rank();
		int size = mpihelper.size();
		Dune::Communication<MPI_Comm> collective_comm = mpihelper.getCommunication();
		T sum = collective_comm.sum(data);
		if (rank == MPI_MASTER_RANK)  { filestr << "<" << title << ">" << sum << "</" << title << ">" << std::endl; }
	}
};

} // namespace CurvGrid

} // namespace Dune

#endif // DUNE_CURVILINEARDIAGNOSTICHELPER_HH
