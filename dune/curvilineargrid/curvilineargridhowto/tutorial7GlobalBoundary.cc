/********************************************
 * Dune-CurvilinearGrid
 * Tutorial 7: Global Boundary Container - Domain Boundary
 *
 * Author: Aleksejs Fomins
 *
 * Description: This tutorial demonstrates the current capabilities for Boundary Integral methods.
 * The boundary container is used to gather all domain boundary segments on each process.
 * Then each process performs normal integration presented in tutorial 5b, and checks that
 * it amounts to approx. zero, thus confirming that the global boundary is closed on each process (all boundary segments are present)
 ********************************************/


#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits.h>


#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>

#include <dune/curvilineargrid/utility/globalboundarycontainer.hh>
#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>



using namespace Dune;
using namespace Dune::CurvGrid;

// Calculates the outer normal to the intersection times the integration element
template<class Grid, int mydim, class SegmentContainer>
struct NormalFunctor
{
	typedef Dune::FieldVector<double, mydim>           LocalCoordinate;
	typedef typename Grid::Traits::LeafIntersection    Intersection;
	typedef Dune::FieldVector<double, 3>  GlobalCoordinate;

    static const unsigned int RETURN_SIZE = 1;
    typedef GlobalCoordinate               ResultValue;
    typedef std::vector<GlobalCoordinate>  ResultType;

    const SegmentContainer & s_;

	NormalFunctor(const SegmentContainer & s) : s_(s)  { }

	// Calculates the outer normal to the intersection times the integration element
	ResultType operator()(const LocalCoordinate & x) const
    {
    	return ResultType(1, s_.unitOuterNormal(x));
    }

	GlobalCoordinate zeroValue(unsigned int rezIndex) const { return GlobalCoordinate(0.0); }
};


// This method iterates over all elements of the grid
// Then over all intersections
// Then for each intersection that is a boundary segment it integrates the normal component by component
// [TODO] The ultimate test: For every vertex, store global index AND its coordinate. When found vertex with same global index from another communicated face, check that the coordinate matches
template<class GridType>
class TestIntegrals
{

	static const int ELEMENT_CODIM = 0;
	static const int FACE_CODIM = 1;
	static const int EDGE_CODIM = 2;
	static const int VERTEX_CODIM = 3;

	typedef unsigned int UInt;

public:
	static const int dim =  GridType::dimension;
	typedef typename GridType::ctype ct;
	typedef typename GridType::LeafGridView LeafGridView;
	typedef typename LeafGridView::IntersectionIterator IntersectionIterator;
	typedef typename LeafGridView::template Codim< ELEMENT_CODIM >::Entity Entity;
	typedef typename IntersectionIterator :: Intersection Intersection;

	typedef typename LeafGridView::template Codim<FACE_CODIM>::Entity			EntityFace;
	typedef typename LeafGridView::template Codim<EDGE_CODIM>::Entity			EntityEdge;
	typedef typename LeafGridView::template Codim<VERTEX_CODIM>::Entity		EntityVertex;

	typedef typename LeafGridView::template Codim<ELEMENT_CODIM>::Iterator EntityLeafIterator;
	typedef typename LeafGridView::template Codim<FACE_CODIM>::Geometry FaceGeometry;
	typedef typename LeafGridView::template Codim<EDGE_CODIM>::Geometry EdgeGeometry;

	typedef typename GridType::template Codim<FACE_CODIM>::EntityGeometryMappingImpl  BaseGeometryFace;
	typedef typename GridType::template Codim<EDGE_CODIM>::EntityGeometryMappingImpl  BaseGeometryEdge;

	typedef typename FaceGeometry::LocalCoordinate  LocalCoordinate;
	typedef typename FaceGeometry::GlobalCoordinate  GlobalCoordinate;

	typedef Dune::ReferenceElements<ct, dim> ReferenceElements3d;
	typedef Dune::ReferenceElements<ct, dim-1> ReferenceElements2d;

	typedef GlobalBoundaryContainer<GridType> BoundaryContainer;
	typedef GlobalBoundaryIterator<GridType> BoundaryIterator;
	//typedef typename BoundaryContainer::BoundaryContainer  SegmentContainer;

	typedef NormalFunctor<GridType, 2, Intersection> Integrand2DVectorLocal;
	typedef NormalFunctor<GridType, 2, BoundaryIterator> Integrand2DVectorParallel;
	typedef QuadratureIntegrator<double, 2>  Integrator2DVector;
	typedef typename Integrator2DVector::template Traits<Integrand2DVectorLocal>::StatInfo  StatInfo;

	typedef std::map<int, int> GIndMap;
	typedef std::pair<int,int> GIndPair;
	typedef typename GIndMap::iterator  GIndMapIter;

public:


	TestIntegrals() {}

	static bool isBoundary(const GridType & grid, const Entity & elem, const Intersection & intr, bool isDB, int volTag = 0, int surfTag = 0) {
		EntityFace faceBra = elem.template subEntity<FACE_CODIM>(intr.indexInInside());
		int elemTag = grid.template entityPhysicalTag<ELEMENT_CODIM>(elem);
		int boundaryTag = grid.template entityPhysicalTag<FACE_CODIM>(faceBra);
		return isDB
			  ? (intr.boundary() == true) && (intr.neighbor() == false)
			  : (elemTag == volTag) && (boundaryTag == surfTag);
	}


	static GlobalCoordinate normalIntegralSelf (const GridType & grid, bool isDB, int volTag = 0, int surfTag = 0, int normalSign = 1) {
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, ":::Testing Normal Integral on own boundary:::");

		// get the instance of the LeafGridView
		LeafGridView leafView = grid.leafGridView();
		const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

		GlobalCoordinate  rez(0.0);
		LocalCoordinate faceCenterLocal = ReferenceElements2d::simplex().position( 0, 0 );

		for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder)) {
			//std::cout << "-accessing entity " << indexSet.index(elemThis) << std::endl;

			for (auto&& intersection : intersections(leafView, elemThis)) {

				// Only process the intersection if it is a valid interior or domain boundary segment
				bool isBS = isBoundary(grid, elemThis, intersection, isDB, volTag, surfTag);
				if (isBS) {
					//Dune::GeometryType gt = intersection.type();
					FaceGeometry geometry = intersection.geometry();
					//GlobalCoordinate normal = intersection.unitOuterNormal(faceCenterLocal);

					//Integrand2DVector integrand(intersection);
					Integrand2DVectorLocal integrand(intersection);

					const double RELATIVE_TOLERANCE = 1.0e-5;
					const double ACCURACY_GOAL = 1.0e-15;
					const int NORM_TYPE = QUADRATURE_NORM_L2;
					StatInfo thisIntegralN = Integrator2DVector::template integrateRecursive<FaceGeometry, Integrand2DVectorLocal, NORM_TYPE>(geometry, integrand, RELATIVE_TOLERANCE, ACCURACY_GOAL);

					//std::cout << "---- adding normal contribution " << thisIntegralN.second[0] << " from " << gt << ". Needed order " << thisIntegralN.first << std::endl;

					rez += thisIntegralN.second[0];
				}
		  }
	  }

	  // Note that the normal may be inner or outer wrt associated volume element, determined by user
	  rez *= normalSign;

	  return rez;
	}


	static GlobalCoordinate normalIntegralOtherDomain(const BoundaryContainer & container, int normalSign = 1) {
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, ":::Testing Normal Integral on global boundary container:::");

		GlobalCoordinate  rez(0.0);

		BoundaryIterator biter(container);
		while (!biter.end()) {
			Integrand2DVectorParallel integrand(biter);

			const double RELATIVE_TOLERANCE = 1.0e-5;
			const double ACCURACY_GOAL = 1.0e-15;
			const int NORM_TYPE = QUADRATURE_NORM_L2;
			StatInfo thisIntegralN = Integrator2DVector::template integrateRecursive<BaseGeometryFace, Integrand2DVectorParallel, NORM_TYPE>(
					biter.geometry(), integrand, RELATIVE_TOLERANCE, ACCURACY_GOAL);

			//std::cout << "---- adding normal contribution " << thisIntegralN.second[0] << " from " << baseGeom.type() << ". Needed order " << thisIntegralN.first << std::endl;

			rez += thisIntegralN.second[0];

			++biter;
		}

		// Note that the normal may be inner or outer wrt associated volume element, determined by user
		rez *= normalSign;
		return rez;
	}


	static ct edgeLengthSelf(const GridType & grid, bool isDB, int volTag = 0, int surfTag = 0) {
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, ":::Testing Edge Length on own boundary:::");

		// get the instance of the LeafGridView
		LeafGridView leafView = grid.leafGridView();
		const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

		ct rez(0.0);

		for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder)) {
			//std::cout << "-accessing entity " << indexSet.index(elemThis) << std::endl;

			for (auto&& intersection : intersections(leafView, elemThis)) {

				// Only process the intersection if it is a valid interior or domain boundary segment
				bool isBS = isBoundary(grid, elemThis, intersection, isDB, volTag, surfTag);
				if (isBS) {
					for (int iEdge = 0; iEdge < 3; iEdge++) {
						// Get subindex of edge inside element
						int edgeIndInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intersection.indexInInside(), FACE_CODIM, iEdge, EDGE_CODIM);
						EdgeGeometry edgeGeom = elemThis.template subEntity<EDGE_CODIM>(edgeIndInElem).geometry();
						rez += edgeGeom.volume();
					}
				}
			}
		}
		return rez * 0.5; // Count each edge only 1/2 since it is accessed twice, once from each neighboring face
	}


	static ct edgeLengthOtherDomain(const BoundaryContainer & container) {
		LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, ":::Testing Edge Length on global boundary container:::" );

		const double RELATIVE_TOLERANCE = 1.0e-5;

		ct  rez(0.0);

		BoundaryIterator biter(container);
		while (!biter.end()) {
			for (int iEdge = 0; iEdge < 3; iEdge++) {
				rez += biter.geometryEdge(iEdge).volume(RELATIVE_TOLERANCE);
			}
			++biter;
		}

		return rez * 0.5; // Count each edge only 1/2 since it is accessed twice, once from each neighboring face
	}


	static void globalIndexSelf(
		  const GridType & grid,
		  GIndMap & gindmapface,
		  GIndMap & gindmapedge,
		  GIndMap & gindmapvertex,
		  bool isDB, int volTag = 0, int surfTag = 0
	) {
	  LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, ":::Testing Global Index on own boundary:::");


	  // get the instance of the LeafGridView
	  LeafGridView leafView = grid.leafGridView();
	  const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

	  int elemCount = 0;
	  unsigned int nElementInterior = grid.numInternal(ELEMENT_CODIM);
	  for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		  //std::cout << "-accessing entity " << indexSet.index(elemThis) << std::endl;
		  LoggingMessage::writePatience("test: globalIndexSelf...", elemCount++, nElementInterior);

		  for (auto&& intersection : intersections(leafView, elemThis)) {

			  // Only process the intersection if it is a valid interior or domain boundary segment
			  bool isBS = isBoundary(grid, elemThis, intersection, isDB, volTag, surfTag);
			  if (isBS) {
				  UInt intrIndexInInside = intersection.indexInInside();
				  EntityFace faceThis = elemThis.template subEntity<FACE_CODIM>(intrIndexInInside);
				  int gindface = grid.template entityGlobalIndex<FACE_CODIM>(faceThis);

				  assert(gindmapface.find(gindface) == gindmapface.end());  // The face global index may not repeat
				  gindmapface.insert(GIndPair(gindface, 1));

				  for (int iEdge = 0; iEdge < 3; iEdge++) {
					UInt indexEdgeInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intrIndexInInside, FACE_CODIM, iEdge, EDGE_CODIM);
					EntityEdge edgeThis = elemThis.template subEntity<EDGE_CODIM>(indexEdgeInElem);
					int gindedge = grid.template entityGlobalIndex<EDGE_CODIM>(edgeThis);

					GIndMapIter edgeptr = gindmapedge.find(gindedge);

					if (edgeptr == gindmapedge.end()) {
						gindmapedge.insert(GIndPair(gindedge, 1));
					} else {
						//std::cout << "found copy 1 with size " << edgeptr->second << std::endl;

						  if (edgeptr->second >= 2) {  // Each edge should be accessed at most twice
							  std::cout << "Bug: nEdge=" << edgeptr->second + 1 << " @ index=" << edgeptr->first << std::endl;
							  DUNE_THROW(Dune::IOError, "bug");
						  }
						  edgeptr->second += 1;
					}
				  }

				  for (int iCorner = 0; iCorner < 3; iCorner++) {
						UInt indexCornerInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intrIndexInInside, FACE_CODIM, iCorner, VERTEX_CODIM);
						EntityVertex cornerThis = elemThis.template subEntity<VERTEX_CODIM>(indexCornerInElem);
						int gindcorner = grid.template entityGlobalIndex<VERTEX_CODIM>(cornerThis);

						GIndMapIter cornerptr = gindmapvertex.find(gindcorner);
						if (cornerptr == gindmapvertex.end()) {	gindmapvertex.insert(GIndPair(gindcorner, 1)); }
						else {															cornerptr->second += 1; }
				  }
			  }
		  }
	  }
	}


	static void globalIndexOtherDomain(
		  const BoundaryContainer & container,
		  GIndMap & gindmapface,
		  GIndMap & gindmapedge,
		  GIndMap & gindmapvertex
	) {
	  LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, ":::Testing GlobalIndex on global boundary container:::");

	  BoundaryIterator biter(container);
	  while (!biter.end()) {
		  int gindface = biter.template globalIndex<FACE_CODIM>(0);
		  assert(gindmapface.find(gindface) == gindmapface.end());  // The face global index may not repeat
		  gindmapface.insert(GIndPair(gindface, 1));

		  for (int iEdge = 0; iEdge < 3; iEdge++) {
			  int gindedge = biter.template globalIndex<EDGE_CODIM>(iEdge);
			  GIndMapIter edgeptr = gindmapedge.find(gindedge);

				if (edgeptr == gindmapedge.end()) {
					gindmapedge.insert(GIndPair(gindedge, 1));
				} else {
					//std::cout << "found copy 2" << std::endl;

					  if (edgeptr->second >= 2) {  // Each edge should be accessed at most twice
						  std::cout << "Bug: nEdge=" << edgeptr->second + 1 << " @ index=" << edgeptr->first << std::endl;
						  DUNE_THROW(Dune::IOError, "bug");
					  }
					  edgeptr->second += 1;
				}
		  }

		  for (int iCorner = 0; iCorner < 3; iCorner++) {
			  int gindcorner = biter.template globalIndex<VERTEX_CODIM>(iCorner);
			  GIndMapIter cornerptr = gindmapvertex.find(gindcorner);

			  if (cornerptr == gindmapvertex.end()) {	gindmapvertex.insert(GIndPair(gindcorner, 1)); }
			  else {															cornerptr->second += 1; }
		  }

		  ++biter;
	  }
	}
};




/** \brief
 * NOTE: ONLY SELECT LINEAR FILES, FOR NOW BOUNDARY COMMUNICATOR CAN NOT DEAL WITH CURVILINEAR
 * NOTE: Currently only works with grid number 7, as boundary tags are hard-coded [TODO] Fix this
 *
 */


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	const bool isCached = true;
	const int ELEMENT_CODIM = 0;  // Codimension of element in 3D
	const int FACE_CODIM = 1;  // Codimension of face in 3D

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, argc, argv);


	typedef GlobalBoundaryContainer<GridType> BoundaryContainer;
	//BoundaryContainer testContainer(*grid, true);

	bool isDomainBoundary = false;
	int volumeTag = 501;
	int surfaceTag = 101;
	int normalSign = -1;
	BoundaryContainer testContainer(*grid, isDomainBoundary, volumeTag, surfaceTag);

	std::cout << ":::Constructed Container:::" << std::endl;

	typedef TestIntegrals<GridType> Integr;
	typename Integr::GlobalCoordinate rezNormalSelf = Integr::normalIntegralSelf(*grid, isDomainBoundary, volumeTag, surfaceTag, normalSign);
	typename Integr::GlobalCoordinate rezNormalOther = Integr::normalIntegralOtherDomain(testContainer, normalSign);
	ctype rezEdgeSelf = Integr::edgeLengthSelf(*grid, isDomainBoundary, volumeTag, surfaceTag);
	ctype rezEdgeOther = Integr::edgeLengthOtherDomain(testContainer);

	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
	MPI_Barrier(comm);  // Wait, to have all processes write next to each other
	std::cout << "Normal: On rank " << grid->comm().rank() << " self=" << rezNormalSelf << ", other=" << rezNormalOther << ", sum=" << rezNormalSelf+rezNormalOther << std::endl;
	MPI_Barrier(comm);  // Wait, to have all processes write next to each other
	std::cout << "EdgeLen: On rank " << grid->comm().rank() << " self=" << rezEdgeSelf << ", other=" << rezEdgeOther << ", sum=" << rezEdgeSelf+rezEdgeOther << std::endl;

	/////////////////////////////
	// Test 1
	/////////////////////////////

	std::map<int, int> globalIndexFace;
	std::map<int, int> globalIndexEdge;
	std::map<int, int> globalIndexVertex;
	Integr::globalIndexSelf(*grid, globalIndexFace, globalIndexEdge, globalIndexVertex, isDomainBoundary, volumeTag, surfaceTag);
	Integr::globalIndexOtherDomain(testContainer, globalIndexFace, globalIndexEdge, globalIndexVertex);

	int faceMin = INT_MAX;
	int faceMax = INT_MIN;
	int edgeMin = INT_MAX;
	int edgeMax = INT_MIN;
	int cornerMin = INT_MAX;
	int cornerMax = INT_MIN;

	typedef typename TestIntegrals<GridType>::GIndMapIter GIndMapIter;
	for (GIndMapIter iter = globalIndexFace.begin(); iter != globalIndexFace.end(); iter++) {
		faceMin = std::min(faceMin, iter->first);
		faceMax = std::max(faceMax, iter->first);
	}
	for (GIndMapIter iter = globalIndexEdge.begin(); iter != globalIndexEdge.end(); iter++) {
		edgeMin = std::min(edgeMin, iter->first);
		edgeMax = std::max(edgeMax, iter->first);
		//std::cout << "Index=" << iter->first << " nEdge=" << iter->second << std::endl;
		assert(iter->second == 2);
	}
	for (GIndMapIter iter = globalIndexVertex.begin(); iter != globalIndexVertex.end(); iter++) {
		cornerMin = std::min(cornerMin, iter->first);
		cornerMax = std::max(cornerMax, iter->first);
	}

	MPI_Barrier(comm);  // Wait, to have all processes write next to each other

	std::cout
		<< "nBoundary=" << globalIndexFace.size() << " min=" << faceMin << " max=" << faceMax
		<< "; nEdge=" << globalIndexEdge.size() << " min=" << edgeMin << " max=" << edgeMax
		<< "; nCorner=" << globalIndexVertex.size() << " min=" << cornerMin << " max=" << cornerMax << std::endl;


	/////////////////////////////
	// Test 2
	/////////////////////////////


	/*
	std::map<int, int> globalIndexFace2self;
	std::map<int, int> globalIndexEdge2self;
	std::map<int, int> globalIndexVertex2self;
	Integr::globalIndexSelf(*grid, globalIndexFace2self, globalIndexEdge2self, globalIndexVertex2self);
	for (int i = 0; i < grid->comm().size(); i++) {
		if (i != grid->comm().rank())  { MPI_Barrier(MPI_COMM_WORLD); }  // If its not your turn, you wait
		else {                                                          // If it is your turn, you write to file
			std::cout << "GInt self on process " << grid->comm().rank() << ": ";
			for (GIndMapIter iter = globalIndexFace2self.begin(); iter != globalIndexFace2self.end(); iter++) { std::cout << iter->first << " "; }
			std::cout << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);  // Note: MPI_Barrier unlocks, when all processes have called it
		}
	} */

	std::map<int, int> globalIndexFace2;
	std::map<int, int> globalIndexEdge2;
	std::map<int, int> globalIndexVertex2;

	Integr::globalIndexOtherDomain(testContainer, globalIndexFace2, globalIndexEdge2, globalIndexVertex2);

	/*
	for (int i = 0; i < grid->comm().size(); i++) {
		if (i != grid->comm().rank())  { MPI_Barrier(MPI_COMM_WORLD); }  // If its not your turn, you wait
		else {                                                          // If it is your turn, you write to file
			std::cout << "GInt on process " << grid->comm().rank() << ": ";
			for (GIndMapIter iter = globalIndexFace2.begin(); iter != globalIndexFace2.end(); iter++) { std::cout << iter->first << " "; }
			std::cout << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);  // Note: MPI_Barrier unlocks, when all processes have called it
		}
	}*/


	int tmpLocalIndex;
	for (GIndMapIter iter = globalIndexFace2.begin(); iter != globalIndexFace2.end(); iter++) {
		//std::cout << "Process "<< grid->comm().rank() << " testing global index " << iter->first << std::endl;

		// Mapped global boundary segments must not exist on this process
		bool local = grid->entityLocalIndex(FACE_CODIM, iter->first, tmpLocalIndex);

		if (local) {
			Dune::PartitionType thisType = grid->gridbase().entityPartitionType(FACE_CODIM, tmpLocalIndex);
			if (thisType != Dune::PartitionType::GhostEntity) {
				std::cout << "...Warning: Unexpexted parallel face of local index " << tmpLocalIndex << " and type " << thisType << std::endl;
			}
		}
	}



	typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
