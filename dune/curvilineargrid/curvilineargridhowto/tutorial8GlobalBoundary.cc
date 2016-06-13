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


// Calculates the outer normal to the intersection times the integration element
template<class Grid, int mydim>
struct NormalFunctorLinear
{
	typedef Dune::FieldVector<double, mydim>           LocalCoordinate;
	typedef typename Grid::Traits::LeafIntersection    Intersection;
	typedef Dune::FieldVector<double, 3>  GlobalCoordinate;

    static const unsigned int RETURN_SIZE = 1;
    typedef GlobalCoordinate               ResultValue;
    typedef std::vector<GlobalCoordinate>  ResultType;

	//Intersection I_;
    GlobalCoordinate n_;

	NormalFunctorLinear(const GlobalCoordinate & n) : n_(n)  { }

	// Calculates the outer normal to the intersection times the integration element
	ResultType operator()(const LocalCoordinate & x) const
    {
    	return ResultType(1, n_);
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

  typedef NormalFunctorLinear<GridType, 2>             Integrand2DVector;
  typedef Dune::QuadratureIntegrator<double, 2>  Integrator2DVector;
  typedef typename Integrator2DVector::template Traits<Integrand2DVector>::StatInfo  StatInfo;

  typedef Dune::ReferenceElements<ct, dim> ReferenceElements3d;
  typedef Dune::ReferenceElements<ct, dim-1> ReferenceElements2d;

  typedef Dune::CurvGrid::GlobalBoundaryContainer<GridType> BoundaryContainer;
  typedef Dune::CurvGrid::GlobalBoundaryIterator<GridType> BoundaryIterator;

  typedef std::map<int, int> GIndMap;
  typedef std::pair<int,int> GIndPair;
  typedef typename GIndMap::iterator  GIndMapIter;

public:


  TestIntegrals() {}


  static GlobalCoordinate normalIntegralSelf (const GridType & grid) {
	  // get the instance of the LeafGridView
	  LeafGridView leafView = grid.leafGridView();
	  const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

	  GlobalCoordinate  rez(0.0);
	  LocalCoordinate faceCenterLocal = ReferenceElements2d::simplex().position( 0, 0 );

	  for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		  //std::cout << "-accessing entity " << indexSet.index(elemThis) << std::endl;

		  for (auto&& intersection : intersections(leafView, elemThis)) {
			  if ( (intersection.boundary() == true) && (intersection.neighbor() == false) ) {
				  Dune::GeometryType gt = intersection.type();
				  FaceGeometry geometry = intersection.geometry();
				  GlobalCoordinate normal = intersection.unitOuterNormal(faceCenterLocal);

				  //Integrand2DVector integrand(intersection);
				  Integrand2DVector integrand(normal);

				  const double RELATIVE_TOLERANCE = 1.0e-5;
				  const double ACCURACY_GOAL = 1.0e-15;
				  const int NORM_TYPE = Dune::QUADRATURE_NORM_L2;
				  StatInfo thisIntegralN = Integrator2DVector::template integrateRecursive<FaceGeometry, Integrand2DVector, NORM_TYPE>(geometry, integrand, RELATIVE_TOLERANCE, ACCURACY_GOAL);

				  //std::cout << "---- adding normal contribution " << thisIntegralN.second[0] << " from " << gt << ". Needed order " << thisIntegralN.first << std::endl;

				  rez += thisIntegralN.second[0];
				}
		  }
	  }
	  return rez;
  }


  static GlobalCoordinate normalIntegralOtherDomain(const BoundaryContainer & container) {

	  GlobalCoordinate  rez(0.0);

	  BoundaryIterator biter(container);
	  while (!biter.end()) {

		  BaseGeometryFace baseGeom = biter.geometry();
		  GlobalCoordinate normal = biter.unitOuterNormal();

		  Integrand2DVector integrand(normal);

		  const double RELATIVE_TOLERANCE = 1.0e-5;
		  const double ACCURACY_GOAL = 1.0e-15;
		  const int NORM_TYPE = Dune::QUADRATURE_NORM_L2;
		  StatInfo thisIntegralN = Integrator2DVector::template integrateRecursive<BaseGeometryFace, Integrand2DVector, NORM_TYPE>(baseGeom, integrand, RELATIVE_TOLERANCE, ACCURACY_GOAL);

		  //std::cout << "---- adding normal contribution " << thisIntegralN.second[0] << " from " << baseGeom.type() << ". Needed order " << thisIntegralN.first << std::endl;

		  rez += thisIntegralN.second[0];

		  ++biter;
	  }

	  return rez;
  }


  static ct edgeLengthSelf(const GridType & grid) {
	  // get the instance of the LeafGridView
	  LeafGridView leafView = grid.leafGridView();
	  const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

	  ct rez(0.0);

	  for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		  //std::cout << "-accessing entity " << indexSet.index(elemThis) << std::endl;

		  for (auto&& intersection : intersections(leafView, elemThis)) {
			  if ( (intersection.boundary() == true) && (intersection.neighbor() == false) ) {
				  for (int iEdge = 0; iEdge < 3; iEdge++) {
					  // Get subindex of edge inside element
					  int edgeIndInElem = ReferenceElements3d::general(elemThis.type()).subEntity(intersection.indexInInside(), FACE_CODIM, iEdge, EDGE_CODIM);
					  EdgeGeometry edgeGeom = elemThis.template subEntity<EDGE_CODIM>(edgeIndInElem).geometry();
					  rez += edgeGeom.volume();
				  }
			  }
		  }
	  }
	  return rez * 0.5;  // Count each edge only 1/2 since it is accessed twice, once from each neighboring face
  }


  static ct edgeLengthOtherDomain(const BoundaryContainer & container) {
	  const double RELATIVE_TOLERANCE = 1.0e-5;

	  ct  rez(0.0);

	  BoundaryIterator biter(container);
	  while (!biter.end()) {
		  for (int iEdge = 0; iEdge < 3; iEdge++) {

			  rez += biter.geometryEdge(iEdge).volume(RELATIVE_TOLERANCE);
		  }
		  ++biter;
	  }

	  return rez * 0.5;
  }


  static void globalIndexSelf(
		  const GridType & grid,
		  GIndMap & gindmapface,
		  GIndMap & gindmapedge,
		  GIndMap & gindmapvertex
  ) {
	  // get the instance of the LeafGridView
	  LeafGridView leafView = grid.leafGridView();
	  const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

	  int elemCount = 0;
	  unsigned int nElementInterior = grid.numInternal(ELEMENT_CODIM);
	  for (auto&& elemThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		  //std::cout << "-accessing entity " << indexSet.index(elemThis) << std::endl;

		  Dune::LoggingMessage::writePatience("test: globalIndexSelf...", elemCount++, nElementInterior);



		  for (auto&& intersection : intersections(leafView, elemThis)) {
			  if ( (intersection.boundary() == true) && (intersection.neighbor() == false) ) {
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
 *
 */


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	const int grid_file_type = 6;  // createGrid procedure provides 7 different example grids numbered 0 to 6
	// NOTE: ONLY SELECT LINEAR FILES, FOR NOW BOUNDARY COMMUNICATOR CAN NOT DEAL WITH CURVILINEAR

	const bool isCached = true;
	const int ELEMENT_CODIM = 0;  // Codimension of element in 3D

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, grid_file_type);


	typedef Dune::CurvGrid::GlobalBoundaryContainer<GridType> BoundaryContainer;
	BoundaryContainer testContainer(*grid);

	typedef TestIntegrals<GridType> Integr;
	typename Integr::GlobalCoordinate rezNormalSelf = Integr::normalIntegralSelf(*grid);
	typename Integr::GlobalCoordinate rezNormalOther = Integr::normalIntegralOtherDomain(testContainer);
	ctype rezEdgeSelf = Integr::edgeLengthSelf(*grid);
	ctype rezEdgeOther = Integr::edgeLengthOtherDomain(testContainer);

	MPI_Comm comm = Dune::MPIHelper::getCommunicator();
	MPI_Barrier(comm);  // Wait, to have all processes write next to each other
	std::cout << "Normal: On rank " << grid->comm().rank() << " self=" << rezNormalSelf << ", other=" << rezNormalOther << ", sum=" << rezNormalSelf+rezNormalOther << std::endl;
	MPI_Barrier(comm);  // Wait, to have all processes write next to each other
	std::cout << "EdgeLen: On rank " << grid->comm().rank() << " self=" << rezEdgeSelf << ", other=" << rezEdgeOther << ", sum=" << rezEdgeSelf+rezEdgeOther << std::endl;


	std::map<int, int> globalIndexFace;
	std::map<int, int> globalIndexEdge;
	std::map<int, int> globalIndexVertex;
	Integr::globalIndexSelf(*grid, globalIndexFace, globalIndexEdge, globalIndexVertex);
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


	typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}