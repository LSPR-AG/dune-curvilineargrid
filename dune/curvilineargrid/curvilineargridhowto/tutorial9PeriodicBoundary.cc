/********************************************
 * Dune-CurvilinearGrid
 * Tutorial 9: Periodic Boundary
 *
 * Author: Aleksejs Fomins
 *
 * Description: This tutorial demonstrates the current Periodic Boundary interface and capabilities
 * A 3D sine function is defined over the grid, propagation direction orthogonal to one of the coordinate axis.
 * Then, for two other axis, the value of the function should be symmetric on opposite faces of the cubic domain.
 * The difference of the function between neighboring periodic faces is integrated. It must amount to 0 for each
 * periodic direction. This tutorial checks that
 * 1) The mesh is periodic-conformal
 * 2) The orientation of the periodic face-neighbors is correct - the same face-local coordinate in each face corresponds to conformal global coordinates
 *
 * (This is not quite true - If viewed from left face, the right ghost-face is internally rotated to match the left one, and vice-versa. So, from point of view
 * of each intersection, both itself and neighbor have the same local index. However, each intersection itself has a different local index. Thus, to not
 * make errors, one must choose one of the periodic faces as principal, and only refer to the other one through its ghost neighbor)
 ********************************************/

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>


using namespace Dune;
using namespace Dune::CurvGrid;

const bool isCached = true;
const int DIM0D = 0;   const int CODIM0D = 3;
const int DIM1D = 1;   const int CODIM1D = 2;
const int DIM2D = 2;   const int CODIM2D = 1;
const int DIM3D = 3;   const int CODIM3D = 0;


template<class Grid, int mydim>
struct GlobalSinFunctor
{
	static const int cdim = Grid::dimension;
	typedef typename Grid::ctype ct;

	typedef Dune::FieldVector<ct, cdim>      GlobalCoordinate;
	typedef Dune::FieldVector<ct, mydim>  LocalCoordinate;
	typedef Dune::FieldVector<ct, cdim>  LocalCoordinateParent;

    static const unsigned int RETURN_SIZE = 1;
    typedef ct                    ResultValue;
    typedef std::vector<ResultValue>  ResultType;

	typedef typename Grid::Traits::LeafIntersection               Intersection;

	const Intersection & I_;
	ct k_;
	GlobalCoordinate evec_;


	// Instantiate with Intersection of interest.
	GlobalSinFunctor(const Intersection & I, ct k, GlobalCoordinate evec) : I_(I), k_(k), evec_(evec)  { }

	ct sineField(const GlobalCoordinate & x) const { return sin(k_ * (x * evec_)); }

	//
    ResultType operator()(const LocalCoordinate & xLocal) const
    {
    	LocalCoordinateParent xLocalInner = I_.geometryInInside().global(xLocal);
    	LocalCoordinateParent xLocalOuter = I_.geometryInOutside().global(xLocal);

    	GlobalCoordinate xGlobalInner = I_.inside().geometry().global(xLocalInner);
    	GlobalCoordinate xGlobalOuter = I_.outside().geometry().global(xLocalOuter);

    	ct fieldInner = sineField(xGlobalInner);
    	ct fieldOuter = sineField(xGlobalOuter);

    	return ResultType{fabs(fieldInner - fieldOuter)};
    	//return ResultType{(xGlobalInner - xGlobalOuter).two_norm()};
    }

    double zeroValue(unsigned int rezIndex) const { return 0.0; }
};


template<class GridType>
void IntegralPeriodicOrientationTest(const GridType & grid, double k, unsigned int kDir) {
	// Define standard grid properties
	const int dim =  GridType::dimension;
	typedef typename GridType::ctype ct;
	typedef typename GridType::LeafGridView LeafGridView;
	typedef typename LeafGridView::template Codim<CODIM2D>::Geometry FaceGeometry;
	typedef Dune::ReferenceElements<ct, DIM2D> ReferenceElements2d;

	// Define Integrator and Functor
	typedef GlobalSinFunctor<GridType, DIM2D> Integrand2D;
	typedef QuadratureIntegrator<ct, DIM2D>  Integrator2DScalar;
	typedef typename Integrator2DScalar::template Traits<Integrand2D>::StatInfo  StatInfo;

	// Choose integrator parameters
	double RELATIVE_TOLERANCE = 1.0e-5;
	double ACCURACY_GOAL = 1.0e-15;
	const int NORM_TYPE = QUADRATURE_NORM_L2;

	// Initialize the integral variable
	double maxIntegralError = 0.0;

	// Generate normal vector for the sinusoid direction
	typename Integrand2D::GlobalCoordinate edir(0.0);
	edir[kDir] = 1.0;



	// Iterate over entities of the grid, then over intersections (faces) of each entity
	LeafGridView leafView = grid.leafGridView();
	for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		for (auto&& intersection : intersections(leafView, elementThis)) {

			// Integrate only over periodic faces
			bool isPeriodic = ((intersection.neighbor() == true) && (intersection.boundary() == true));
			if (isPeriodic) {
				typename Integrand2D::GlobalCoordinate intersectionOuterNormal = intersection.unitOuterNormal( ReferenceElements2d::simplex().position( 0, 0 ) );

				int normalDirDim = -1;
				if (fabs(fabs(intersectionOuterNormal[0]) - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { normalDirDim = 0;}
				if (fabs(fabs(intersectionOuterNormal[1]) - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { normalDirDim = 1;}
				if (fabs(fabs(intersectionOuterNormal[2]) - 1.0) < NUMERICS_RELATIVE_TOLERANCE) { normalDirDim = 2;}
				assert(normalDirDim != -1);

				// We only want to integrate boundaries that are orthogonal to the plane wave, because for them the field should be symmetric for both periodic faces
				if (normalDirDim != kDir) {
					FaceGeometry geometry = intersection.geometry();
					Dune::GeometryType gt = intersection.type();

					Integrand2D sineFunctor(intersection, k, edir);

					StatInfo thisIntegralRez = Integrator2DScalar::template integrateRecursive<FaceGeometry, Integrand2D, NORM_TYPE>(geometry, sineFunctor, RELATIVE_TOLERANCE, ACCURACY_GOAL);

					std::stringstream logsstr;
					logsstr << "---- from entity " << gt;
					logsstr << " adding integral error contribution " << thisIntegralRez.second[0];
					logsstr << ". Needed order " << thisIntegralRez.first;
					LoggingMessage::template write<LOG_MSG_DVERB>(__FILE__, __LINE__, logsstr.str());

					maxIntegralError = std::max(maxIntegralError, thisIntegralRez.second[0]);
				}
			}
		}
	}

	// The actual integral is the sum over all processors
	double rez = grid.comm().max(maxIntegralError);
	if (grid.comm().rank() == 0)  { std::cout << "Maximal integral error in direction " << kDir << " amounted to " << rez << std::endl; }
}










int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

	// Create Grid
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, argc, argv);


    // *********************************************************
	// Check if grid recognises it is parallel
    // *********************************************************
	std::stringstream loggingStream;
	loggingStream << "Is curvilinear grid periodic = " << grid->withPeriodic();
	LoggingMessage::template write<LOG_MSG_PRODUCTION>(__FILE__, __LINE__, loggingStream.str());

	// [FIXME] Only run tests on those dimensions that are specified to be periodic
	double k = 1;
	IntegralPeriodicOrientationTest(*grid, k, 0);
	IntegralPeriodicOrientationTest(*grid, k, 1);
	IntegralPeriodicOrientationTest(*grid, k, 2);

	// Report the parallel timing statistics
	LoggingTimer<LoggingMessage>::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}
