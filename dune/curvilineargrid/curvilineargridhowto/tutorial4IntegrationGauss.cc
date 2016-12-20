// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>

#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>



using namespace Dune;
using namespace Dune::CurvGrid;

const bool isCached = true;
const int DIM0D = 0;   const int CODIM0D = 3;
const int DIM1D = 1;   const int CODIM1D = 2;
const int DIM2D = 2;   const int CODIM2D = 1;
const int DIM3D = 3;   const int CODIM3D = 0;


// This Functor creates a field of a single charge at the origin of global coordinates
// Then calculates the projection of this field on the local surface times the integration element
// This functor is meant to be directly integrated over reference element to calculate the contribution
// to Gauss law from a boundary surface
template<class Grid, int mydim>
struct GaussFunctor
{
	static const int cdim =  Grid::dimension;
	typedef typename Grid::ctype ct;

	typedef Dune::FieldVector<ct, cdim>      GlobalCoordinate;
	typedef Dune::FieldVector<ct, mydim>  LocalCoordinate;

    static const unsigned int RETURN_SIZE = 1;
    typedef double                    ResultValue;
    typedef std::vector<ResultValue>  ResultType;

	typedef typename Grid::Traits::LeafIntersection               Intersection;

	const Intersection & I_;
	const GlobalCoordinate & x0_;

	// Instantiate with Intersection of interest.
    GaussFunctor(const Intersection & I, const GlobalCoordinate & x0) : I_(I), x0_(x0)  { }

    // Calculates EM field vector of a single charge assuming the charge is at the origin of global coordinates
    // Also assuming dielectric permittivity epsilon = 1
    GlobalCoordinate ChargeField(const GlobalCoordinate & x) const
    {
    	GlobalCoordinate rez = x - x0_;
    	rez /= pow(rez.two_norm2(), 1.5);
    	return rez;
    }

    // For a local coordinate on a surface of boundary segment
    // 1) Finds integration outer normal at that point
    // 2) Finds the associated global coordinate
    // 3) Finds the EM field at that coordinate
    // 4) returns scalar product between the normal and the field, thus the integrand for reference integration
    ResultType operator()(const LocalCoordinate & x) const
    {
    	GlobalCoordinate normal = I_.unitOuterNormal(x);
    	GlobalCoordinate global = I_.geometry().global(x);
    	GlobalCoordinate field = ChargeField(global);

    	return ResultType(1, normal * field);
    }

    double zeroValue(unsigned int rezIndex) const { return 0.0; }
};


// This method iterates over all elements of the grid
// Then over all intersections
// Then for each intersection that is a boundary segment it integrates the Surface projection of EM field of a single charge
// The output of this function should be 4pi for whatever closed geometry that includes the origin
template<class GridType, class GlobalCoordinate>
void Integrate (GridType& grid, const GlobalCoordinate & x0)
{
	// Define standard grid properties
	const int dim =  GridType::dimension;
	typedef typename GridType::ctype ct;
	typedef typename GridType::LeafGridView LeafGridView;
	typedef typename LeafGridView::template Codim<CODIM2D>::Geometry FaceGeometry;

	// Define Integrator and Functor
	typedef GaussFunctor<GridType, DIM2D>              Integrand2D;
	typedef QuadratureIntegrator<ct, DIM2D>  Integrator2DScalar;
	typedef typename Integrator2DScalar::template Traits<Integrand2D>::StatInfo  StatInfo;

	// Choose integrator parameters
	double RELATIVE_TOLERANCE = 1.0e-5;
	double ACCURACY_GOAL = 1.0e-15;
	const int NORM_TYPE = QUADRATURE_NORM_L2;

	// Initialize the integral variable
	double gaussintegral = 0.0;


	// Iterate over entities of the grid, then over intersections (faces) of each entity
	LeafGridView leafView = grid.leafGridView();
	for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		for (auto&& intersection : intersections(leafView, elementThis)) {
			// Check if intersection is a domain boundary segment - has no neighbor element
			if (!intersection.neighbor()) {
			  FaceGeometry geometry = intersection.geometry();
			  Dune::GeometryType gt = intersection.type();

			  Integrand2D gaussf(intersection, x0);

			  StatInfo thisIntegralG = Integrator2DScalar::template integrateRecursive<FaceGeometry, Integrand2D, NORM_TYPE>(geometry, gaussf, RELATIVE_TOLERANCE, ACCURACY_GOAL);

			  std::stringstream logsstr;
			  logsstr << "---- from entity " << gt;
			  logsstr << " adding Gauss integral contribution " << thisIntegralG.second[0];
			  logsstr << ". Needed order " << thisIntegralG.first;
			  LoggingMessage::template write<LOG_MSG_PRODUCTION>(__FILE__, __LINE__, logsstr.str());

			  gaussintegral += thisIntegralG.second[0];
			}
		}
	}

  // The actual integral is the sum over all processors
  double rez = grid.comm().sum(gaussintegral);
  if (grid.comm().rank() == 0)  { std::cout << "Gauss integral amounted to " << rez << std::endl; }

}



int main (int argc , char **argv) {
	typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, argc, argv);

	// Allow parallel output over all processes for the remainder of the program
	LoggingMessage::setPVerbose(true);

	// Place a charge slightly off of the origin, to make the field on the spherical boundary asymmetric
	Dune::FieldVector<ctype, dim> x0;  x0[0] = 0.2;

	// Perform the integration
	LoggingTimerDev::time("Integrating Gauss integral");
	Integrate(*grid, x0);
	LoggingTimerDev::time("Integrating Gauss integral");

	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}

