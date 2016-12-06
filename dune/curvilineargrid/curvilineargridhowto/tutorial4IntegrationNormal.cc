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




const bool isCached = true;
const int DIM0D = 0;   const int CODIM0D = 3;
const int DIM1D = 1;   const int CODIM1D = 2;
const int DIM2D = 2;   const int CODIM2D = 1;
const int DIM3D = 3;   const int CODIM3D = 0;


// Calculates the outer normal to the intersection times the integration element
template<class Grid, int mydim>
struct NormalFunctor
{
	static const int cdim =  Grid::dimension;
	typedef typename Grid::ctype ct;

	typedef Dune::FieldVector<ct, cdim>      GlobalCoordinate;
	typedef Dune::FieldVector<ct, mydim>  LocalCoordinate;


	typedef typename Grid::Traits::LeafIntersection    Intersection;

    static const unsigned int RETURN_SIZE = 1;
    typedef GlobalCoordinate               ResultValue;
    typedef std::vector<GlobalCoordinate>  ResultType;

	Intersection I_;

	NormalFunctor(const Intersection & I) : I_(I)  { }

	// Calculates the outer normal to the intersection times the integration element
	ResultType operator()(const LocalCoordinate & x) const
    {
    	return ResultType(1, I_.unitOuterNormal(x));
    }

	GlobalCoordinate zeroValue(unsigned int rezIndex) const { return GlobalCoordinate(0.0); }
};


// This method iterates over all elements of the grid
// Then over all intersections
// Then for each intersection that is a boundary segment it integrates the normal component by component
template<class GridType>
void Integrate (GridType & grid) {
	// Define standard grid properties
	const int cdim =  GridType::dimension;
	typedef typename GridType::ctype ct;
	typedef typename GridType::LeafGridView LeafGridView;
	typedef typename LeafGridView::template Codim<1>::Geometry FaceGeometry;
	typedef Dune::FieldVector<ct, cdim>      GlobalCoordinate;

	// Define Functor and Integrator
	typedef NormalFunctor<GridType, DIM2D>             Integrand2DVector;
	typedef Dune::QuadratureIntegrator<ct, DIM2D>  Integrator2DVector;
	typedef typename Integrator2DVector::template Traits<Integrand2DVector>::StatInfo  StatInfo;

	// Define Integrator parameters
	double RELATIVE_TOLERANCE = 1.0e-5;
	double ACCURACY_GOAL = 1.0e-15;
	const int NORM_TYPE = Dune::QUADRATURE_NORM_L2;

	// Initialize the integral result
	GlobalCoordinate  normalintegral(0.0);

	// Iterate over entities, then over intersections (faces) of the grid
	LeafGridView leafView = grid.leafGridView();
	const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

	for (auto&& elementThis : elements(leafView, Dune::Partitions::interiorBorder)) {
		std::cout << "-accessing element " << indexSet.index(elementThis) << std::endl;
		for (auto&& intersection : intersections(leafView, elementThis)) {
			if (!intersection.neighbor())
			{
			  Dune::GeometryType gt = intersection.type();
			  FaceGeometry geometry = intersection.geometry();

			  Integrand2DVector n(intersection);

			  StatInfo thisIntegralN = Integrator2DVector::template integrateRecursive<FaceGeometry, Integrand2DVector, NORM_TYPE>(geometry, n, RELATIVE_TOLERANCE, ACCURACY_GOAL);

			  std::stringstream logsstr;
			  logsstr << "---- from entity " << gt;
			  logsstr << " adding normal integral contribution " <<  thisIntegralN.second[0];
			  logsstr << ". Needed order " << thisIntegralN.first;
			  Dune::LoggingMessage::template write<Dune::CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, logsstr.str());

			  normalintegral += thisIntegralN.second[0];
			}
		}
	}


	// The actual integral is the sum over all processors
	// Unfortunately, DynamicVector can not be directly communicated since it is dynamic
	GlobalCoordinate rez = grid.comm().sum(normalintegral);
	if (grid.comm().rank() == 0)  { std::cout << "Normal integral amounted to " << rez << std::endl; }
}



int main (int argc , char **argv) {
	typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	const int grid_file_type = 1;  // createGrid procedure provides 6 different example grids numbered 0 to 5

	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;
	GridType * grid = createGrid<GridType>(mpihelper, grid_file_type);

	// Allow parallel output over all processes for the remainder of the program
	Dune::LoggingMessage::setPVerbose(true);

	// Traverse all entities of the grid and write information about each entity
	LoggingTimerDev::time("Integrating normal integral");
	Integrate(*grid);
	LoggingTimerDev::time("Integrating normal integral");

	// Use barrier to clean up the output a little bit. Note that this slows down the program
	grid->comm().barrier();

	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}

