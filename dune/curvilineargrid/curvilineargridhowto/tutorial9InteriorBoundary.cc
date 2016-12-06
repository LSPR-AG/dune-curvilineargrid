// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/curvilineargeometry/integration/quadratureintegrator.hh>

#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>




const bool isCached = true;
typedef Dune::FieldVector<double, 3>  GlobalCoordinate;


// This Functor creates a field of a single charge at the origin of global coordinates
// Then calculates the projection of this field on the local surface times the integration element
// This functor is meant to be directly integrated over reference element to calculate the contribution
// to Gauss law from a boundary surface
template<class Grid, int mydim>
struct GaussFunctor
{
	typedef Dune::FieldVector<double, mydim>       LocalCoordinate;

    static const unsigned int RETURN_SIZE = 1;
    typedef double                    ResultValue;
    typedef std::vector<ResultValue>  ResultType;

	typedef typename Grid::Traits::LeafIntersection               Intersection;

	const Intersection & I_;
	GlobalCoordinate x0_;     // Displace the charge a bit on x-axis to make the integrals a bit more non-trivial

	// Instantiate with Intersection of interest.
    GaussFunctor(const Intersection & I, const GlobalCoordinate & x0) : I_(I), x0_(x0)  {
    }

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
template<class GridType>
std::vector<double> Integrate (GridType& grid, const std::vector<GlobalCoordinate> & chargePos, int tagVol, int tagSurf)
{
  const int dim =  GridType::dimension;
  typedef typename GridType::ctype ct;
  typedef typename GridType::LeafGridView LeafGridView;
  typedef typename LeafGridView::IntersectionIterator IntersectionIterator;
  typedef typename LeafGridView::template Codim< 0 >::Entity Element;
  typedef typename LeafGridView::template Codim< 1 >::Entity Face;
  typedef typename IntersectionIterator :: Intersection Intersection;

  // get the instance of the LeafGridView
  LeafGridView leafView = grid.leafGridView();
  const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

  typedef typename LeafGridView::template Codim<0>::Iterator EntityLeafIterator;
  typedef typename LeafGridView::template Codim<1>::Geometry FaceGeometry;

  typedef typename GridType::PhysicalTagType         PhysicalTagType;

  typedef GaussFunctor<GridType, 2>              Integrand2D;
  typedef Dune::QuadratureIntegrator<double, 2>  Integrator2DScalar;
  typedef typename Integrator2DScalar::template Traits<Integrand2D>::StatInfo  StatInfo;


  std::vector<double> gaussintegral(chargePos.size(), 0.0);
  double RELATIVE_TOLERANCE = 1.0e-5;
  double ACCURACY_GOAL = 1.0e-15;
  const int NORM_TYPE = Dune::QUADRATURE_NORM_L2;

  for (auto&& elem : elements(leafView, Dune::Partitions::interiorBorder)) {

	  // Only count faces neighboring the desired subdomain.
	  // The subdomain is necessary to determine the direction of the unit outer normal
	  PhysicalTagType physicalTagElem = grid.template entityPhysicalTag<0>(elem);
	  if (physicalTagElem == tagVol) {
		 //std::cout << "-accessing entity " << indexSet.index(elem) << std::endl;

		  for(auto && intr : intersections(leafView, elem) )
		  {
			  // Only count faces of a boundary of the selected type, as the subdomain can have multiple boundaries
			  Face face = elem.template subEntity<1>(intr.indexInInside());
			  PhysicalTagType physicalTagFace = grid.template entityPhysicalTag<1>(face);
			  if (physicalTagFace == tagSurf) {
				  FaceGeometry geometry = intr.geometry();
				  Dune::GeometryType gt = intr.type();

				  for (int iCharge = 0; iCharge < chargePos.size(); iCharge++) {
					  Integrand2D g(intr, chargePos[iCharge]);

					  StatInfo thisIntegralG = Integrator2DScalar::template integrateRecursive<FaceGeometry, Integrand2D, NORM_TYPE>(geometry, g, RELATIVE_TOLERANCE, ACCURACY_GOAL);
					  std::cout << "---- adding gauss contribution from " << gt << "  " << thisIntegralG.second[0] << ". Needed order " << thisIntegralG.first << std::endl;

					  gaussintegral[iCharge] += thisIntegralG.second[0];
				  }
			  }
		  }
	  }
  }

  // The actual integral is the sum over all processors
  std::vector<double> rez(chargePos.size());
  for (int iCharge = 0; iCharge < chargePos.size(); iCharge++) {
	  rez[iCharge] = grid.comm().sum(gaussintegral[iCharge]);
  }

  return rez;
}



int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

	// NOTE: This test only makes sense for geometries with several interior boundaries
	const int grid_file_type = 7;  // createGrid procedure provides 8 different example grids numbered 0 to 7

	static const int INNER_VOLUME_TAG = 501;
	static const int OUTER_VOLUME_TAG = 503;
	static const int INNER_SURFACE_TAG = 101;
	static const int OUTER_SURFACE_TAG = 102;

	typedef Dune::CurvilinearGrid<ctype, dim, isCached> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, grid_file_type);

	std::vector<double> xCoord{0.0,
		-0.01, -0.02, -0.03, -0.04, -0.06, -0.4, -0.6, -0.9, -1.1, -2.0,
		0.01, 0.02, 0.03, 0.04, 0.06, 0.4, 0.6, 0.9, 1.1, 100.0
	};

	// Expected integral result
	double A = 4.0 * M_PI;
	std::vector<double> rez1 {A,
		A, A, A, A, 0, 0, 0, 0, 0, 0,
		A, A, A, A, 0, 0, 0, 0, 0, 0
	};

	std::vector<double> rez2 {A,
		A, A, A, A, A, A, A, A, 0, 0,
		A, A, A, A, A, A, A, A, 0, 0
	};


	std::vector<GlobalCoordinate> chargePos(xCoord.size(), GlobalCoordinate(0.0));
	for (int i = 0; i < xCoord.size(); i++)  { chargePos[i][0] = xCoord[i]; }

	// Perform the integration
	std::vector<double> I1 = Integrate(*grid, chargePos, INNER_VOLUME_TAG, INNER_SURFACE_TAG);
	std::vector<double> I2 = Integrate(*grid, chargePos, OUTER_VOLUME_TAG, OUTER_SURFACE_TAG);

	// Report result on the master
	if (grid->comm().rank() == Dune::CurvGrid::MPI_MASTER_RANK) {
		for (int i = 0; i < xCoord.size(); i++) {
			std::cout
				<< "For charge at pos (" << chargePos[i]
				<<") the integration results: F_in=" << I1[i] <<", err=" << I1[i] - rez1[i]
				<<", F_out=" << I2[i] <<", err=" << I2[i] - rez2[i] << std::endl;
		}
	}

	typedef Dune::LoggingTimer<Dune::LoggingMessage>                 LoggingTimerDev;
	LoggingTimerDev::reportParallel();

    // Delete the grid
    delete grid;


    return 0;
}

