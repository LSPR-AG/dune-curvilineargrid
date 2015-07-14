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

	// Instantiate with Intersection of interest.
    GaussFunctor(const Intersection & I) : I_(I)  { }

    // Calculates EM field vector of a single charge assuming the charge is at the origin of global coordinates
    // Also assuming dielectric permittivity epsilon = 1
    static GlobalCoordinate ChargeField(const GlobalCoordinate & x)
    {
    	GlobalCoordinate rez = x;
    	rez /= pow(x.two_norm2(), 1.5);
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
void Integrate (GridType& grid)
{
  const int dim =  GridType::dimension;
  typedef typename GridType::ctype ct;
  typedef typename GridType::LeafGridView LeafGridView;
  typedef typename LeafGridView::IntersectionIterator IntersectionIterator;
  typedef typename LeafGridView::template Codim< 0 >::Entity Entity;
  typedef typename IntersectionIterator :: Intersection Intersection;




  // get the instance of the LeafGridView
  LeafGridView leafView = grid.leafGridView();
  const typename GridType::LeafIndexSet & indexSet = grid.leafIndexSet();

  typedef typename LeafGridView::template Codim<0>::Iterator EntityLeafIterator;
  typedef typename LeafGridView::template Codim<1>::Geometry FaceGeometry;


  typedef GaussFunctor<GridType, 2>              Integrand2D;
  typedef Dune::QuadratureIntegrator<double, 2>  Integrator2DScalar;
  typedef typename Integrator2DScalar::template Traits<Integrand2D>::StatInfo  StatInfo;


  double gaussintegral = 0.0;
  double rel_tol = 1.0e-5;

  // Iterate over entities of this codimension
  EntityLeafIterator ibegin = leafView.template begin<0>();
  EntityLeafIterator iend   = leafView.template end<0>();

  for (EntityLeafIterator it = ibegin; it != iend; ++it)
  {
	  const Entity &entity = *it;

	  std::cout << "-accessing entity " << indexSet.index(entity) << std::endl;

	  const IntersectionIterator nbegin = leafView.ibegin(entity);
	  const IntersectionIterator nend = leafView.iend(entity);

	  for( IntersectionIterator nit = nbegin; nit != nend; ++nit )
	  {
		  const Intersection &intersection = *nit;
		  FaceGeometry geometry = intersection.geometry();


		  if (!intersection.neighbor())
		  {
			  Dune::GeometryType gt = intersection.type();

			  Integrand2D g(intersection);

			  StatInfo thisIntegralG = Integrator2DScalar::integrateRecursive(geometry, g, rel_tol);
			  std::cout << "---- adding gauss contribution from " << gt << "  " << thisIntegralG.second[0] << ". Needed order " << thisIntegralG.first << std::endl;

			  gaussintegral += thisIntegralG.second[0];
		  }
	  }
  }

  std::cout << "Gauss integral amounted to " << gaussintegral << std::endl;
}



int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

	typedef Dune::CurvilinearGrid<ctype, dim, isCached, Dune::LoggingMessage> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper);

	// Traverse all entities of the grid and write information about each entity
	Integrate(*grid);

    // Delete the grid
    delete grid;


    return 0;
}

