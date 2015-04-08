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


// Calculates the outer normal to the intersection times the integration element
template<class Grid, int mydim>
struct NormalFunctor
{
	typedef Dune::FieldVector<double, mydim>       LocalCoordinate;
	typedef typename Grid::Traits::LeafIntersection               Intersection;

	Intersection I_;

	NormalFunctor(const Intersection & I) : I_(I)  { }

	// Calculates the outer normal to the intersection times the integration element
	GlobalCoordinate operator()(const LocalCoordinate & x) const
    {
    	return I_.integrationOuterNormal(x);
    }
};


// This method iterates over all elements of the grid
// Then over all intersections
// Then for each intersection that is a boundary segment it integrates the normal component by component
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

  typedef Dune::QuadratureIntegrator<double, 2, 1>  Integrator2DScalar;
  typedef Dune::QuadratureIntegrator<double, 2, 3>  Integrator2DVector;
  typedef Dune::FieldVector<double, 1>              ResultType;

  GlobalCoordinate  normalintegral(0.0);
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

		  if (!intersection.neighbor())
		  {
			  Dune::GeometryType gt = intersection.type();

			  NormalFunctor<GridType, 2> n(intersection);

			  Integrator2DVector::StatInfo thisIntegralN = Integrator2DVector::integrateRecursive(gt, n, rel_tol);
			  std::cout << "---- adding normal contribution from " << gt << "  " << thisIntegralN.second << ". Needed order " << thisIntegralN.first << std::endl;

			  normalintegral += thisIntegralN.second;
		  }
	  }
  }
  std::cout << "Normal integral amounted to " << normalintegral << std::endl;
}



int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;

    // Instantiation of the logging message and loggingtimer
    typedef Dune::LoggingMessage<Dune::LoggingMessageHelper::Phase::DEVELOPMENT_PHASE>   LoggingMessageDev;
    typedef Dune::LoggingTimer<LoggingMessageDev>                                        LoggingTimerDev;
    LoggingMessageDev::getInstance().init(mpihelper, true, true);
    LoggingTimerDev::getInstance().init(false);

	typedef Dune::CurvilinearGrid<ctype, dim, isCached, LoggingMessageDev> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper);

	// Traverse all entities of the grid and write information about each entity
	Integrate(*grid);

    // Delete the grid
    delete grid;


    return 0;
}

