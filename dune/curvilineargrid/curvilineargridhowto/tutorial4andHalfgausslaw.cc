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
#include <dune/curvilineargrid/curvilineargrid/factory.hh>

#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>




const bool isCached = true;


typedef Dune::FieldVector<double, 3>  GlobalCoordinate;


// This Functor creates a field of a single charge at the origin of global coordinates
// Then calculates the projection of this field on the local surface times the integration element
// This functor is meant to be directly integrated over reference element to calculate the contribution
// to Gauss law from a boundary surface
template<class Grid, int mydim>
struct GaussFunctor
{
	typedef Dune::FieldVector<double, mydim>  LocalCoordinate;

	typedef typename Grid::Traits::LeafIntersection               Intersection;
	typedef typename Grid::Traits::template Codim<1>::Geometry    FaceGeometry;

	Intersection I_;
	FaceGeometry G_;

	// Instantiate with Intersection of interest. Also store geometry,
	// so that need not recalculate it at every sample point
    GaussFunctor(const Intersection & I)
	  : I_(I),
	    G_(I_.geometry())
	{

	}

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
    double operator()(const LocalCoordinate & x) const
    {
    	GlobalCoordinate integrnormal = I_.integrationOuterNormal(x);
    	GlobalCoordinate global = G_.global(x);
    	GlobalCoordinate field = ChargeField(global);

    	double rez = 0;
    	for (int i = 0; i < 3; i++) { rez += integrnormal[i] * field[i]; }
    	return rez;
    }
};



template <class GridType>
GridType * createGrid(Dune::MPIHelper & mpihelper)
{
	// Obtain path to the mesh folder from a CMake constant
    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

    // A choice of 32 tetrahedra spheres with interpolation orders 1 to 5.
    const std::string GMSH_FILE_NAME[5] {"sphere32.msh", "sphere32ord2.msh", "sphere32ord3.msh", "sphere32ord4.msh", "sphere32ord5.msh"};

    // Choice of file name
    int interpOrder = 3;
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME[interpOrder - 1];

    // Additional constants
    bool insertBoundarySegment = true;  // If boundary segments will be inserted from GMSH. At the moment MUST BE true
    bool withGhostElements = true;      // to create Ghost elements
    bool verbose = true;                // to write logging output on master process
    bool processVerbose = true;         // to write logging output on all processes
    bool writeReaderVTKFile = false;    // to write mesh to VTK during reading stage

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, verbose, processVerbose, mpihelper);

    // Factory requires total vertex and element number in the mesh for faster performance
    int nVertexTotal;
    int nElementTotal;

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory,
                                                  filename,
                                                  mpihelper,
                                                  verbose,
                                                  processVerbose,
                                                  writeReaderVTKFile,
                                                  insertBoundarySegment);

    // Create the grid
    return factory.createGrid();
}



// This method iterates over all elements of the grid
// Then over all intersections
// Then for each intersection that is a boundary segment it integrates the Surface projection of EM field of a single charge
// The output of this function should be 4pi for whatever closed geometry that includes the origin
template<class GridType>
void gaussIntegral (GridType& grid)
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


  double gaussintegral = 0.0;
  double rel_tol = 1.0e-5;

  // Iterate over entities of this codimension
  EntityLeafIterator ibegin = leafView.template begin<0>();
  EntityLeafIterator iend   = leafView.template end<0>();

  for (EntityLeafIterator it = ibegin; it != iend; ++it)
  {
	  const Entity &entity = *it;

	  std::cout << "-accessing entity " << indexSet.index(entity) << std::endl;

	  const IntersectionIterator nend = leafView.iend(entity);
	  for( IntersectionIterator nit = leafView.ibegin(entity); nit != nend; ++nit )
	  {
		  const Intersection &intersection = *nit;

		  if (!intersection.neighbor())
		  {
			  Dune::GeometryType gt = intersection.type();
			  GaussFunctor<GridType, 2> f(intersection);
			  Dune::QuadratureIntegrator<double, 2> qInt;
			  double thisIntegral = qInt.integrateRecursive(gt, f, rel_tol).second;

			  std::cout << "-- adding contribution from " << gt << "  " << thisIntegral << std::endl;

			  gaussintegral += thisIntegral;
		  }
	  }
  }

  std::cout << "Gauss integral amounted to " << gaussintegral << std::endl;
}



int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	const int dimworld = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<dim, dimworld, ctype, isCached> GridType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper);

	// Traverse all entities of the grid and write information about each entity
	gaussIntegral(*grid);

    // Delete the grid
    delete grid;


    return 0;
}

