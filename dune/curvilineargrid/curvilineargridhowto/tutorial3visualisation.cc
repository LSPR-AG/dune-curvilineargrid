// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>




const bool isCached = true;


// example for a generic algorithm that traverses
// the entities of a given mesh, extracts the relevant information
// and inserts it to the VTK writer
template<int codim, class GridType>
void writeVTKentities (
	GridType& grid,
	Dune::MPIHelper & mpihelper,
	Dune::CurvilinearVTKWriter<GridType> & vtkCurvWriter,
	int N_DISCRETIZATION_POINTS,
	bool interpolate,
	bool explode,
	std::vector<bool> writeCodim,
	std::set<typename GridType::StructuralType > typeset)
{
  const int dim =  GridType::dimension;
  const int mydim = dim - codim;
  typedef typename GridType::ctype ct;
  typedef typename GridType::LeafGridView LeafGridView;

  typedef typename GridType::LocalIndexType          LocalIndexType;
  typedef typename GridType::GlobalIndexType         GlobalIndexType;
  typedef typename GridType::PhysicalTagType         PhysicalTagType;
  typedef typename GridType::StructuralType          StructuralType;
  typedef typename GridType::InterpolatoryOrderType  InterpolatoryOrderType;

  // get the instance of the LeafGridView
  LeafGridView leafView = grid.leafGridView();

  typedef typename LeafGridView::template Codim<codim>::Iterator EntityLeafIterator;
  //typedef typename LeafGridView::template Codim<codim>::Geometry EntityGeometry;
  typedef typename GridType::template Codim<codim>::EntityGeometryMappingImpl  BaseGeometry;
  typedef typename BaseGeometry::GlobalCoordinate                              GlobalCoordinate;

  // Iterate over entities of this codimension
  EntityLeafIterator ibegin = leafView.template begin<codim>();
  EntityLeafIterator iend   = leafView.template end<codim>();

  for (EntityLeafIterator it = ibegin; it != iend; ++it)
  {
	  // Constructing a geometry is quite expensive, do it only once
	  //EntityGeometry geom = it->geometry();
	  BaseGeometry geom = grid.template entityBaseGeometry<codim>(*it);
	  std::vector<GlobalCoordinate>  interpVertices = geom.vertexSet();

	  Dune::GeometryType gt              = it->type();
	  Dune::PartitionType ptype          = it->partitionType();
	  LocalIndexType  localIndex         = grid.leafIndexSet().index(*it);
	  GlobalIndexType globalIndex        = grid.template entityGlobalIndex<codim>(*it);
	  PhysicalTagType physicalTag        = grid.template entityPhysicalTag<codim>(*it);
	  InterpolatoryOrderType interpOrder = grid.template entityInterpolationOrder<codim>(*it);

	  // If we requested to output entities of this type, we will write them to VTK
	  if (typeset.find(ptype) != typeset.end())
	  {
		std::vector<int>         tags  { physicalTag, ptype, mpihelper.rank() };

    	vtkCurvWriter.template addCurvilinearElement<mydim>(
    			gt,
    			interpVertices,
    			tags,
    			interpOrder,
    			N_DISCRETIZATION_POINTS,
    			interpolate,
    			explode,
    			writeCodim);
	  }
  }

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

	// Define curvilinear grid
    typedef Dune::CurvilinearGrid<ctype, dim, isCached, LoggingMessageDev> GridType;
	typedef typename GridType::GridStorageType         GridStorageType;
	typedef typename GridType::StructuralType          StructuralType;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper);


	// Additional VTK writer constants
	// Note that these can be specified uniquely for each element added to the writer if necessary
	int N_DISCRETIZATION_POINTS = 5;    // Number of linear points to subdivide a curvilinear line (min=2 - linear)
	bool interpolate = true;            // If interpolate=false, vtk writer just uses the interpolatory vertices as discretization vertices and ignores above constant
	bool explode = true;                // Artificially increase distance between all entities, by shrinking all entities a bit wrt their center of mass
	std::vector<bool> writeCodim {true, false, false, false};  // For now, only use linear elements to discretize inserted entities


    // Construct the VTK writer
	Dune::CurvilinearVTKWriter<GridType> vtkCurvWriter(mpihelper);


	// write elements
	// *************************************************
	std::set<StructuralType>  writeElements;
	writeElements.insert(Dune::PartitionType::InteriorEntity);
	writeElements.insert(Dune::PartitionType::GhostEntity);
	writeVTKentities<0, GridType>(
			*grid,
			mpihelper,
			vtkCurvWriter,
			N_DISCRETIZATION_POINTS,
			interpolate,
			explode,
			writeCodim,
			writeElements);






	// Writing the PVTU and VTU files
	vtkCurvWriter.writeParallelVTU("./curvreader_output");

    // Delete the grid
    delete grid;


    return 0;
}
