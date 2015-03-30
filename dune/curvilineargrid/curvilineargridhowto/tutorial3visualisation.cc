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
#include <dune/curvilineargrid/curvilineargrid/factory.hh>

#include <dune/curvilineargrid/io/file/curvilineargmshreader.hh>




const bool isCached = true;


template <class GridType>
GridType * createGrid(Dune::MPIHelper & mpihelper, Dune::LoggingMessage & loggingmessage)
{
	// Obtain path to the mesh folder from a CMake constant
    const std::string CURVILINEARGRID_TEST_GRID_PATH = std::string(DUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH) + "curvilinear/";

    // A choice of 32 tetrahedra spheres with interpolation orders 1 to 5.
    const std::string GMSH_FILE_NAME[5] {"sphere32.msh", "sphere32ord2.msh", "sphere32ord3.msh", "sphere32ord4.msh", "sphere32ord5.msh"};

    // Choice of file name
    int interpOrder = 1;
    std::string filename = CURVILINEARGRID_TEST_GRID_PATH + GMSH_FILE_NAME[interpOrder - 1];

    // Additional constants
    bool insertBoundarySegment = true;  // If boundary segments will be inserted from GMSH. At the moment MUST BE true
    bool withGhostElements = true;      // to create Ghost elements
    bool writeReaderVTKFile = false;    // to write mesh to VTK during reading stage

    // Construct the grid factory
    Dune::CurvilinearGridFactory<GridType> factory(withGhostElements, mpihelper, loggingmessage);

    // Factory requires total vertex and element number in the mesh for faster performance
    int nVertexTotal;
    int nElementTotal;

    // Read the mesh into the factory using Curvilinear GMSH Reader
    Dune::CurvilinearGmshReader< GridType >::read(factory,
                                                  filename,
                                                  mpihelper,
                                                  loggingmessage,
                                                  writeReaderVTKFile,
                                                  insertBoundarySegment);

    // Create the grid
    return factory.createGrid();
}


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
	bool WRITE_VTK_EDGES,
	bool WRITE_VTK_TRIANGLES,
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
	  LocalIndexType  localIndex         = grid.leafIndexSet().index(*it);
	  GlobalIndexType globalIndex        = grid.template entityGlobalIndex<codim>(*it);
	  PhysicalTagType physicalTag        = grid.template entityPhysicalTag<codim>(*it);
	  StructuralType  structType         = grid.template entityStructuralType<codim>(*it);
	  InterpolatoryOrderType interpOrder = grid.template entityInterpolationOrder<codim>(*it);

	  // If we requested to output entities of this type, we will write them to VTK
	  if (typeset.find(structType) != typeset.end())
	  {
		std::vector<int>         tags  { physicalTag, structType, mpihelper.rank() };

    	vtkCurvWriter.template addCurvilinearElement<mydim>(
    			gt,
    			interpVertices,
    			tags,
    			interpOrder,
    			N_DISCRETIZATION_POINTS,
    			interpolate,
    			explode,
    			WRITE_VTK_EDGES,
    			WRITE_VTK_TRIANGLES);
	  }
  }

}


int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	const int dim = 3;
	const int dimworld = 3;
	typedef  double    ctype;
	typedef Dune::CurvilinearGrid<dim, dimworld, ctype, isCached> GridType;
	typedef typename GridType::GridStorageType         GridStorageType;
	typedef typename GridType::StructuralType          StructuralType;

	// Verbose constants
    bool verbose = true;                // to write logging output on master process
    bool processVerbose = true;         // to write logging output on all processes
    Dune::LoggingMessage loggingmessage(verbose, processVerbose, mpihelper);

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, loggingmessage);


	// Define the structural types used by the VTK writer
	const int InternalType        = GridStorageType::PartitionType::Internal;
	const int GhostType           = GridStorageType::PartitionType::Ghost;
	const int DomainBoundaryType  = GridStorageType::PartitionType::DomainBoundary;
	const int ProcessBoundaryType = GridStorageType::PartitionType::ProcessBoundary;


	// Additional VTK writer constants
	// Note that these can be specified uniquely for each element added to the writer if necessary
	int N_DISCRETIZATION_POINTS = 5;    // Number of linear points to subdivide a curvilinear line (min=2 - linear)
	bool interpolate = true;            // If interpolate=false, vtk writer just uses the interpolatory vertices as discretization vertices and ignores above constant
	bool explode = true;                // Artificially increase distance between all entities, by shrinking all entities a bit wrt their center of mass
	bool WRITE_VTK_EDGES;               // If edges should be used to discretize this element
	bool WRITE_VTK_TRIANGLES;           // If triangles should be used to discretize this element


    // Construct the VTK writer
	Dune::CurvilinearVTKWriter<GridType> vtkCurvWriter(mpihelper, loggingmessage);


	// write elements
	// *************************************************
	WRITE_VTK_EDGES     = false;       // If edges should be used to discretize this element
	WRITE_VTK_TRIANGLES = true;        // If triangles should be used to discretize this element
	std::set<StructuralType>  writeElements;
	writeElements.insert(InternalType);
	writeElements.insert(GhostType);
	writeVTKentities<0, GridType>(
			*grid,
			mpihelper,
			vtkCurvWriter,
			N_DISCRETIZATION_POINTS,
			interpolate,
			explode,
			WRITE_VTK_EDGES,
			WRITE_VTK_TRIANGLES,
			writeElements);






	// Writing the PVTU and VTU files
	vtkCurvWriter.writeParallelVTU("./curvreader_output");

    // Delete the grid
    delete grid;


    return 0;
}
