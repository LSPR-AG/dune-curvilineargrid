// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/curvilineargrid/curvilineargrid/grid.hh>
#include <dune/curvilineargrid/io/file/curvilinearvtkgridwriter.hh>
#include <dune/curvilineargrid/curvilineargridhowto/creategrid.hh>




const bool isCached = true;



// Generates Sinusoid as function of local coordinates
// Note: This method is very naive, as it does not even consider the global orientation of the element
// Normally you can find out the global orientation by considering the global coordinates and/or global ID's of the associated geometry
template <typename Grid>
class VTKFunctionLocalSinusoid	:	public Dune::VTKFunction<typename Grid::LeafGridView>
{
protected:
	/***************************************************************************************************************/
	/* extract information from the Grid and                                                                       */
	/***************************************************************************************************************/
    static const int dimension = Grid::dimensionworld;
	typedef typename Grid::ctype CoordinateType;

	typedef typename Grid::template Codim<0>::Entity EntityElement;

    typedef Dune::FieldVector< CoordinateType, dimension >  LocalCoordinate;
    typedef Dune::FieldVector< CoordinateType, dimension >  GlobalCoordinate;

public:
    VTKFunctionLocalSinusoid()  {  }

    // Number of components of the field (3 in 3d)
	virtual int ncomps() const  { return(dimension); }

	// Component of the field
	virtual double evaluate(int comp, const EntityElement & element, const LocalCoordinate &xi) const
	{
		return sin(5 * xi[comp]);
	}

	/**********************************************************************/
	virtual std::string name() const { return std::string("Local-Sinusoid"); }
};


/** \brief Describes a sinusoid function in the global coordinates of the world
 *  The only difference to local sinusoid method is having to map from local to global coordinates
 *    */
template <typename Grid>
class VTKFunctionGlobalSinusoid	:	public Dune::VTKFunction<typename Grid::LeafGridView>
{
protected:
	/***************************************************************************************************************/
	/* extract information from the Grid and                                                                       */
	/***************************************************************************************************************/
    static const int dimension = Grid::dimensionworld;
	typedef typename Grid::ctype CoordinateType;

	typedef typename Grid::template Codim<0>::Entity EntityElement;

    typedef Dune::FieldVector< CoordinateType, dimension >  LocalCoordinate;
    typedef Dune::FieldVector< CoordinateType, dimension >  GlobalCoordinate;

public:
    VTKFunctionGlobalSinusoid(const Grid & grid) : grid_(grid)
	{

	}

    // Number of components of the field (3 in 3d)
	virtual int ncomps() const  { return(dimension); }

	// Component of the field
	virtual double evaluate(int comp, const EntityElement & element, const LocalCoordinate &xi) const
	{
		return sin(5 * element.geometry().global(xi)[comp]);
	}

	/**********************************************************************/
	virtual std::string name() const { return std::string("Global-Sinusoid"); }

private:
	const Grid & grid_;
};









int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);


	// Define curvilinear grid
	const int dim = 3;
	typedef  double    ctype;
	const int grid_file_type = 4;  // createGrid procedure provides 6 different example grids numbered 0 to 5

	typedef Dune::CurvilinearGrid<ctype, dim, isCached, Dune::LoggingMessage> GridType;
	typedef typename GridType::GridStorageType         GridStorageType;
	typedef typename GridType::StructuralType          StructuralType;
	typedef typename Dune::VTKFunction<typename GridType::LeafGridView>   DuneVTKFunction;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, grid_file_type);


    // Construct the VTK writer
	Dune::CurvilinearVTKGridWriter<GridType> writer(*grid);

	// Set fixed virtual refinement to better resolve the fine sinusoidal fields
	// This action is not recommended for large grids because it will severely slow down the writing procedure.
	// If the function polynomial order is not higher than the curvature order of the element, then the writer
	// will automatically adapt the correct virtual refinement order.
	const int userDefinedVirtualRefinement = 15;
	writer.useFixedVirtualRefinement(userDefinedVirtualRefinement);

	// We will attempt to attach two fields to the grid and plot them
	std::vector<DuneVTKFunction *> vtkFuncSet_; /** \brief Stores all functors that will be plotted over the mesh during diagnostics phase **/
    vtkFuncSet_.push_back(new VTKFunctionLocalSinusoid<GridType>());        // Sinusoid in local coordinates
    vtkFuncSet_.push_back(new VTKFunctionGlobalSinusoid<GridType>(*grid));  // Sinusoid in global coordinates

    writer.write("./", "tutorial-3-output", vtkFuncSet_);


    // Delete the grid
    delete grid;


    return 0;
}
