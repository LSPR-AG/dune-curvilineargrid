/********************************************
 * Dune-CurvilinearGrid
 * Tutorial 3: Visualisation
 *
 * Author: Aleksejs Fomins
 *
 * Description: This tutorial samples 3D sine functions over the grid elements (3D)
 * and boundary segments (2D), using local and global coordinates. It then writes
 * the grid and the associated fields to the VTU/PVTU file set.
 ********************************************/

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



using namespace Dune;
using namespace Dune::CurvGrid;

const int DIM0D = 0;
const int DIM1D = 1;
const int DIM2D = 2;
const int DIM3D = 3;
const bool isCached = true;



template <typename Grid>
class VTKFunctionElementIndex	:	public VTKScalarFunction <Grid, DIM3D>
{
	typedef VTKScalarFunction <Grid, DIM3D>  Base;

protected:
	/***************************************************************************************************************/
	/* extract information from the Grid and                                                                       */
	/***************************************************************************************************************/
    static const int dimension   = Grid::dimensionworld;
    static const int mydimension = dimension;
	typedef typename Grid::LeafIndexSet     LeafIndexSet;

	typedef typename Base::Entity   Entity;
	typedef typename Base::Domain   Domain;
	typedef typename Base::Range    Range;

public:


	VTKFunctionElementIndex(const Grid & grid) : indexset_(grid.leafIndexSet()) {  }

    virtual void init (const Entity & entity) { thisIndex_ = indexset_.index(entity); }

	// Component of the field
    virtual Range evaluate(const Domain & x) const  { return thisIndex_; }

    virtual std::string name() const  { return "localIndex3D"; }

private:
    const LeafIndexSet & indexset_;
    Range thisIndex_;
};


template <typename Grid>
class VTKFunctionBoundarySegmentIndex	:	public VTKScalarFunction <Grid, DIM2D>
{
	typedef VTKScalarFunction <Grid, DIM2D>  Base;

protected:
	/***************************************************************************************************************/
	/* extract information from the Grid and                                                                       */
	/***************************************************************************************************************/
    static const int dimension   = Grid::dimensionworld;
    static const int mydimension = dimension;
	typedef typename Grid::LeafIndexSet     LeafIndexSet;
	typedef typename Grid::LocalIndexType   LocalIndexType;

	typedef typename Base::Entity   Entity;
	typedef typename Base::Domain   Domain;
	typedef typename Base::Range    Range;

public:


	VTKFunctionBoundarySegmentIndex(const Grid & grid) : grid_(grid) {  }

	// Note: In this case entity is an Intersection
    virtual void init (const Entity & entity) {
		//LocalIndexType faceIndex = grid_.leafIndexSet().index(entity);
		//thisIndex_ = grid_.gridbase().boundarySegmentIndex(faceIndex);   // Convert from face local index to boundary segment local index
    	thisIndex_ = entity.boundarySegmentIndex();
    }

	// Component of the field
    virtual Range evaluate(const Domain & x) const  { return thisIndex_; }

    virtual std::string name() const  { return "localIndex2D"; }

private:
    const Grid & grid_;
    Range thisIndex_;
};



// Generates Sinusoid as function of local coordinates
// Note: This method is very naive, as it does not even consider the global orientation of the element
// Normally you can find out the global orientation by considering the global coordinates and/or global ID's of the associated geometry
template <typename Grid, int mydim>
class VTKFunctionLocalSinusoid	:	public VTKVectorFunction <Grid, mydim>
{
	typedef VTKVectorFunction <Grid, mydim>  Base;

protected:
	/***************************************************************************************************************/
	/* extract information from the Grid and                                                                       */
	/***************************************************************************************************************/
    static const int dimension   = Grid::dimensionworld;
    static const int mydimension = mydim;
	typedef typename Base::CoordinateType   CoordinateType;

	typedef typename Base::Entity   Entity;
	typedef typename Base::Domain   LocalCoordinate;
	typedef typename Base::Range    GlobalCoordinate;

public:


    VTKFunctionLocalSinusoid()  {  }

    virtual void init (const Entity & entity) { }

	// Component of the field
    virtual GlobalCoordinate evaluate(const LocalCoordinate & x) const
	{
    	GlobalCoordinate rez;
    	for (unsigned int i = 0; i < mydimension; i++)  { rez[i] = sin(5 * x[i]); }
		return rez;
	}

    virtual std::string name() const  { return "localSinusoid" + std::to_string(mydim) + "D"; }
};


/** \brief Describes a sinusoid function in the global coordinates of the world
 *  The only difference to local sinusoid method is having to map from local to global coordinates
 *    */
template <typename Grid, int mydim>
class VTKFunctionGlobalSinusoid	:	public VTKVectorFunction <Grid, mydim>
{
	typedef VTKVectorFunction <Grid, mydim>  Base;

protected:
	/***************************************************************************************************************/
	/* extract information from the Grid and                                                                       */
	/***************************************************************************************************************/
    static const int dimension   = Grid::dimensionworld;
    static const int mydimension = mydim;
	typedef typename Base::CoordinateType   CoordinateType;

	typedef typename Base::Entity   Entity;
	typedef typename Base::Domain   LocalCoordinate;
	typedef typename Base::Range    GlobalCoordinate;

public:

	VTKFunctionGlobalSinusoid() : entity_(nullptr)  {  }

    virtual void init (const Entity & entity)  { entity_ = &entity; }

	// Component of the field
    virtual GlobalCoordinate evaluate(const LocalCoordinate & x) const
	{
    	GlobalCoordinate rez;
    	for (unsigned int i = 0; i < mydimension; i++)  { rez[i] = sin(5 * entity_->geometry().global(x)[i]); }
		return rez;
	}

    virtual std::string name() const  { return "globalSinusoid" + std::to_string(mydim) + "D"; }

private:

    const Entity * entity_;
};




int main (int argc , char **argv) {
	static Dune::MPIHelper & mpihelper = Dune::MPIHelper::instance(argc, argv);

	// Define curvilinear grid
	typedef  double    ctype;

	typedef Dune::CurvilinearGrid<ctype, DIM3D, isCached>   GridType;
	typedef typename GridType::GridStorageType              GridStorageType;
	typedef typename GridType::StructuralType               StructuralType;
	typedef CurvilinearVTKGridWriter<GridType>        GridWriter;
	typedef typename GridWriter::VTKScalarFunction2D        BaseVTKScalarFunction2D;
	typedef typename GridWriter::VTKScalarFunction3D        BaseVTKScalarFunction3D;
	typedef typename GridWriter::VTKVectorFunction2D        BaseVTKVectorFunction2D;
	typedef typename GridWriter::VTKVectorFunction3D        BaseVTKVectorFunction3D;

	typedef VTKFunctionLocalSinusoid<GridType, DIM2D>       VTKFunctionLocalSinusoidFace;
	typedef VTKFunctionLocalSinusoid<GridType, DIM3D>       VTKFunctionLocalSinusoidElem;
	typedef VTKFunctionGlobalSinusoid<GridType, DIM2D>      VTKFunctionGlobalSinusoidFace;
	typedef VTKFunctionGlobalSinusoid<GridType, DIM3D>      VTKFunctionGlobalSinusoidElem;

	// Create Grid
	GridType * grid = createGrid<GridType>(mpihelper, argc , argv);

    // Construct the VTK writer
	GridWriter writer(*grid);

	/***************************************************************
	 * In the following block we set up the parameters of the writer
	 * Note that this is entirely optional, as all parameters have default values
	 ***************************************************************/

	// Set fixed virtual refinement to better resolve the fine sinusoidal fields
	// This action is not recommended for large grids because it will severely slow down the writing procedure.
	// If the function polynomial order is not higher than the curvature order of the element, then the writer
	// will automatically adapt the correct virtual refinement order.
	const int userDefinedVirtualRefinement = 7;
	writer.useFixedVirtualRefinement(userDefinedVirtualRefinement);
	writer.writeProcessBoundary(true);
	writer.writeDomainBoundary(true);
	writer.writeInteriorBoundary(true);
	writer.writePeriodicBoundary(grid->withPeriodic());  // Note: Only write periodic fields if there are at all any periodic elements
	writer.writePeriodicBind(grid->withPeriodic());  // Note: Only write periodic fields if there are at all any periodic elements
	writer.writeGhost(true);
	writer.writePatience(true);
	writer.writeInterpolate(true);
	writer.writeExplode(false);

	/***************************************************************
	 * We will attempt to attach two fields to the grid and plot them
	 ***************************************************************/

	/** \brief Stores all functors that will be plotted over the mesh during diagnostics phase **/
	std::vector<BaseVTKScalarFunction2D *> vtkFuncScalarSet2D_;
	std::vector<BaseVTKScalarFunction3D *> vtkFuncScalarSet3D_;
	std::vector<BaseVTKVectorFunction2D *> vtkFuncVectorSet2D_;
	std::vector<BaseVTKVectorFunction3D *> vtkFuncVectorSet3D_;


	vtkFuncScalarSet2D_.push_back(new VTKFunctionBoundarySegmentIndex<GridType>(*grid));   // Local index of boundary segments, sampled over faces
	vtkFuncScalarSet3D_.push_back(new VTKFunctionElementIndex<GridType>(*grid));           // Local index of elements, sampled over elements

	vtkFuncVectorSet2D_.push_back(new VTKFunctionLocalSinusoidFace());   // Sinusoid in local  coordinates, sampled over faces
	vtkFuncVectorSet2D_.push_back(new VTKFunctionGlobalSinusoidFace());  // Sinusoid in global coordinates, sampled over faces
	vtkFuncVectorSet3D_.push_back(new VTKFunctionLocalSinusoidElem());   // Sinusoid in local  coordinates, sampled over elements
	vtkFuncVectorSet3D_.push_back(new VTKFunctionGlobalSinusoidElem());  // Sinusoid in global coordinates, sampled over elements

    writer.addFieldSet(vtkFuncScalarSet2D_);
    writer.addFieldSet(vtkFuncScalarSet3D_);
    writer.addFieldSet(vtkFuncVectorSet2D_);
    writer.addFieldSet(vtkFuncVectorSet3D_);
    writer.write("./", "tutorial-3-output");

    typedef LoggingTimer<LoggingMessage>                 LoggingTimerDev;
    LoggingTimerDev::reportParallel();


    // Delete the grid
    delete grid;


    return 0;
}
