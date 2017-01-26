#ifndef CURVILINEAR_GRID_BASE_PROPERTY
#define CURVILINEAR_GRID_BASE_PROPERTY

#include <dune/common/exceptions.hh>

namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseProperty {

	typedef typename GridBaseType::GridStorageType  GridStorageType;
	typedef typename GridStorageType::StructuralType            StructuralType;

	static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;

public:

	CurvilinearGridBaseProperty(GridBaseType & gridbase) :
		gridbase_(gridbase),
		gridstorage_(gridbase.gridstorage())
	{

	}



	/* ***************************************************************************
	     * Section: Setting user constants
	     * ***************************************************************************/

	    /** \brief Sets the geometry tolerance. Geometry tolerance determines the precision of curvilinear volumes
	     *  Returned by the grid
	     *  */
	    void   geometryRelativeTolerance(double tolerance)  { gridstorage_.GEOMETRY_TOLERANCE = tolerance; }

	    /** \brief Retrieves the geometry tolerance */
	    double geometryRelativeTolerance() const            { return gridstorage_.GEOMETRY_TOLERANCE; }

	    bool withGhostElements() const { return gridstorage_.withGhostElements_; }

	    bool withPeriodicBoundaries() const { return gridstorage_.periodicCuboidDimensions_.size() != 0; }

	    bool isPeriodicDimension(int dim) const { return gridstorage_.periodicCuboidDimensions_[dim]; }



	    /****************************************************************************
	     * Section: (Fake) Refinement Methods of the mesh
	     *
	     * \note that there is no refinement functionality implemented at the moment
	     ****************************************************************************/

	    // Checks if the grid is located on a single process
	    bool isSerial () { return gridstorage_.mpihelper_.size() == 1; }


	    /* ***************************************************************************
	     * Section: Methods of the mesh
	     * ***************************************************************************/

	    /** Get total number of entities in a mesh  */
	    int nEntityTotal(int codim) const  { return gridstorage_.nEntityTotal_[codim]; }


	    /** Get total number of entities on this process
	     *
	     * \note Currently interpolation points are not counted towards this number.
	     * One should not use this number to loop over entities
	     * */
	    int nEntity(int codim) const  { return gridstorage_.entityAllIndexSet_[codim].size(); }


	    /** Get total number of entities of specific type on this process  */
	    int nEntity(int codim, PartitionType ptype, StructuralType btype = NO_BOUNDARY_TYPE) const
	    {
	    	return gridbase_.indexset().entityIndexSetSelect(codim, ptype, btype).size();
	    }


	    /* ***************************************************************************
	     * Section: Public Auxiliary Methods
	     * ***************************************************************************/

	    // Checks if entities of a given codim are allowed to be of a given structural type
	    // If not throws an error
	    void assertValidCodimStructuralType(int codim, StructuralType ptype) const
	    {
	    	bool pass = false;

	    	pass |= (ptype == Dune::PartitionType::InteriorEntity);
	    	pass |= (ptype == Dune::PartitionType::GhostEntity);

	    	if (codim > 0)
	    	{
	    		// Elements are not allowed to be boundaries
	    		pass |= (ptype == Dune::PartitionType::BorderEntity);
	    	}


	    	if (!pass)  {
	    		std::stringstream logstr;
	    		logstr << "CurvilinearGridBase: Unexpected codim-structtype pair codim=" << codim;
	    		logstr << " ptype=" << ptype;
	    		LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, logstr.str());
	    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected codim-structtype pair");
	    	}
	    }

private:
	    GridBaseType & gridbase_;
	    GridStorageType & gridstorage_;


};


}

}


#endif //CURVILINEAR_GRID_BASE_PROPERTY
