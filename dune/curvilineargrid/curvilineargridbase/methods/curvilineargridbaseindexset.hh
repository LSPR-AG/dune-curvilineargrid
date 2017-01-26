#ifndef CURVILINEAR_GRID_BASE_INDEXSET
#define CURVILINEAR_GRID_BASE_INDEXSET

#include <dune/common/exceptions.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseproperty.hh>


namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseIndexSet {

	typedef typename GridBaseType::GridStorageType  GridStorageType;

    typedef typename GridStorageType::ctype  ctype;
	static const int   dimension  = GridStorageType::dimension;

    typedef typename GridStorageType::StructuralType            StructuralType;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;
    static const unsigned int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;
    static const unsigned int INTERIOR_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::InteriorBoundary;
    static const unsigned int PERIODIC_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::PeriodicBoundary;

    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;

    typedef CurvilinearGridBaseProperty<GridStorageType>  GridProperty;


public:

	CurvilinearGridBaseIndexSet(GridBaseType & gridbase) :
		gridbase_(gridbase),
		gridstorage_(gridbase.gridstorage())
	{

	}


    /* ***************************************************************************
     * Section: Selector methods for shorthand access of specific arrays
     *
     * Iterators over local indices of the mesh
     * NOTE: There are no iterators over entities because there is no entity object in the core mesh
     * There will be generic entity in the wrapper because the wrapper will define an entity object
     * ***************************************************************************/



    // Returns a link to the set of all local indices of entities of a given codimension
    // This construction allows fast iteration over entities of specific structural type
    const LocalIndexSet & entityIndexSetSelect(int codim) const {
    	return gridstorage_.entityAllIndexSet_[codim];
    }

    // Returns a link to the set of all local indices of entities of a given codimension and specific partition type and boundary type
    // This construction allows fast iteration over entities of specific structural type
    const LocalIndexSet & entityIndexSetSelect(int codim, PartitionType ptype, StructuralType boundaryType = NO_BOUNDARY_TYPE) const
    {
    	gridbase_.property().assertValidCodimStructuralType(codim, ptype);

    	// Check if this is a request for a boundary index set
    	if (boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary)  {
    		assert(codim == FACE_CODIM);  // According to convention, only faces can be boundarySegments
    		assert(ptype == Dune::PartitionType::InteriorEntity);
    		return gridstorage_.faceDomainBoundaryIndexSet_;
    	} else if (boundaryType == GridStorageType::FaceBoundaryType::InteriorBoundary)  {
    		assert(codim == FACE_CODIM);  // According to convention, only faces can be boundarySegments
    		return gridstorage_.faceInteriorBoundaryIndexSet_;
    	} else if (boundaryType == GridStorageType::FaceBoundaryType::PeriodicBoundary)  {
    		assert(codim == FACE_CODIM);  // According to convention, only faces can be boundarySegments
    		assert(ptype == Dune::PartitionType::InteriorEntity);
    		return gridstorage_.facePeriodicBoundaryIndexSet_;
    	} else if (boundaryType == NO_BOUNDARY_TYPE)  {
        	// Otherwise, this must be a request for a standard dune entity iterator
        	switch(ptype)
        	{
        	case Dune::PartitionType::InteriorEntity   : return gridstorage_.entityInternalIndexSet_[codim];          break;
        	case Dune::PartitionType::BorderEntity     : return gridstorage_.entityProcessBoundaryIndexSet_[codim];   break;
        	case Dune::PartitionType::GhostEntity      : return gridstorage_.entityGhostIndexSet_[codim];             break;
        	}
    	}

    	DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Invalid partition type / structural type");
    }


    // Returns a link to the set of all local indices of entities of a given codimension, based on Dune-convention partition type
    // This construction allows fast iteration over entities of specific structural type
    const LocalIndexSet & entityIndexSetDuneSelect(int codim, Dune::PartitionIteratorType pitype) const
    {
//    	std::cout << " Requested Dune Iter codim=" << codim << " type=" << pitype << std::endl;
//    	std::cout << " Iterator set sizes for this codim are "
//    			<< gridstorage_.entityDuneInteriorIndexSet_[codim].size() << " "
//				<< gridstorage_.entityDuneInteriorBorderIndexSet_[codim].size() << " "
//				<< gridstorage_.entityGhostIndexSet_[codim].size() << " "
//				<< gridstorage_.entityAllIndexSet_[codim].size() << std::endl;

    	const int DuneIPartition   = Dune::PartitionIteratorType::Interior_Partition;
    	const int DuneIBPartition  = Dune::PartitionIteratorType::InteriorBorder_Partition;
    	const int DuneGPartition   = Dune::PartitionIteratorType::Ghost_Partition;
    	const int DuneAllPartition = Dune::PartitionIteratorType::All_Partition;

    	switch(pitype)
    	{
    	case DuneIPartition     : return gridstorage_.entityDuneInteriorIndexSet_[codim];          break;
    	case DuneIBPartition    : return gridstorage_.entityDuneInteriorBorderIndexSet_[codim];    break;
    	case DuneGPartition     : return gridstorage_.entityGhostIndexSet_[codim];                 break;
    	case DuneAllPartition   : return gridstorage_.entityAllIndexSet_[codim];                   break;
    	default:
    	{
    		std::cout << "CurvilinearGridBase: Unexpected dune-pitype" << std::endl;
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected dune-pitype");         break;
    	}
    	}
    }


private:
    	GridBaseType & gridbase_;
	    GridStorageType & gridstorage_;

};


}

}


#endif //CURVILINEAR_GRID_BASE_INDEXSET
