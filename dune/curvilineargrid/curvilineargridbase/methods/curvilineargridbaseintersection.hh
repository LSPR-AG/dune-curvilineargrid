#ifndef CURVILINEAR_GRID_BASE_INTERSECTION
#define CURVILINEAR_GRID_BASE_INTERSECTION

#include <assert.h>

#include <dune/common/exceptions.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>

#include <dune/curvilineargrid/curvilineargridbase/methods/curvilineargridbaseproperty.hh>

namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseIntersection {

	typedef typename GridBaseType::GridStorageType  GridStorageType;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

    static const unsigned int NO_BOUNDARY_TYPE     = GridStorageType::FaceBoundaryType::None;
    static const unsigned int DOMAIN_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::DomainBoundary;
    static const unsigned int INTERIOR_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::InteriorBoundary;
    static const unsigned int PERIODIC_BOUNDARY_TYPE = GridStorageType::FaceBoundaryType::PeriodicBoundary;

    typedef CurvilinearGridBaseProperty<GridStorageType> GridProperty;
    typedef CurvilinearGridBaseEntity<GridStorageType> GridEntity;

public:

	CurvilinearGridBaseIntersection(GridBaseType & gridbase) :
		gridbase_(gridbase),
		gridstorage_(gridbase.gridstorage())
	{

	}



    /** Get boundary segment index of this entity if it is a Domain Boundary Face */
    LocalIndexType boundarySegmentIndex(LocalIndexType localIndex) const
    {
    	bool isDB = (gridstorage_.face_[localIndex].boundaryType == GridStorageType::FaceBoundaryType::DomainBoundary);
    	bool isPeriodic = (gridstorage_.face_[localIndex].boundaryType == GridStorageType::FaceBoundaryType::PeriodicBoundary);

    	assert(isDB || isPeriodic);  // Periodic boundaries are also domain boundaries
    	return gridstorage_.boundarySegmentIndexMap_.at(localIndex);
    }


    /** Get DomainBoundaryFace index given its boundary segment index */
    // [TODO] The map in this case can be replaced by a simple array/vector, since index is continuous
    LocalIndexType boundarySegment2LocalIndex(LocalIndexType boundarySegmentIndex) const
    {
    	assert(boundarySegmentIndex < gridstorage_.boundarySegment2LocalIndexMap_.size());  // Boundary segment index is continuous
    	return gridstorage_.boundarySegment2LocalIndexMap_.at(boundarySegmentIndex);
    }




    /** Get Boundary type of an face */
    StructuralType boundaryType(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.face_.size());
    	return gridstorage_.face_[localIndex].boundaryType;
    }



    /** \brief Check if there is an outer neighbor defined for the given face local index
     *  - expect false for domain boundaries and serial process boundaries, true for everything else
     * */
    bool checkOuterNeighbor(LocalIndexType localIndex) const {
    	assert(localIndex < gridstorage_.face_.size());
    	return gridstorage_.face_[localIndex].element2Index >= 0;
    }

    /** \brief local index of the element that is neighbor to this face
     *
     * \param[in]    localIndex               local index of this face
     * \param[in]    internalNeighborIndex    {0,1} - determines which of the two neighbors to address
     *
     * \return local index of the element that is neighbor to this face.
     * NOTE: If the neighbor is a GHOST ELEMENT, the returned int will be for the GhostElementLocalIndex
     * and not the ElementIndex
     *
     * Conventions of internalNeighborIndex for face types
     * * For Domain Boundary there is only one neighbor
     * * For Process Boundary 2nd neighbor is always the Ghost Element
     * * For Internal Face there is no convention on order of the neighbors
     *
     *  */
    LocalIndexType neighborElement(LocalIndexType localIndex, InternalIndexType internalNeighborIndex) const
    {
    	LocalIndexType rez;

    	assert(localIndex < gridstorage_.face_.size());

        switch(internalNeighborIndex)
        {
        case 0 : rez = gridstorage_.face_[localIndex].element1Index;  break;
        case 1 : rez = gridstorage_.face_[localIndex].element2Index;  break;
        default:
        {
        	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearPostConstructor: Unexpected neighbor subentity index =" + std::to_string(internalNeighborIndex));
        	DUNE_THROW(Dune::IOError, "CurvilinearGrid: faceNeighbor() unexpected neighbor index");  break;
        }
        }

        // User should not request non-existing local index
        if (rez < 0) {
        	std::stringstream logstr;
        	logstr << "CurvilinearGrid: Face does not appear to have a neighbor. localIndex=" << localIndex
        			<< ", internal=" << internalNeighborIndex
					<< ", resultIndex=" << rez
					<< ", facePType=" << gridstorage_.face_[localIndex].ptype
					<< ", boundaryType=" << gridstorage_.face_[localIndex].boundaryType;
        	DUNE_THROW(Dune::IOError, logstr.str());
        }

        return rez;
    }


    // IMPORTANT NOTE: IN CASE OF PERIODIC FACES, THE INTERSECTION AS SEEN FROM INSIDE AND OUTSIDE ARE DIFFERENT ENTITIES
    InternalIndexType subIndexInNeighbor(LocalIndexType localIndex, InternalIndexType internalNeighborIndex) const {
    	InternalIndexType rez;

    	assert(localIndex < gridstorage_.face_.size());

        switch(internalNeighborIndex)
        {
        case 0 : rez = gridstorage_.face_[localIndex].element1SubentityIndex;  break;
        case 1 : rez = gridstorage_.face_[localIndex].element2SubentityIndex;  break;
        default:
        {
        	LoggingMessage::template write<LOG_MSG_PERSISTENT>( __FILE__, __LINE__, "CurvilinearPostConstructor: Unexpected neighbor subentity index =" + std::to_string(internalNeighborIndex));
        	DUNE_THROW(Dune::IOError, "CurvilinearGrid: faceNeighbor() unexpected neighbor index");  break;
        }
        }

        // User should not request details of non-existing neighbors. Use checkFaceOuterNeighbor()
        assert(rez >= 0);

        return rez;
    }


    // Given the periodic face index, finds the associated ghost face index
    // [TODO] This approach is cumbersome and error-prone. Implement periodicGhost indexSet and check if index is in it
    LocalIndexType periodicNeighborFace(LocalIndexType periodicfaceIndex) const {
    	assert(gridbase_.property().withGhostElements());
    	assert(boundaryType(periodicfaceIndex) == GridStorageType::FaceBoundaryType::PeriodicBoundary);
    	assert(checkOuterNeighbor(periodicfaceIndex));

		InternalIndexType subind = subIndexInNeighbor(periodicfaceIndex, 1);
		LocalIndexType elemind = neighborElement(periodicfaceIndex, 1);
		LocalIndexType neighborFaceIndex = gridbase_.entity().subentityLocalIndex (elemind, ELEMENT_CODIM, FACE_CODIM, subind);

		bool isGhost = gridbase_.entity().partitionType(FACE_CODIM, neighborFaceIndex) == Dune::PartitionType::GhostEntity;
		bool isPeriodic = boundaryType(neighborFaceIndex) == PERIODIC_BOUNDARY_TYPE;

		// Given ghost elements, for each periodic face there must be an associated (ghost/local periodic) face
		assert(isGhost || isPeriodic);
		return neighborFaceIndex;
    }


    // Given the index of ghost face, looks for the associated periodic face index. If there isnt any, returns false. If there is, returns true, and the periodic index
    // [TODO] This approach is cumbersome and error-prone. Implement periodicGhost indexSet and check if index is in it
    bool findGhostPeriodicNeighborFace(LocalIndexType ghostfaceindex, LocalIndexType & periodicfaceIndex) const {
    	assert(gridbase_.property().withGhostElements());
    	assert(gridbase_.entity().partitionType(FACE_CODIM, ghostfaceindex) == Dune::PartitionType::GhostEntity);

    	if (checkOuterNeighbor(ghostfaceindex)) {
    		InternalIndexType subind = subIndexInNeighbor(ghostfaceindex, 1);
    		LocalIndexType elemind = neighborElement(ghostfaceindex, 1);
    		periodicfaceIndex = gridbase_.entity().subentityLocalIndex (elemind, ELEMENT_CODIM, FACE_CODIM, subind);

    		// If the ghost face at all has a neighbor defined, it must be because it is a periodic ghost. Other ghosts do not have neighbors
    		assert(boundaryType(periodicfaceIndex) == GridStorageType::FaceBoundaryType::PeriodicBoundary);
    		return true;
    	} else {
    		// This is a regular ghost face, it is not a periodic boundary ghost
    		return false;
    	}
    }


    // Extracts the permutation index by which the indices of the periodic face ON THIS PROCESS have to be permuted to be conformal to the periodic face ON THE NEIGHBOR PROCESS
    unsigned int periodicPermutationInner(LocalIndexType faceIndex) const {
    	auto periodicFaceIter = gridstorage_.periodicBoundaryIndexMap_.find(faceIndex);
    	assert(periodicFaceIter != gridstorage_.periodicBoundaryIndexMap_.end());  // Should not run this method on a non-periodic face
    	assert(periodicFaceIter->second < gridstorage_.periodicFaceMatchPermutationIndexInner_.size());
    	return gridstorage_.periodicFaceMatchPermutationIndexInner_[periodicFaceIter->second];
    }


    // Extracts the permutation index by which the indices of the periodic face ON THE NEIGHBOR PROCESS have to be permuted to be conformal to the periodic face ON THIS PROCESS
    // [TODO] These two permutations are inverses of each other. Having both is excessive
    //   Only communicate the Inner permutation, and obtain the outer by simple permutation inversion routine. Implement routine in CurvilinearGeometryHelper
    unsigned int periodicPermutationOuter(LocalIndexType faceIndex) const {
    	auto periodicFaceIter = gridstorage_.periodicBoundaryIndexMap_.find(faceIndex);
    	assert(periodicFaceIter != gridstorage_.periodicBoundaryIndexMap_.end());  // Should not run this method on a non-periodic face
    	assert(periodicFaceIter->second < gridstorage_.periodicFaceMatchPermutationIndexOuter_.size());
    	return gridstorage_.periodicFaceMatchPermutationIndexOuter_[periodicFaceIter->second];
    }


private:
    	GridBaseType & gridbase_;
	    GridStorageType & gridstorage_;

};


}

}


#endif //CURVILINEAR_GRID_BASE_INTERSECTION
