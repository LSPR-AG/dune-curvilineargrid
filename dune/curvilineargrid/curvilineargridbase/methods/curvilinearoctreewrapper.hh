#ifndef CURVILINEAR_GRID_BASE_OCTREE_WRAPPER
#define CURVILINEAR_GRID_BASE_OCTREE_WRAPPER

#include <dune/common/exceptions.hh>

namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseOctreeWrapper {

	typedef typename GridBaseType::GridStorageType  GridStorageType;

	typedef typename GridStorageType::GlobalCoordinate                 GlobalCoordinate;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

public:

	CurvilinearGridBaseOctreeWrapper(GridBaseType & gridbase) :
		gridbase_(gridbase),
		gridstorage_(gridbase.gridstorage())
	{

	}


    /** \brief Minimal bounding box for set of elements on this process */
    void processBoundingBox(GlobalCoordinate & center, GlobalCoordinate & extent) const
    {
        center = gridstorage_.boundingBoxCenter_;
        extent = gridstorage_.boundingBoxExtent_;
    }


    /** Searches for elements containing this global coordinate
     *
     * \param[in]    globalC         the global coordinate of the point looked for
     * \param[in]    containerIndex  returns all element indices in which the point was found.
     * In principle there can be more than one if the point is close to the boundary between elements.
     * \param[in]    localC          returns local coordinates associated with this point within all elements it is found in
     *
     * \returns      Whether the point found at all on this process. If not, return arrays should be empty
     *
     *
     * Algorithm:
     *
     * 1) Check if this point is inside the bounding box
     * 2) If yes, check if this point is inside OCTree
     * 3) OCTree returns a list of candidate elements - those located in the lowest level octant the point is found in
     * 4) Loop over candidate elements, use CurvilinearGeometry search functionality
     * 5) If the point is found inside by CurvilinearGeometry, it also returns the associated local coordinate at no extra cost
     * 6) All the elements in which the point is found, as well as the associated local coordinates are returned
     *
     * TODO: PBE file allows femaxx to count time. Use alternative in Dune?
     * */
    bool locateCoordinate(const GlobalCoordinate & globalC, std::vector<int> & containerIndex, std::vector<GlobalCoordinate> & localC) const {
    	DUNE_THROW(Dune::IOError, "Using non-initialised OCTree");

        // If the point is not even in the process bounding box, then the coordinate is definitely not on this process
        if (isInsideProcessBoundingBoxGracious(globalC))
        {
            //pbe_start(132, "Octree traversal");

            // Get list of indices of elements which may contain the coordinate
            // This corresponds to the elements at the lowest level of the OCTree
            // Note that no internal element search is performed by the OCTree itself
            int nNodeVisited = 0;
            std::vector<int> elementIndices;
            gridstorage_.octree_->findNode(globalC, elementIndices, nNodeVisited, &isInsideProcessBoundingBoxGracious);



            // Loop over candidate elements and search for local coordinate using internal Geometry mechanism
            for (unsigned int i = 0; i < elementIndices.size(); i++) {
            	auto thisGeometry = gridbase_.entity().template geometry<ELEMENT_CODIM>(elementIndices[i]);

                GlobalCoordinate thisLocalC;
                bool isInside = thisGeometry.local( globalC, thisLocalC );

                // If the point is found inside the geometry, add the element index and the found local coordinate to the output
                if (isInside)
                {
                    containerIndex.push_back(elementIndices[i]);
                    localC.push_back(thisLocalC);
                }
            }

            // pbe_stop(132);
            // rDebug("find_tets_by_point: nof_found=%d, nNodeVisited=%d", static_cast<int>(nodes.size()), nNodeVisited);

            return (containerIndex.size() > 0);
        } else {
            return false;
        }
    }

protected:

    // Checks if given point fits into the bounding box of the mesh
    bool isInsideProcessBoundingBoxGracious(const GlobalCoordinate & point) const {
        const double grace_tolerance = 1e-13;
        GlobalCoordinate center, extent;

        bool isInside = true;

        isInside &= fabs(gridstorage_.boundingBoxCenter_[0] - point[0]) <= gridstorage_.boundingBoxExtent_[0] * (1.0 + grace_tolerance);
        isInside &= fabs(gridstorage_.boundingBoxCenter_[1] - point[1]) <= gridstorage_.boundingBoxExtent_[1] * (1.0 + grace_tolerance);
        isInside &= fabs(gridstorage_.boundingBoxCenter_[2] - point[2]) <= gridstorage_.boundingBoxExtent_[2] * (1.0 + grace_tolerance);

        return isInside;
    }



private:
	    GridBaseType & gridbase_;
	    GridStorageType & gridstorage_;


};


}

}


#endif //CURVILINEAR_GRID_BASE_OCTREE_WRAPPER
