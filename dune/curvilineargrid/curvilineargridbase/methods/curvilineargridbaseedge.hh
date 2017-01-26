#ifndef CURVILINEAR_GRID_BASE_EDGE
#define CURVILINEAR_GRID_BASE_EDGE

#include <dune/common/exceptions.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseEdge {

	typedef typename GridBaseType::GridStorageType  GridStorageType;

    typedef typename GridStorageType::ctype  ctype;
	static const int   dimension  = GridStorageType::dimension;

    typedef typename GridStorageType::LocalIndexType            LocalIndexType;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;





public:

	CurvilinearGridBaseEdge(GridBaseType & gridbase) :
		gridstorage_(gridbase.gridstorage())
	{

	}


    /** \brief local index of one of the elements that is neighbour to this edge. */
    LocalIndexType parent(LocalIndexType localIndex) const
    {
    	assert(localIndex < gridstorage_.edge_.size());
    	return gridstorage_.edge_[localIndex].elementIndex;
    }


    /** Check if edge is a complex edge
     *   Complex edges are shared by more than two processes.
     * */
    bool isComplex(LocalIndexType localIndex) const
    {
    	auto tmp = gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].find(localIndex);
    	if (tmp == gridstorage_.processBoundaryIndexMap_[EDGE_CODIM].end())  {
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected local edge index");
    	}

    	LocalIndexType edgePBIndex = tmp->second;
    	return gridstorage_.PB2PBNeighborRank_[EDGE_CODIM][edgePBIndex].size() > 1;
    }


private:
	    GridStorageType & gridstorage_;


};


}

}


#endif //CURVILINEAR_GRID_BASE_EDGE
