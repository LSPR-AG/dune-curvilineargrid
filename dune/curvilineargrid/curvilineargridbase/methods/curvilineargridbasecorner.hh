#ifndef CURVILINEAR_GRID_BASE_CORNER
#define CURVILINEAR_GRID_BASE_CORNER

#include <dune/common/exceptions.hh>

#include <dune/curvilineargeometry/interpolation/curvilineargeometryhelper.hh>


namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseCorner {

    typedef typename GridBaseType::ctype  ctype;
	static const int   dimension  = GridBaseType::dimension;

	typedef typename GridBaseType::GridStorageType  GridStorageType;

    typedef typename GridStorageType::LocalIndexType            LocalIndexType;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;





public:

	CurvilinearGridBaseCorner(GridBaseType & gridbase) :
		gridstorage_(gridbase.gridstorage())
	{

	}


    /** Get unique local index that enumerates only corners, disregarding other interpolatory vertices
     * The only purpose of this index is to satisfy dune-standard for consecutive corner index, it is not used internally in the grid
     * */
    LocalIndexType uniqueLocalIndex(LocalIndexType localVertexIndex) const {

    	const auto tmp = gridstorage_.cornerIndexMap_.find(localVertexIndex);
    	if (tmp == gridstorage_.cornerIndexMap_.end())
    	{
    		std::cout << "Error: CurvilinearGridBase: CornerUniqueIndex: Vertex with localIndex=" << localVertexIndex << " is not marked as a corner" << std::endl;
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected vertex local index for a corner");
    	}

    	return tmp->second;
    }


    // The reverse of cornerUniqueLocalIndex. Given the unique corner in
    LocalIndexType unique2LocalIndex(LocalIndexType localCornerIndex) const {


    	const auto tmp = gridstorage_.cornerIndexMapRev_.find(localCornerIndex);

    	if (tmp == gridstorage_.cornerIndexMapRev_.end())
    	{
    		std::cout << "Error: CurvilinearGridBase: CornerUnique2LocalIndex: Vertex with cornerIndex=" << localCornerIndex << " is not marked as a corner" << std::endl;
    		DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected vertex local index for a corner");
    	}

    	return tmp->second;
    }



private:
	    GridStorageType & gridstorage_;


};


}

}


#endif //CURVILINEAR_GRID_BASE_CORNER
