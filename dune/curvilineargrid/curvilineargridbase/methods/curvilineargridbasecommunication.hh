#ifndef CURVILINEAR_GRID_BASE_COMMUNICATION
#define CURVILINEAR_GRID_BASE_COMMUNICATION

#include <assert.h>

#include <dune/common/exceptions.hh>


namespace Dune {

namespace CurvGrid {


template<typename GridBaseType>
class CurvilinearGridBaseCommunication {

	typedef typename GridBaseType::GridStorageType  GridStorageType;

    typedef typename GridStorageType::GlobalIndexType           GlobalIndexType;
    typedef typename GridStorageType::LocalIndexType            LocalIndexType;
	typedef typename GridStorageType::InternalIndexType	        InternalIndexType;
	typedef typename GridStorageType::StructuralType            StructuralType;
	typedef typename GridStorageType::PhysicalTagType           PhysicalTagType;

    typedef typename GridStorageType::Local2LocalMap            Local2LocalMap;
    typedef typename GridStorageType::Global2LocalMap           Global2LocalMap;
    typedef typename GridStorageType::Local2LocalIterator       Local2LocalIterator;
    typedef typename GridStorageType::Global2LocalIterator      Global2LocalIterator;
    typedef typename GridStorageType::Global2LocalConstIterator      Global2LocalConstIterator;
    typedef typename GridStorageType::LocalIndexSet             LocalIndexSet;
    typedef typename GridStorageType::IndexSetIterator          IndexSetIterator;

    typedef typename GridStorageType::EntityNeighborRankVector  EntityNeighborRankVector;

    // Codimensions of entity types for better code readability
    static const int   VERTEX_CODIM   = GridStorageType::VERTEX_CODIM;
    static const int   EDGE_CODIM     = GridStorageType::EDGE_CODIM;
    static const int   FACE_CODIM     = GridStorageType::FACE_CODIM;
    static const int   ELEMENT_CODIM  = GridStorageType::ELEMENT_CODIM;

public:

	CurvilinearGridBaseCommunication(GridBaseType & gridbase) :
		gridstorage_(gridbase.gridstorage())
	{

	}


	/** Get the neighbour ranks of this communication entity  */
	std::vector<int> & entityNeighborRankSet(
		int codim, LocalIndexType localIndex, PartitionType ptypeSend, PartitionType ptypeRecv, StructuralType boundarytype
	)
	{
		//typedef typename std::map<LocalIndexType,LocalIndexType>::const_iterator Local2LocalConstIter;

		Local2LocalMap & tmpMap = selectCommMap(codim, ptypeSend, boundarytype);
		const auto tmpiter = tmpMap.find(localIndex);

		if (tmpiter == tmpMap.end())  { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected local index"); }

		LocalIndexType entityLocalSubsetIndex = tmpiter->second;

		EntityNeighborRankVector & rezMapSet = selectCommRankVector(codim, ptypeSend, ptypeRecv, boundarytype);

		if (entityLocalSubsetIndex >= rezMapSet.size()) {
			std::cout << gridstorage_.mpihelper_.rank() << " fail to find comm protocol for codim=" << codim
					<< " pSend=" << ptypeSend
					<< " pRecv=" << ptypeRecv
					<< " btype=" << boundarytype
					<< " ::: real_size =" << rezMapSet.size()
					<< " requested=" << entityLocalSubsetIndex << std::endl;
		}

		assert(entityLocalSubsetIndex < rezMapSet.size());			// Communicating entity must have defined comm set
		assert(rezMapSet[entityLocalSubsetIndex].size() != 0);		// If entity is communicating, it must have at least one destination

		return rezMapSet[entityLocalSubsetIndex];
	}


	Local2LocalMap & selectCommMap(int codim, Dune::PartitionType ptype, StructuralType boundarytype)
	{
		switch (ptype)
		{
		case Dune::PartitionType::InteriorEntity :  {
			if (boundarytype == GridStorageType::FaceBoundaryType::PeriodicBoundary) {
				assert(codim == FACE_CODIM);
				return gridstorage_.periodicBoundaryIndexMap_;
			} else {
				return gridstorage_.boundaryInternalEntityIndexMap_[codim];
			}
		} break;
		case Dune::PartitionType::BorderEntity   :  return gridstorage_.processBoundaryIndexMap_[codim];         break;
		case Dune::PartitionType::GhostEntity    :  return gridstorage_.ghostIndexMap_[codim];                   break;
		default                                  :  DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected comm structural type");  break;
		}
	}


	EntityNeighborRankVector & selectCommRankVector(
			int codim, Dune::PartitionType ptypesend, Dune::PartitionType ptyperecv, StructuralType boundarytype)
	{
		// Can only communicate over these 3 PartitionTypes
		//assertValidCodimStructuralType(codim, ptypesend);
		//assertValidCodimStructuralType(codim, ptyperecv);

		switch (ptypesend)
		{
			case Dune::PartitionType::InteriorEntity :   // Internal -> Ghost protocol
			{
				if (boundarytype == GridStorageType::FaceBoundaryType::PeriodicBoundary) {
					assert(codim == FACE_CODIM);  // Other codim not implemented for this protocol yet
					return gridstorage_.PERB2PERBNeighborRank_[codim];
				} else {
					assert(ptyperecv == Dune::PartitionType::GhostEntity);
					return gridstorage_.BI2GNeighborRank_[codim];
				}
			} break;
			case Dune::PartitionType::BorderEntity   :   // PB -> PB and PB -> Ghost protocols
			{
				if      (ptyperecv == Dune::PartitionType::BorderEntity)  { return gridstorage_.PB2PBNeighborRank_[codim]; }
				else if (ptyperecv == Dune::PartitionType::GhostEntity)   { return gridstorage_.PB2GNeighborRank_[codim]; }
				else { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected comm partition type"); }
			} break;
			case Dune::PartitionType::GhostEntity    :   // Ghost -> (Internal & PB) and Ghost -> Ghost protocols
			{
				if      (ptyperecv == Dune::PartitionType::InteriorEntity) { return gridstorage_.G2BIPBNeighborRank_[codim]; }
				else if (ptyperecv == Dune::PartitionType::BorderEntity)   { return gridstorage_.G2BIPBNeighborRank_[codim]; }
				else if (ptyperecv == Dune::PartitionType::GhostEntity)    { return gridstorage_.G2GNeighborRank_[codim]; }
				else { DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected comm partition type"); }
			} break;
			default: DUNE_THROW(Dune::IOError, "CurvilinearGridBase: Unexpected comm partition type");  break;
		}
	}

private:
	    GridStorageType & gridstorage_;

};


}

}


#endif //CURVILINEAR_GRID_BASE_COMMUNICATION
