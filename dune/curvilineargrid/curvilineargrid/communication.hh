// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_COMMUNICATION_HH
#define DUNE_CURVGRID_COMMUNICATION_HH

#include <dune/curvilineargrid/utility/allcommunication.hh>



namespace Dune
{

  namespace CurvGrid
  {

    // Stores communication data as requested by Dune DataHandleIF
    template <class DataType>
    class CurvilinearMessageBuffer
    {
    public:
    	CurvilinearMessageBuffer()  { readpos_ = 0; }

    	// Stores data
        void write(const DataType & val)  { buff_.push_back(val); }

        // Retrieves data from storage
        // if (readpos_ > buff_.size())  { THROW ERROR; }
        void read(DataType & val) const   { val = buff_[readpos_++]; }

    private:
    	int readpos_;
    	std::vector<DataType> buff_;
    };





    // Performs communication requested from gridview.communication()
    template<class Grid>
    class Communication
    {
    	typedef typename Grid::Traits  Traits;

    	typedef typename Traits::LocalIndexType   LocalIndexType;
    	typedef typename Traits::GlobalIndexType  GlobalIndexType;

    	typedef typename Traits::GridBaseType     GridBaseType;
    	typedef typename Traits::GridStorageType  GridStorageType;


    	struct InterfaceSubsetType
    	{
        	enum
        	{
        		Internal_Ghost = 0,
        		ProcessBoundary_ProcessBoundary = 1,
        		ProcessBoundary_Ghost = 2,
        		Ghost_Internal = 3,
        		Ghost_ProcessBoundary = 4,
        		Ghost_Ghost = 5
        	};
    	};

    	typedef int          InterfaceSubType;






    public:

    	Communication(const GridBaseType & gridbase, MPIHelper & mpihelper)
    		: mpihelper_(mpihelper),
    		  gridbase_(gridbase)
    	{
    		rank_ = mpihelper.rank();
    		size_ = mpihelper.size();
    	}


    	// Wrapper Communication Algorithm
    	// 1) Loop over all codim
    	// 1.1) Check if this codim is allowed by DataHandle
    	// 1.2) For each InterfaceSubset, check if it is consistent with InterfaceType and CommunicationDirection
    	// 1.3) If it is, call main communication protocol main_communicate(codim, mapSend, ranklistSend)
    	template<class DataHandle, int codim>
    	void communicateWrapper(DataHandle& datahandle, InterfaceType iftype, CommunicationDirection dir, int level) const
    	{
    		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::ProcessBoundary_ProcessBoundary))
    		{
    			communicateMain<DataHandle, codim>(
    				datahandle,
    				gridstorage_.processBoundaryIndexMap_[codim],
    				gridstorage_.PB2PBNeighborRank_[codim]
    			);
    		}

    		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::ProcessBoundary_Ghost))
    		{
    			communicateMain<DataHandle, codim>(
    				datahandle,
    				gridstorage_.processBoundaryIndexMap_[codim],
    				gridstorage_.PB2GNeighborRank_[codim]
    			);
    		}

    		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::Internal_Ghost))
    		{
    			communicateMain<DataHandle, codim>(
    				datahandle,
    				gridstorage_.boundaryInternalEntityIndexMap_[codim],
    				gridstorage_.BI2GNeighborRank_[codim]
    			);
    		}

    		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::Ghost_Internal))
    		{
    			communicateMain<DataHandle, codim>(
    				datahandle,
    				gridstorage_.ghostIndexMap_[codim],
    				gridstorage_.G2BIPBNeighborRank_[codim]
    			);
    		}

    		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::Ghost_Ghost))
    		{
    			communicateMain<DataHandle, codim>(
    				datahandle,
    				gridstorage_.ghostIndexMap_[codim],
    				gridstorage_.G2GNeighborRank_[codim]
    			);
    		}
    	}


    	// Main Communication Algorithm
    	// 1) Iterate over provided map
    	// 1.1) mapSend.first -> LocalIndex, GlobalIndex -> Entity -> Data from DataHandle
    	// 1.2) mapSend.second -> LocalStructIndex -> neighbor ranks from ranklist
    	// 1.3) Data Structures to communicate
    	//        * nEntitiesPerProcess
    	//        * nDataPerProcess
    	//        * nDataPerEntity
    	//        * entityGlobalIndex
    	//        * Data - Communicated via Allcommunicate
    	// 2) On receiving end
    	// 2.1) GlobalIndex -> LocalIndex -> Entity -> scatter via DataHandle
    	template<class DataHandle, int codim>
    	void communicationMain(DataHandle& datahandle,
    		Local2LocalMap mapSend,
			std::vector< std::vector<int> > ranklist
    	)
    	{

    	}










    	// [FIXME] Choose if we want entity or entityImpl
    	// [FIXME] Pass on IteratorRange somehow

    	template<class DataHandle, int codim>
    	void communicate (DataHandle& datahandle, InterfaceType iftype, CommunicationDirection dir, int level) const
    	{
    		typedef typename DataHandle::DataType DataType;             // Type of data to be communicated
    		typedef typename Traits::Codim<codim>::Entity  EntityType;  // Type of the entity

    		Dune::CurvGrid::AllCommunication allcommunicate(mpihelper_);

    		CurvilinearMessageBuffer<DataType> gathermessagebuffer;
    		CurvilinearMessageBuffer<DataType> scattermessagebuffer;

    		int nDataSend = 0;
    		int nEntitySend = 0;

    		std::vector<int> nEntityPerProcessSend(size_);
    		std::vector<int> nEntityPerProcessRecv(size_);

    		std::vector<GlobalIndexType> globalIndexSend;
    		std::vector<GlobalIndexType> globalIndexRecv;

    		std::vector<int> nDataPerEntitySend;
    		std::vector<int> nDataPerEntityRecv;

    		std::vector<int> nDataPerProcessSend(size_, 0);
    		std::vector<int> nDataPerProcessRecv(size_);

    		std::vector<int> displEntityPerProcessSend(size_);
    		std::vector<int> displDataPerProcessSend(size_);

    		std::vector<DataType> dataSend;
    		std::vector<DataType> dataRecv;


    		// 1) First loop over all entities of the interface
    		// Compute total number of entities to be communicated to each process
    		// ********************************************************


    		for (int i = 0; i < 10; i++)
    		{
    			const EntityType & thisEntity = entity(i);

    			LocalIndexType thisLocalIndex = thisEntity.localIndex();
    			std::vector<int> thisNeighborRankSet = gridbase_.processBoundaryNeighborRankSet(codim, thisLocalIndex);

    			for (int iProc = 0; iProc < thisNeighborRankSet.size(); iProc++)
    			{
    				nEntityPerProcessSend[thisNeighborRankSet[iProc]]++;
    				nDataPerProcessSend[thisNeighborRankSet[iProc]] += datahandle.size(thisEntity);
    			}
    		}

    		// 1.1) Construct send-displacements
    		for (int i = 0; i < size_; i++)  {
    			displEntityPerProcessSend[i]  = (i == 0) ? 0 : displEntityPerProcessSend[i - 1] + nEntityPerProcessSend[i - 1];
    			displDataPerProcessSend[i]    = (i == 0) ? 0 : displDataPerProcessSend[i - 1]   + nDataPerProcessSend[i - 1];

    			nEntitySend += nEntityPerProcessSend[i];
    			nDataSend += nDataPerProcessSend[i];
    		}
    		globalIndexSend.resize(nEntitySend);
    		nDataPerEntitySend.resize(nEntitySend);
    		dataSend.resize(nDataSend);


    		// 2) Second loop over all entities of the interface
    		// Read and fill in data
    		// ********************************************************

    		for (int i = 0; i < 10; i++)
    		{
    			const EntityType & thisEntity = entity(i);

    			LocalIndexType thisLocalIndex = thisEntity.localIndex();
    			GlobalIndexType thisGlobalIndex;
    			if (!gridbase_.findEntityGlobalIndex(codim, thisLocalIndex, thisGlobalIndex))  {  /* THROW ERROR */ }
    			std::vector<int> thisNeighborRankSet = gridbase_.processBoundaryNeighborRankSet(codim, thisLocalIndex);

    			datahandle.gather(gathermessagebuffer, thisEntity);


    			// Straight up read all the data that was added by the user to this entity
    			for (int iData = 0; iData < datahandle.size(thisEntity); iData++)
    			{
    				DataType thisData;
    				gathermessagebuffer.read(thisData);

    				// Add this data to send array for each process neighboring this entity
        			for (int iProc = 0; iProc < thisNeighborRankSet.size(); iProc++)
        			{
        				dataSend[displDataPerProcessSend[thisNeighborRankSet[iProc]]++] = thisData;
        			}
    			}

    			// Record data number and global index of this entity for each process neighboring this entity
        		for (int iProc = 0; iProc < thisNeighborRankSet.size(); iProc++)
        		{
        			int entityDisplIndex = displEntityPerProcessSend[thisNeighborRankSet[iProc]]++;

        			globalIndexSend[entityDisplIndex]    = thisGlobalIndex;
        			nDataPerEntitySend[entityDisplIndex] = datahandle.size(thisEntity);
    			}
    		}


    		// 3) Communicate
    		// **************************************************************

    		// Communicate entity global index
    		allcommunicate.communicate(globalIndexSend, nEntityPerProcessSend, globalIndexRecv, nEntityPerProcessRecv);

    		// Communicate data per entity
    		allcommunicate.communicate(nDataPerEntitySend, nEntityPerProcessSend, nDataPerEntityRecv, nEntityPerProcessRecv);

    		// Communicate data per entity
    		allcommunicate.communicate(dataSend, nDataPerProcessSend, dataRecv, nDataPerProcessRecv);


    		// 4) Refill
    		// **************************************************************

    		int iData = 0;
    		for (int i = 0; i < globalIndexRecv.size(); i++)
    		{
    			GlobalIndexType thisGlobalIndex = globalIndexRecv[i];
    			int thisNData = nDataPerEntityRecv[i];

    			const EntityType & thisEntity = entity(thisGlobalIndex);

    			for (int j = 0; j < thisNData; j++)  { scattermessagebuffer.write(dataRecv[iData++]); }

    			datahandle.scatter(scattermessagebuffer, thisEntity);
    		}
    	}



    protected:

    	bool allowedInterfaceSubset(InterfaceType iftype, CommunicationDirection dir, InterfaceSubType istype)
    	{
    		if (InterfaceType == Dune::InteriorBorder_InteriorBorder_Interface)
    		{
    			return (istype == InterfaceSubsetType::ProcessBoundary_ProcessBoundary);
    		}
    		else if ((iftype == Dune::InteriorBorder_All_Interface)&&(dir == Dune::ForwardCommunication) )
    		{
    			return (istype == InterfaceSubsetType::ProcessBoundary_ProcessBoundary) ||
    				   (istype == InterfaceSubsetType::Internal_Ghost) ||
    				   (istype == InterfaceSubsetType::ProcessBoundary_Ghost);
    		}
    		else if ((iftype == Dune::InteriorBorder_All_Interface)&&(dir == Dune::BackwardCommunication) )
    		{
    			return (istype == InterfaceSubsetType::ProcessBoundary_ProcessBoundary) ||
    				   (istype == InterfaceSubsetType::Ghost_Internal) ||
    				   (istype == InterfaceSubsetType::Ghost_ProcessBoundary);
    		}
    		else if (iftype == Dune::All_All_Interface)  { return true; } // Assuming istype has any allowed value

    		// Otherwise an unsupported iftype is provided, so no communication will be done
    		return false;
    	}



    private:
    	MPIHelper & mpihelper_;
    	int rank_;
    	int size_;

    	const GridStorageType & gridstorage_;

    	const GridBaseType & gridbase_;


    }; // Class



  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_COMMUNICATION_HH
