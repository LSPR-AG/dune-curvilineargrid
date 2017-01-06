// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_COMMUNICATION_HH
#define DUNE_CURVGRID_COMMUNICATION_HH

#include <queue>
#include <dune/curvilineargrid/utility/allcommunication.hh>
#include <dune/curvilineargrid/curvilineargrid/entity.hh>



namespace Dune
{

  namespace CurvGrid
  {

    // Stores communication data as requested by Dune DataHandleIF
    template <class DataType>
    class CurvilinearMessageBuffer
    {
    public:
    	CurvilinearMessageBuffer()  { }

    	// Stores data
        void write(const DataType & val)  { buff_.push(val); }

        // Retrieves data from storage
        void read(DataType & val)
        {
        	assert(!buff_.empty());    // Should not attempt to read empty buffer
        	val = buff_.front();       // Retrieve first element of the queue
        	buff_.pop();               // Remove read element from the queue
        }

    private:
    	std::queue<DataType> buff_;
    };





    // Performs communication requested from gridview.communication()
    template<class Grid>
    class Communication
    {
    	typedef typename Grid::Traits Traits;
    	typedef typename Grid::ctype ctype;

    	static const int dimension   = Grid::dimension;		 //! dimension of the grid

  	    typedef typename Grid::GridStorageType  GridStorageType;
  	    typedef typename Grid::GridBaseType     GridBaseType;

    	typedef typename GridBaseType::LocalIndexType        LocalIndexType;
    	typedef typename GridBaseType::GlobalIndexType       GlobalIndexType;
    	typedef typename GridBaseType::StructuralType        StructuralType;

    	typedef typename GridBaseType::Local2LocalMap        Local2LocalMap;
    	typedef typename GridBaseType::Local2LocalIterator   Local2LocalIterator;


    	struct InterfaceSubsetType
    	{
        	enum
        	{
        		Interior_to_Ghost = 0,
        		ProcessBoundary_to_ProcessBoundary = 1,
        		ProcessBoundary_Ghost = 2,
        		Ghost_to_Interior = 3,
        		Ghost_to_ProcessBoundary = 4,
        		Ghost_to_Ghost = 5,
				PeriodicBoundary_to_PeriodicBoundary = 6
        	};
    	};

    	typedef int          InterfaceSubType;


    public:

    	Communication(GridBaseType & gridbase, MPIHelper & mpihelper)
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
    	template<class DataHandle, class Data, int codim>
    	void communicate(CommDataHandleIF< DataHandle, Data> & datahandle, InterfaceType iftype, CommunicationDirection dir) const
    	{
    		std::stringstream logstr;
    		logstr << "Started communication for codim " << codim;
    		logstr << " using interface " << interface2string(iftype, dir);
    		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, logstr.str());


    		if (codim > 0)
    		{
        		// Communication protocol ProcessBoundary -> ProcessBoundary
        		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::ProcessBoundary_to_ProcessBoundary))
        		{
        			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Using communication protocol PB->PB");
        			communicateMain<DataHandle, Data, codim>(
        				datahandle,
        				gridbase_.selectCommMap(codim, PartitionType::BorderEntity),
        				gridbase_.selectCommRankVector(codim, PartitionType::BorderEntity, PartitionType::BorderEntity)
        			);
        		}

        		// Communication protocol ProcessBoundary -> Ghost
        		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::ProcessBoundary_Ghost))
        		{
        			if (gridbase_.withGhostElements()) {
            			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Using communication protocol PB->G");
            			communicateMain<DataHandle, Data, codim>(
            				datahandle,
            				gridbase_.selectCommMap(codim, PartitionType::BorderEntity),
            				gridbase_.selectCommRankVector(codim, PartitionType::BorderEntity, PartitionType::GhostEntity)
            			);
        			}
        		}

        		// Communication protocol Periodic Boundary -> Periodic Boundary
        		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::PeriodicBoundary_to_PeriodicBoundary))
        		{
        			if (gridbase_.withPeriodicBoundaries()) {
            			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Using communication protocol Periodic->Periodic");
            			communicateMain<DataHandle, Data, codim>(
            				datahandle,
            				gridbase_.selectCommMap(codim, PERIODIC_BOUNDARY_PARTITION_TYPE),
            				gridbase_.selectCommRankVector(codim, PERIODIC_BOUNDARY_PARTITION_TYPE, PERIODIC_BOUNDARY_PARTITION_TYPE)
            			);
        			}
        		}
    		} else {
    			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Skipping PB communication for codim=" + std::to_string(codim));
    		}


    		if (gridbase_.withGhostElements()) {
        		// Communication protocol BoundaryInternal -> Ghost
        		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::Interior_to_Ghost))
        		{
        			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Using communication protocol I->G");
        			communicateMain<DataHandle, Data, codim>(
        				datahandle,
        				gridbase_.selectCommMap(codim, PartitionType::InteriorEntity),
        				gridbase_.selectCommRankVector(codim, PartitionType::InteriorEntity, PartitionType::GhostEntity)
        			);
        		}

        		// Communication protocol Ghost -> BoundaryInternal + ProcessBoundary
        		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::Ghost_to_Interior))
        		{
        			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Using communication protocol G->I");
        			communicateMain<DataHandle, Data, codim>(
        				datahandle,
        				gridbase_.selectCommMap(codim, PartitionType::GhostEntity),
        				gridbase_.selectCommRankVector(codim, PartitionType::GhostEntity, PartitionType::InteriorEntity)
        			);
        		}

        		// Communication protocol Ghost -> Ghost
        		if (allowedInterfaceSubset(iftype, dir, InterfaceSubsetType::Ghost_to_Ghost))
        		{
        			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Using communication protocol G->G");
        			communicateMain<DataHandle, Data, codim>(
        				datahandle,
        				gridbase_.selectCommMap(codim, PartitionType::GhostEntity),
        				gridbase_.selectCommRankVector(codim, PartitionType::GhostEntity, PartitionType::GhostEntity)
        			);
        		}
    		}
    	}



    	// Auxiliary method to convert interface to string for output
    	std::string interface2string(InterfaceType iftype, CommunicationDirection dir) const
    	{
    		if (iftype == Dune::InteriorBorder_InteriorBorder_Interface)  { return "IB->IB"; }
    		else if (iftype == Dune::InteriorBorder_All_Interface)
    		{
    			if (dir == Dune::ForwardCommunication)  { return "IB->ALL"; }
    			else                                    { return "ALL->IB"; }
    		}
    		else if (iftype == Dune::All_All_Interface)  { return "ALL->ALL"; }
    		else                                         { return "(unknown interface)"; }
    	}

    protected:

    	bool allowedInterfaceSubset(InterfaceType iftype, CommunicationDirection dir, InterfaceSubType istype) const
    	{
    		if (iftype == Dune::InteriorBorder_InteriorBorder_Interface)
    		{
    			return
    					(istype == InterfaceSubsetType::ProcessBoundary_to_ProcessBoundary) ||
						(istype == InterfaceSubsetType::PeriodicBoundary_to_PeriodicBoundary);
    		}
    		else if ((iftype == Dune::InteriorBorder_All_Interface)&&(dir == Dune::ForwardCommunication) )
    		{
    			return
    					(istype == InterfaceSubsetType::ProcessBoundary_to_ProcessBoundary) ||
						(istype == InterfaceSubsetType::PeriodicBoundary_to_PeriodicBoundary) ||
    				    (istype == InterfaceSubsetType::Interior_to_Ghost) ||
    				    (istype == InterfaceSubsetType::ProcessBoundary_Ghost);
    		}
    		else if ((iftype == Dune::InteriorBorder_All_Interface)&&(dir == Dune::BackwardCommunication) )
    		{
    			return
    					(istype == InterfaceSubsetType::ProcessBoundary_to_ProcessBoundary) ||
						(istype == InterfaceSubsetType::PeriodicBoundary_to_PeriodicBoundary) ||
						(istype == InterfaceSubsetType::Ghost_to_Interior) ||
						(istype == InterfaceSubsetType::Ghost_to_ProcessBoundary);
    		}
    		else if (iftype == Dune::All_All_Interface)  { return true; } // Assuming istype has any allowed value

    		// Otherwise an unsupported iftype is provided, so no communication will be done
    		DUNE_THROW(Dune::IOError, "Unexpected interface type");
    		return false;
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

    	template<class DataHandle, class Data, int codim>
    	void communicateMain(CommDataHandleIF< DataHandle, Data> & datahandle,
    		Local2LocalMap & mapSend,
			std::vector< std::vector<int> > & ranklist
    	) const
    	{
    		typedef typename DataHandle::DataType DataType;                                 // Type of data to be communicated
    		typedef typename Traits::template Codim<codim>::Entity             Entity;      // Type of the entity
    		typedef Dune::CurvGrid::CurvEntity<codim, dimension, const Grid>   EntityImpl;  // Type of the real entity

    		Dune::CurvGrid::AllCommunication allcommunicate(mpihelper_);

    		// [FIXME] May have error due to non-refreshing buffer
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
    		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Computing entities to be sent");

    		for (const auto & entityIndexPair : mapSend)
    		{
    			LocalIndexType thisEntityLocalIndex = entityIndexPair.first;
    			LocalIndexType thisEntityLocalSubIndex = entityIndexPair.second;

    			// Get Entity
    			Entity thisEntity( EntityImpl(thisEntityLocalIndex, gridbase_, Dune::PartitionIteratorType::All_Partition));

    			// Get data sizes
    			for (int iProc = 0; iProc < ranklist[thisEntityLocalSubIndex].size(); iProc++)
    			{
    				nEntityPerProcessSend[ranklist[thisEntityLocalSubIndex][iProc]]++;
    				nDataPerProcessSend[ranklist[thisEntityLocalSubIndex][iProc]] += datahandle.size(thisEntity);
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
    		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Constructing arrays for communication");


    		std::stringstream logstr;
    		logstr << "---- contents of this comm map: total=" << mapSend.size() << " ";
    		for (Local2LocalIterator iter = mapSend.begin(); iter != mapSend.end(); iter++)  {
    			logstr << "(" << (*iter).first << "," << gridbase_.entityPartitionType(codim, (*iter).first) << ") ";
    		}
    		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, logstr.str());


    		for (const auto & entityIndexPair : mapSend)
    		{
    			// Get local, sublocal and global indices of this entity (sublocal is local index of all entities of same codim and structural type)
    			LocalIndexType thisEntityLocalIndex = entityIndexPair.first;
    			LocalIndexType thisEntityLocalSubIndex = entityIndexPair.second;
    			GlobalIndexType thisEntityGlobalIndex;
    			if (!gridbase_.findEntityGlobalIndex(codim, thisEntityLocalIndex, thisEntityGlobalIndex))  {
    				std::cout << " Communication: Global index of communicating entity with local index  " << thisEntityLocalIndex << " not found" << std::endl;
    				DUNE_THROW( GridError, " Communication: Global index of communicating entity with local index  " << thisEntityLocalIndex << " not found" );
    			}

    			// Get Entity
    			//StructuralType thisEntityStructType = gridbase_.entityStructuralType(codim, thisEntityLocalIndex);
    			Entity thisEntity = Entity(EntityImpl(thisEntityLocalIndex, gridbase_, Dune::PartitionIteratorType::All_Partition));
    			datahandle.gather(gathermessagebuffer, thisEntity);


    			// Straight up read all the data that was added by the user to this entity
    			for (int iData = 0; iData < datahandle.size(thisEntity); iData++)
    			{
    				DataType thisData;
    				gathermessagebuffer.read(thisData);

    				// Add this data to send array for each process neighboring this entity
    				for (int iProc = 0; iProc < ranklist[thisEntityLocalSubIndex].size(); iProc++)
        			{
    					int thisNeighborRank = ranklist[thisEntityLocalSubIndex][iProc];
        				dataSend[displDataPerProcessSend[thisNeighborRank]++] = thisData;
        			}
    			}

    			// Record data number and global index of this entity for each process neighboring this entity
        		for (int iProc = 0; iProc < ranklist[thisEntityLocalSubIndex].size(); iProc++)
        		{
        			int thisNeighborRank = ranklist[thisEntityLocalSubIndex][iProc];
        			int entityDisplIndex = displEntityPerProcessSend[thisNeighborRank]++;

        			globalIndexSend[entityDisplIndex]    = thisEntityGlobalIndex;
        			nDataPerEntitySend[entityDisplIndex] = datahandle.size(thisEntity);
    			}
    		}


    		// 3) Communicate
    		// **************************************************************

    		// Communicate number of entities to communicate to each process
    		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Communicating");
    		allcommunicate.all2all(globalIndexSend, nEntityPerProcessSend, globalIndexRecv, nEntityPerProcessRecv);

    		// Communicate data per entity
    		allcommunicate.all2all(nDataPerEntitySend, nEntityPerProcessSend, nDataPerEntityRecv, nEntityPerProcessRecv);

    		// Communicate data
    		allcommunicate.all2all(dataSend, nDataPerProcessSend, dataRecv, nDataPerProcessRecv);


    		// 4) Scatter
    		// **************************************************************
    		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>( __FILE__, __LINE__, " -- Scattering received data");

    		int iData = 0;
    		for (int i = 0; i < globalIndexRecv.size(); i++)
    		{
    			// Get Entity
    			GlobalIndexType thisEntityGlobalIndex = globalIndexRecv[i];
    			LocalIndexType thisEntityLocalIndex;
    			if (!gridbase_.findEntityLocalIndex(codim, thisEntityGlobalIndex, thisEntityLocalIndex))  {
    				std::cout << " Communication: Local index of received entity with global index " << thisEntityGlobalIndex << " not found" << std::endl;
    				DUNE_THROW( GridError, " Communication: Local index of received entity with global index " << thisEntityGlobalIndex << " not found" );
    			}

    			//StructuralType thisEntityStructType = gridbase_.entityStructuralType(codim, thisEntityLocalIndex);
    			Entity thisEntity = Entity(EntityImpl(thisEntityLocalIndex, gridbase_, Dune::PartitionIteratorType::All_Partition));

    			// Scatter data
    			int thisNData = nDataPerEntityRecv[i];
    			for (int j = 0; j < thisNData; j++)  { scattermessagebuffer.write(dataRecv[iData++]); }

    			datahandle.scatter(scattermessagebuffer, thisEntity, thisNData);
    		}
    	}

    private:
    	MPIHelper & mpihelper_;
    	int rank_;
    	int size_;

    	GridBaseType & gridbase_;
    }; // Class



  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_COMMUNICATION_HH
