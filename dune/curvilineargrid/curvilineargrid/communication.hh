// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CURVGRID_DECLARATION_HH
#define DUNE_CURVGRID_DECLARATION_HH

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
    struct Communication
    {

    	template<class DataHandle, int codim>
    	void communicate (DataHandle& datahandle, InterfaceType iftype, CommunicationDirection dir, int level) const
    	{
    		// Type of data to be communicated
    		typedef typename DataHandle::DataType DataType;

    		CurvilinearMessageBuffer<DataType> gathermessagebuffer;


    		// rank -> communicationIndex -> neighborEntityGlobalIndex
    		std::vector< std::vector< GlobalIndexType> > neighborEntityGlobalIndex(size_);

    		// rank -> communicationIndex -> std::vector<data to send>
    		std::vector< std::vector< std::vector< DataType > > > neighborEntityData(size_);



    		// 1) Loop over all entities of the interface
    		for (int i = 0; i < 10; i++)
    		{
    			const EntityType & thisEntity = entity(i);

    			datahandle.gather(gathermessagebuffer, thisEntity);

    			int neighborEntityRank = entity_neighbor_rank(i);
    			int neighborEntityGlobalIndex = entity_neighbor_index(i);

    			std::vector< DataType > dataToSend;

    			// Straight up read all the data that was added by the user to this entity
    			for (int iData = 0; iData < datahandle.size(thisEntity); iData++)
    			{
    				DataType thisData;
    				messagebuffer.read(thisData);
    				dataToSend.push_back(thisData);
    			}

    			neighborEntityGlobalIndex[neighborEntityRank].push_back(neighborEntityGlobalIndex);
    			neighborEntityData[neighborEntityRank].push_back(dataToSend);
    		}

    		// 2) Communicate data
    		// 2.1) Communicate number of entities with non-zero data to be communicated
    		std::vector<int> entitiesPerRankSend;
    		std::vector<int> entitiesPerRankRecv(size_);

    		for (int iProc = 0; iProc < size_; iProc++)  { entitiesPerRankSend.push_back(neighborEntityData[iProc].size()); }
    		int MPI_Alltoall(entitiesPerRankSend.data(), 1, MPI_INT, reinterpret_cast<int*>(entitiesPerRankRecv.data()), 1, MPI_INT, comm);

    		// 2.2) Communicate globalindex for each entity
    		// 2.3) Communicate datasize for each entity
    		int nEntityRecv = 0;
    		std::vector< GlobalIndexType > sendEntityGlobalIndex;
    		std::vector< GlobalIndexType > recvEntityGlobalIndex;
    		std::vector< int > sendEntityDataSize;
    		std::vector< int > recvEntityDataSize;
    		std::vector<int> sdispls, rdispls;

    		for (int iProc = 0; iProc < size_; iProc++)
    		{
    			nEntityRecv += entitiesPerRankRecv[iProc];

    			sdispls.push_back( iProc == 0 ? 0 : entitiesPerRankSend[iProc - 1] - sdispls[iProc - 1] );
    			rdispls.push_back( iProc == 0 ? 0 : entitiesPerRankRecv[iProc - 1] - rdispls[iProc - 1] );

    			for (int iElem = 0; iElem < entitiesPerRankSend[iProc]; iElem++)
    			{
    				sendEntityGlobalIndex.push_back(neighborEntityGlobalIndex[iProc][iElem]);
    				sendEntityDataSize.push_back(neighborEntityData[iProc][iElem].size());
    			}
    		}

    		recvEntityGlobalIndex.resize(nEntityRecv);
    		recvEntityDataSize.resize(nEntityRecv);

    		MPI_Alltoallv (sendEntityGlobalIndex.data(), entitiesPerRankSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvEntityGlobalIndex.data()), entitiesPerRankRecv.data(), rdispls.data(), MPI_INT, comm );
    		MPI_Alltoallv (sendEntityDataSize.data(), entitiesPerRankSend.data(), sdispls.data(), MPI_INT, reinterpret_cast<int*>(recvEntityDataSize.data()), entitiesPerRankRecv.data(), rdispls.data(), MPI_INT, comm );


    		// 2.4) Serialize all data and communicate serialization lengths

    		// 2.5) Communicate package all data for each entity - send as integer pointers.

    		// 2.6) Unserialize data

    		// 3) Refill
    		// 3.1) Fill the scattermessagebuffer
    		// 3.2) Fill the recepient list ()







    	}


    };

  } // namespace CurvGrid

} // namespace Dune

#endif // #ifndef DUNE_CURVGRID_DECLARATION_HH
