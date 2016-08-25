#ifndef DUNE_PARALLEL_DATA_WRITER_HH_
#define DUNE_PARALLEL_DATA_WRITER_HH_


#include <vector>
#include <string>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>





namespace Dune
{


// This class writes parallel vectors to file by gathering them on MPI_MASTER_RANK
// [TODO] This code could be optimized by using MPI_File_Write_Ordered
template<class Grid, class IndexType, class DataType>
class ParallelDataWriter
{
	typedef std::vector<IndexType>     IndexVector;
	typedef std::vector<DataType>      DataVector;
	typedef std::vector<int>  SizeVector;

	static const int MPI_MASTER_RANK = 0;
	static const int ELEMENT_CODIM = 0;


protected:


	struct ElementHolder
	{
		IndexType  globalIndex_;
		DataVector dof_;

		bool operator < (const ElementHolder & other) const
		{
			return (globalIndex_ < other.globalIndex_);
		}
	};

public:

	// Sortable index given by array of integers of arbitrary length
	struct FancyIndex
	{
		std::vector<int> index;

		// Particularly useful constructor for index pairs
		FancyIndex(int a, int b)  {
			index = std::vector<int> {a, b};
		}


		bool operator < (const FancyIndex & other) const
		{
			assert(index.size() == other.index.size());
			for (int i = 0; i < index.size(); i++)  { if (index[i] != other.index[i])  { return index[i] < other.index[i]; }  }
			return true;  // If both are equal, it is sort of irrelevant what we return
		}

	};



	static void writeParallelData2File(
	  std::string filename,
	  IndexVector & interiorElementGlobalIndex,
	  SizeVector  & interiorElementNDof,
	  DataVector  & data,
	  const Grid & grid)
	{
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "Started writing the parallel data to file [[...");

		int rank = grid.comm().rank();
		int size = grid.comm().size();

		// ******************************************************
		// Self-test of consistency of input data
		// ******************************************************
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Testing");

		int nElement = grid.numInternal(ELEMENT_CODIM);
		int nData    = data.size();

		assert(interiorElementGlobalIndex.size() == nElement);  // Both index and nDof arrays must be of size nElement
		assert(interiorElementNDof.size()        == nElement);

		int nDataTest = 0;
		for (int i = 0; i < nElement; i++)  { nDataTest += interiorElementNDof[i]; }
		assert(nDataTest == nData);  // The number of DoF on this process should be the sum of Dof of all its interior elements


		// ******************************************************
		// Compute total number of elements and data points
		// ******************************************************
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Computing total data size");

		int nElementTot = grid.comm().sum(nElement);
		int nDataTot    = grid.comm().sum(nData);

		// ******************************************************
		// Communicate data sizes over all processes to Master process
		// ******************************************************
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Communicating sizes to master process");

		SizeVector nElementProcess(size, 0);
		SizeVector nDataProcess(size, 0);
		grid.comm().gather (&nElement, static_cast<int *>(nElementProcess.data()), 1, MPI_MASTER_RANK);
		grid.comm().gather (&nData,    static_cast<int *>(nDataProcess.data()),    1, MPI_MASTER_RANK);


		// ******************************************************
		// Compute data displacements
		// ******************************************************
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Computing displacements");

		SizeVector displElementProcess(1, 0);
		SizeVector displDataProcess(1, 0);

		for (int i = 0; i < size - 1; i++)  {
			displElementProcess.push_back(displElementProcess[i] + nElementProcess[i]);
			displDataProcess.push_back(displDataProcess[i] + nDataProcess[i]);
		}


		// ******************************************************
		// Gather all data on master process
		// ******************************************************
		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Communicating data to master process");

		IndexVector  processElementGlobalIndex(nElementTot, 0);
		SizeVector   processElementNDof       (nElementTot, 0);
		DataVector   processData              (nDataTot,    DataType(0));

		grid.comm().gatherv (interiorElementGlobalIndex.data(), nElement, static_cast<IndexType *> (processElementGlobalIndex.data()), nElementProcess.data(), displElementProcess.data(), MPI_MASTER_RANK);
		grid.comm().gatherv (interiorElementNDof.data(),        nElement, static_cast<int *>       (processElementNDof.data()),        nElementProcess.data(), displElementProcess.data(), MPI_MASTER_RANK);
		grid.comm().gatherv (data.data(),                       nData,    static_cast<DataType *>  (processData.data()),               nDataProcess.data(),    displDataProcess.data(),    MPI_MASTER_RANK);


		// ******************************************************
		// Package all data into sortable structure
		// ******************************************************

		if (rank == MPI_MASTER_RANK)
		{
			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Packaging data to sortable structure");

			std::vector<ElementHolder> elemHolder(nElementTot);

			int dataIndex = 0;
			for (int i = 0; i < nElementTot; i++)
			{
				elemHolder[i].globalIndex_ = processElementGlobalIndex[i];
				elemHolder[i].dof_.clear();

				for (int j = 0; j < processElementNDof[i]; j++) {
					elemHolder[i].dof_.push_back(processData[dataIndex++]);
				}
			}

			assert(dataIndex == nDataTot);  // If everything went right, the data index should have spanned exactly all of the data


			// ******************************************************
			// Sort packaged data
			// ******************************************************
			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Sorting data");

			std::sort(elemHolder.begin(), elemHolder.end());


			// ******************************************************
			// Write structure to file
			// ******************************************************
			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "  -- Writing data to file");

			std::ofstream outfile;
			outfile.open(filename);

			for (int i = 0; i < nElementTot; i++)
			{
				for (int j = 0; j < elemHolder[i].dof_.size(); j++)
				{
					outfile << elemHolder[i].globalIndex_ << " " << j << " " << elemHolder[i].dof_[j] << std::endl;
				}
			}

			outfile.close();
		}

		LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, "Finished writing the parallel data to file ...]]");
	}




};


};

#endif //DUNE_PARALLEL_DATA_WRITER_HH_
