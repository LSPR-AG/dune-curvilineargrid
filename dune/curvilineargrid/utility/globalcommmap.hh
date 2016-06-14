#ifndef DUNE_GLOBALCOMMMAP_HH
#define DUNE_GLOBALCOMMMAP_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>
#include <numeric>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>


namespace Dune {

namespace CurvGrid {


template <class Key, class Data>
class GlobalCommMap {
public:
	typedef std::map<Key, Data> DataMap;
	typedef std::pair<Key, Data> DataPair;
	typedef typename DataMap::const_iterator DataMapConstIter;

public:

	GlobalCommMap() {}

	void init(MPIHelper &mpihelper, const std::vector<DataPair> & thisSrc) {
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();
		int mpi_err;
		int rank = mpihelper.rank();
		int size = mpihelper.size();

		// Construct the number of communicated data
		int structSize = sizeof(DataPair);	// the size of the communiated structure
		int nCommThis = thisSrc.size();
		int sizeCommThis = nCommThis * structSize;
		std::vector<int> nCommRecv(size, 0);			// Number of elements that are sent by each of the processes
		std::vector<int> sizeCommRecv(size, 0);		// Total size in bytes sent by each process
		std::vector<int> displCommRecv(size, 0);		// Displacement in bytes
		mpi_err = MPI_Allgather(&nCommThis, 1, MPI_INT, reinterpret_cast<int *>(nCommRecv.data()), 1, MPI_INT, comm);
		int nCommRecvTot = std::accumulate(nCommRecv.begin(), nCommRecv.end(), 0);		// The total number of structrures that will be received by any one process
		for (int i = 0; i < size; i++) { sizeCommRecv[i] = nCommRecv[i] * structSize; }
		for (int i = 1; i < size; i++) { displCommRecv[i] = displCommRecv[i-1] + sizeCommRecv[i-1]; }

		// Communicate
		std::vector<DataPair> globalSrc(nCommRecvTot);
		mpi_err = MPI_Allgatherv(
				reinterpret_cast<const DataPair *>(thisSrc.data()),
				sizeCommThis,
				MPI_BYTE,
				reinterpret_cast<DataPair *>(globalSrc.data()),
				reinterpret_cast<int *>(sizeCommRecv.data()),
				reinterpret_cast<int *>(displCommRecv.data()),
				MPI_BYTE,
				comm
		);

		// Fill in the data
		for (int i = 0; i < nCommRecvTot; i++) {
			assert(globalmap_.find(globalSrc[i].first) == globalmap_.end());  // Make sure that there are no two data communicated with the same key
			globalmap_.insert(globalSrc[i]);
		}
	}

	const DataMap & map() { return globalmap_; }


private:

	DataMap globalmap_;
};










} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_GLOBALCOMMMAP_HH
