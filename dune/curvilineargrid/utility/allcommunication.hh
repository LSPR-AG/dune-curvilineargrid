#ifndef DUNE_ALLCOMMUNICATION_HH
#define DUNE_ALLCOMMUNICATION_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <dune/curvilineargrid/common/loggingmessage.hh>


namespace Dune {

namespace CurvGrid {



class AllCommunication
{


public:

	// Public Typedefs
	// *****************************************************


	// Implementation
	// *****************************************************

	AllCommunication(MPIHelper & mpihelper) :
		mpihelper_(mpihelper)
	{
		rank_ = mpihelper_.rank();
		size_ = mpihelper_.size();
	}

	// A wrapper for MPI_Alltoallv communication
	// Each process sends a vector of T=POD to each other process
	template <typename T>
	void communicate(
		std::vector<T> & in,
		std::vector<int> & lengthIn,
		std::vector<T> & out,
		std::vector<int> & lengthOut
	)
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		assert(in.size()  == size_);
		assert(lengthIn.size()  == size_);


		// 1) Communicate lengths of the arrays to be communicated
		lengthOut.resize(size_);
		MPI_Alltoall(lengthIn.data(), size_, MPI_INT, reinterpret_cast<int*>(lengthOut.data()), size_, MPI_INT, comm);


		// 2) Assemble and communicate data vectors
		int dataSize = sizeof(T);
		int lengthDataOut = 0;
		std::vector<int> lengthByteIn(size_),  displByteIn(size_);
		std::vector<int> lengthByteOut(size_), displByteOut(size_);

		for (int iProc = 0; iProc < size_; iProc++)
		{
			lengthByteIn[iProc]  = dataSize * lengthIn[iProc];
			lengthByteOut[iProc] = dataSize * lengthOut[iProc];
			displByteIn[iProc]  = (iProc == 0) ? 0 : displByteIn[iProc - 1]  + lengthByteIn[iProc - 1];
			displByteOut[iProc] = (iProc == 0) ? 0 : displByteOut[iProc - 1] + lengthByteOut[iProc - 1];
			lengthDataOut += lengthOut[iProc];
		}

		out.resize(lengthDataOut);
		MPI_Alltoallv (in.data(), lengthByteIn.data(), displByteIn.data(), MPI_BYTE,
				reinterpret_cast<int*>(out.data()), lengthByteOut.data(), displByteOut.data(), MPI_BYTE, comm );
	}



private:

	MPIHelper & mpihelper_;
	int rank_;
	int size_;


};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_ALLCOMMUNICATION_HH
