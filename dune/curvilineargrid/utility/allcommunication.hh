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


namespace Dune {

namespace CurvGrid {



class AllCommunication
{

public:

	AllCommunication(MPIHelper & mpihelper) :
		mpihelper_(mpihelper)
	{
		rank_ = mpihelper_.rank();
		size_ = mpihelper_.size();
	}

	// A wrapper for MPI_Alltoallv communication
	// Each process sends a vector of T=POD to each other process
	// NOTE: IT IS ASSUMED THAT SUFFICIENT MEMORY FOR OUT AND LENGTHOUT IS RESERVED
	template <typename T>
	void communicate(
		const T * in,
		const int * lengthIn,
		T * out,
		int * lengthOut
	)
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// 1) Communicate lengths of the arrays to be communicated
		communicate_lengths(lengthIn, lengthOut, comm);

		// 2) Assemble and communicate data vectors
		communicate_data(in, lengthIn, out, lengthOut, comm);
	}


	// NOTE: out and lengthOut need not be reserved before init
	template <typename T>
	void communicate(
		const std::vector<T> & in,
		const std::vector<int> & lengthIn,
		std::vector<T> & out,
		std::vector<int> & lengthOut
	)
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// 1) Communicate lengths of the arrays to be communicated
		lengthOut.resize(size_);
		communicate_lengths(lengthIn.data(), reinterpret_cast<int *>(lengthOut.data()), comm);


		// 2) Assemble and communicate data vectors
		int lengthDataOut = 0;
		for (int iProc = 0; iProc < size_; iProc++)  { lengthDataOut += lengthOut[iProc]; }

		out.resize(lengthDataOut);
		communicate_data(in.data(), lengthIn.data(), reinterpret_cast<T *>(out.data()), lengthOut.data(), comm);
	}


protected:
	void communicate_lengths(const int * lengthIn, int * lengthOut, MPI_Comm comm)
	{
		MPI_Alltoall(lengthIn, 1, MPI_INT, lengthOut, 1, MPI_INT, comm);
	}


	template<typename T>
	void communicate_data(
		const T * in,
		const int * lengthIn,
		T * out,
		int * lengthOut,
		MPI_Comm comm
	)
	{
		int dataSize = sizeof(T);
		int lengthByteIn[size_], displByteIn[size_];
		int lengthByteOut[size_], displByteOut[size_];

		for (int iProc = 0; iProc < size_; iProc++)
		{
			lengthByteIn[iProc]  = dataSize * lengthIn[iProc];
			lengthByteOut[iProc] = dataSize * lengthOut[iProc];
			displByteIn[iProc]  = (iProc == 0) ? 0 : displByteIn[iProc - 1]  + lengthByteIn[iProc - 1];
			displByteOut[iProc] = (iProc == 0) ? 0 : displByteOut[iProc - 1] + lengthByteOut[iProc - 1];
		}

		MPI_Alltoallv (in, lengthByteIn, displByteIn, MPI_BYTE, out, lengthByteOut, displByteOut, MPI_BYTE, comm );
	}




private:

	template<class T>
	void print_array(T * arr, int length, std::string name)
	{
		std::cout << "process_" << rank_ << "printing array " << name << "= (";
		for (int i = 0; i < length; i++) { std::cout << arr[i] << " "; }
		std::cout << ")" << std::endl;
	}


	MPIHelper & mpihelper_;
	int rank_;
	int size_;


};

} // namespace CurvGrid

} // namespace Dune

#endif  // DUNE_ALLCOMMUNICATION_HH
