#ifndef DUNE_ALLCOMMUNICATION_HH
#define DUNE_ALLCOMMUNICATION_HH

#include <limits>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>
#include <numeric>

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/curvilineargrid/common/constant.hh>


namespace Dune {

namespace CurvGrid {


/** \brief This class implements wrappers for MPI_Alltoallv communication, as well as sparse MPI_Neighbor_Alltoallv,
 * allowing the communication of abstract templated POD arrays among processes. Note that main restriction on the
 * templated datatype is that its sizeof() must be constant and completely determined at compile-time.
 *
 * */
class AllCommunication
{

public:

	AllCommunication(MPIHelper & mpihelper) :
		mpihelper_(mpihelper)
	{
		rank_ = mpihelper_.rank();
		size_ = mpihelper_.size();
	}




	/** \brief A wrapper for Gatherv communication. Non scalable - do not use on very large architectures.
	 * For sparse communication use gathev_neighbor
	 * Works optimal if each process has non-zero data to contribute to the Master process
	 *
	 *  \param[in] in				buffer with data to send to Master process
	 *  \param[in] lengthIn		array of sizes of data to be sent to each process (size = mpihelper.size())
	 *  \param[in] out				buffer with data to be received by Master process. For all other processes it is untouched
	 *
	 *  \note This operation does not reserve memory for output arrays. It is assumed that sufficient memory is already reserved
	 *
	 * */
	template <typename T>
	void gatherv(
			const T * in,
			int lengthIn,
			T * out,
			int * lengthOut) const
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// 1) Calculate the lengths to be sent/received
		MPI_Gather(&lengthIn, 1, MPI_INT, lengthOut, 1, MPI_INT, MPI_MASTER_RANK, comm);

//		if (rank_ == MPI_MASTER_RANK) {
//			std::stringstream ss;
//			for (int i = 0; i < size_; i++) { ss << lengthOut[i] << " "; }
//			std::cout << "on " << rank_ << " stuff " << ss.str() << std::endl;
//		}

		// 2) Gather data
		gatherv_comm_data(in, lengthIn, out, lengthOut, comm);
	}


	/** \brief A wrapper for Gatherv communication. Non scalable - do not use on very large architectures.
	 * For sparse communication use gathev_neighbor
	 * Works optimal if each process has non-zero data to contribute to the Master process
	 *
	 *  \param[in] in				buffer with data to send to Master process
	 *  \param[in] lengthIn		array of sizes of data to be sent to each process (size = mpihelper.size())
	 *  \param[in] out				buffer with data to be received by Master process. For all other processes it is untouched
	 *
	 *  \note This operation does not reserve memory for output arrays. It is assumed that sufficient memory is already reserved
	 *
	 * */
	template <typename T>
	void gatherv(
			const std::vector<T> & in,
			std::vector<T> & out,
			std::vector<int> & lengthOut) const
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// 1) Calculate the lengths to be sent/received
		int lengthIn = in.size();
		if (rank_ == MPI_MASTER_RANK) { lengthOut.resize(size_); }
		MPI_Gather(&lengthIn, 1, MPI_INT, reinterpret_cast<int *>(lengthOut.data()), 1, MPI_INT, MPI_MASTER_RANK, comm);

		// 2) Assemble and communicate data vectors
		if (rank_ == MPI_MASTER_RANK) {
			int lengthTot = std::accumulate(lengthOut.begin(), lengthOut.end(), 0);
			out.resize(lengthTot);
		}

		gatherv_comm_data(in.data(), lengthIn, reinterpret_cast<T *>(out.data()),  lengthOut.data(), comm);
	}

	/** \brief A wrapper for MPI_Alltoallv communication. Non scalable - do not use on very large architectures.
	 * Works optimal if each process has non-zero communication to each other. For sparse communication use
	 * communicate_neighbor
	 *
	 *  \param[in] in            buffer with data to send to neighbouring processes
	 *  \param[in] lengthIn      array of sizes of data to be sent to each in-neighbour (size = mpihelper.size())
	 *  \param[in] out           buffer with data to be received by this process
	 *  \param[in] lengthOut     array of sizes of data to be received from each out-neighbor (size = mpihelper.size())
	 *
	 *  \note This operation does not reserve memory for output arrays. It is assumed that sufficient memory is already reserved
	 *
	 * */
	template <typename T>
	void all2all(
		const T * in,
		const int * lengthIn,
		T * out,
		int * lengthOut) const
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// 1) Communicate lengths of the arrays to be communicated
		all2all_comm_len(lengthIn, lengthOut, false, comm);

		// 2) Assemble and communicate data vectors
		all2all_comm_data(in, lengthIn, out, lengthOut, size_, comm);
	}


	/** \brief A wrapper for MPI_Alltoallv communication. Non scalable - do not use on very large architectures.
	 * Works optimal if each process has non-zero communication to each other. For sparse communication use
	 * communicate_neighbor
	 *
	 *  \param[in] in            buffer with data to send to neighbouring processes
	 *  \param[in] lengthIn      vector of sizes of data to be sent to each in-neighbour (size = mpihelper.size())
	 *  \param[in] out           buffer with data to be received by this process
	 *  \param[in] lengthOut     vector of sizes of data to be received from each out-neighbor (size = mpihelper.size())
	 *
	 *  \note This operation automatically reserves memory for all output vectors, so no a priori knowledge
	 *  of output sizes required. It is assumed that out and lengthOut are empty vectors.
	 *
	 * */
	template <typename T>
	void all2all(
		const std::vector<T> & in,
		const std::vector<int> & lengthIn,
		std::vector<T> & out,
		std::vector<int> & lengthOut) const
	{
		MPI_Comm comm = Dune::MPIHelper::getCommunicator();

		// 1) Communicate lengths of the arrays to be communicated
		lengthOut.resize(size_);
		all2all_comm_len(lengthIn.data(), reinterpret_cast<int *>(lengthOut.data()), false, comm);


		// 2) Assemble and communicate data vectors
		int lengthDataOut = 0;
		for (int iProc = 0; iProc < size_; iProc++)  { lengthDataOut += lengthOut[iProc]; }

		out.resize(lengthDataOut);
		all2all_comm_data(in.data(), lengthIn.data(), reinterpret_cast<T *>(out.data()), lengthOut.data(), size_, comm);
	}


	/** \brief Function allowing scalable MPI-Alltoallv communication
	 *  \param[in] in            buffer with data to send to neighbouring processes
	 *  \param[in] nNeighborIn   number of neighbour processes this process sends to (in-neighbours)
	 *  \param[in] ranksIn       array of ranks of in-neighbours of this process
	 *  \param[in] lengthIn      array of sizes of data to be sent to each in-neighbour
	 *  \param[in] out           buffer with data to be received by this process
	 *  \param[in] nNeighborOut  number of processes that will send something to this process (out-neighbours)
	 *  \param[in] ranksOut      array of ranks of out-neighbours of this process
	 *  \param[in] lengthOut     array of sizes of data to be received from each out-neighbor
	 *
	 * */
	template <typename T>
	void all2allnb(
		const T * in,
		int nNeighborIn,
		const int * ranksIn,
		const int * lengthIn,
		T * out,
		int & nNeighborOut,
		int * ranksOut,
		int * lengthOut) const
	{
		// 0) Construct neighbor communicator
		// *******************************************************************

		MPI_Comm comm = Dune::MPIHelper::getCommunicator();
		MPI_Comm comm_neighbor;

		const int ranksOutInput[1] = {rank_};
		std::vector<int> weightsIn(nNeighborIn, 1);

		MPI_Dist_graph_create(
			comm, 1, ranksOutInput, &nNeighborIn, ranksIn, weightsIn.data(), MPI_INFO_NULL, true, &comm_neighbor);

		/*
		MPI_Dist_graph_create_adjacent(
			comm,
		    nNeighborOut, ranksOut, reinterpret_cast<const int*>(weightsOut.data()),
		    nNeighborIn, ranksIn, reinterpret_cast<const int*>(weightsIn.data()),
		    MPI_INFO_NULL, true, comm_neighbor);
		*/

		// Find out-neighbour length
		int weighted = (int) true;
		int nNeighborInTmp;
		MPI_Dist_graph_neighbors_count(comm_neighbor, &nNeighborInTmp, &nNeighborOut, &weighted);
		assert(nNeighborInTmp == nNeighborIn);

		// Find out-neighbour ranks
		int neighborRanksInTmp[nNeighborIn];
		int neighborWeightsInTmp[nNeighborIn];
		int neighborWeightsOutTmp[nNeighborOut];
		MPI_Dist_graph_neighbors(comm_neighbor, nNeighborOut, ranksOut, neighborWeightsOutTmp, nNeighborIn, neighborRanksInTmp, neighborWeightsInTmp);


		// 1) Communicate lengths of the arrays to be communicated
		// *******************************************************************
		all2all_comm_len(lengthIn, lengthOut, true, comm_neighbor);


		// 2) Assemble and communicate data vectors
		// *******************************************************************
		all2all_comm_data(in, lengthIn, out, lengthOut, nNeighborOut, comm_neighbor);
	}


	/** \brief Function allowing scalable MPI-Alltoallv communication [vector form]
	 *  \param[in] in            buffer with data to send to neighbouring processes
	 *  \param[in] ranksIn       vector of ranks of in-neighbours of this process
	 *  \param[in] lengthIn      vector of sizes of data to be sent to each in-neighbour
	 *  \param[in] out           buffer with data to be received by this process
	 *  \param[in] ranksOut      array of ranks of out-neighbours of this process
	 *  \param[in] lengthOut     vector of sizes of data to be received from each out-neighbor
	 *
	 *  \note space for all output vectors is reserved automatically
	 *
	 * */
	template <typename T>
	void all2allnb(
		const std::vector<T>   & in,
		const std::vector<int> & ranksIn,
		const std::vector<int> & lengthIn,
		std::vector<T> & out,
		std::vector<int> & ranksOut,
		std::vector<int> & lengthOut) const
	{
		// 0) Construct neighbor communicator
		// *******************************************************************

		MPI_Comm comm = Dune::MPIHelper::getCommunicator();
		MPI_Comm comm_neighbor;

		int nNeighborIn = ranksIn.size();
		int nNeighborOut;
		const int ranksOutInput[1] = {rank_};
		std::vector<int> weightsIn(nNeighborIn, 1);

		MPI_Dist_graph_create(
			comm, 1, ranksOutInput, &nNeighborIn, ranksIn.data(), weightsIn.data(), MPI_INFO_NULL, true, &comm_neighbor);

		/*
		MPI_Dist_graph_create_adjacent(
			comm,
		    nNeighborOut, ranksOut, reinterpret_cast<const int*>(weightsOut.data()),
		    nNeighborIn, ranksIn, reinterpret_cast<const int*>(weightsIn.data()),
		    MPI_INFO_NULL, true, comm_neighbor);
		*/

		// Find out-neighbour length
		int weighted = (int) true;
		int nNeighborInTmp;
		MPI_Dist_graph_neighbors_count(comm_neighbor, &nNeighborInTmp, &nNeighborOut, &weighted);
		assert(nNeighborInTmp == nNeighborIn);

		ranksOut.resize(nNeighborOut);

		// Find out-neighbour ranks
		int neighborRanksInTmp[nNeighborIn];
		int neighborWeightsInTmp[nNeighborIn];
		int neighborWeightsOutTmp[nNeighborOut];
		MPI_Dist_graph_neighbors(comm_neighbor, nNeighborOut,
				reinterpret_cast<int *>(ranksOut.data()),
				neighborWeightsOutTmp,
				nNeighborIn,
				neighborRanksInTmp,
				neighborWeightsInTmp);


		// 1) Communicate lengths of the arrays to be communicated
		// *******************************************************************
		lengthOut.resize(nNeighborOut);
		all2all_comm_len(lengthIn.data(), reinterpret_cast<int *>(lengthOut.data()), true, comm_neighbor);


		// 2) Assemble and communicate data vectors
		// *******************************************************************
		int lengthDataOut = 0;
		for (int iProc = 0; iProc < nNeighborOut; iProc++)  { lengthDataOut += lengthOut[iProc]; }

		out.resize(lengthDataOut);
		all2all_comm_data(in.data(), lengthIn.data(), reinterpret_cast<T *>(out.data()), lengthOut.data(), nNeighborOut, comm_neighbor);
	}




protected:

	template<typename T>
	void gatherv_comm_data(const T * in, int lengthIn, T * out, const int * lengthOut, MPI_Comm comm) const {

		int dataSize = sizeof(T);
		int lengthByteIn = lengthIn * dataSize;

		if (rank_ == MPI_MASTER_RANK) {
			int lengthByteOut[size_], displByteOut[size_];
			for (int iProc = 0; iProc < size_; iProc++) {
				lengthByteOut[iProc] = dataSize * lengthOut[iProc];
				displByteOut[iProc] = (iProc == 0) ? 0 : displByteOut[iProc - 1] + lengthByteOut[iProc - 1];
			}

			MPI_Gatherv(in, lengthByteIn, MPI_BYTE, out, lengthByteOut, displByteOut, MPI_BYTE, MPI_MASTER_RANK, comm);
		} else {
			MPI_Gatherv(in, lengthByteIn, MPI_BYTE, out, nullptr, nullptr, MPI_BYTE, MPI_MASTER_RANK, comm);
		}
	}


	void all2all_comm_len(const int * lengthIn, int * lengthOut, bool neighbour, MPI_Comm comm) const
	{
		if (neighbour)  { MPI_Neighbor_alltoall (lengthIn, 1, MPI_INT, lengthOut, 1, MPI_INT, comm); }
		else            { MPI_Alltoall          (lengthIn, 1, MPI_INT, lengthOut, 1, MPI_INT, comm); }

	}


	template<typename T>
	void all2all_comm_data(
		const T * in,
		const int * lengthIn,
		T * out,
		int * lengthOut,
		int nNeighborOut,
		MPI_Comm comm) const
	{
		int dataSize = sizeof(T);
		int lengthByteIn[nNeighborOut], displByteIn[nNeighborOut];
		int lengthByteOut[nNeighborOut], displByteOut[nNeighborOut];

		for (int iProc = 0; iProc < nNeighborOut; iProc++)
		{
			lengthByteIn[iProc]  = dataSize * lengthIn[iProc];
			lengthByteOut[iProc] = dataSize * lengthOut[iProc];
			displByteIn[iProc]  = (iProc == 0) ? 0 : displByteIn[iProc - 1]  + lengthByteIn[iProc - 1];
			//displByteOut[iProc] = (iProc == 0) ? 0 : displByteOut[iProc - 1] + lengthByteOut[iProc - 1];
		}

		for (int iProc = 0; iProc < nNeighborOut; iProc++)
		{
			displByteOut[iProc] = (iProc == 0) ? 0 : displByteOut[iProc - 1] + lengthByteOut[iProc - 1];
		}

		if (nNeighborOut == size_)  { MPI_Alltoallv          (in, lengthByteIn, displByteIn, MPI_BYTE, out, lengthByteOut, displByteOut, MPI_BYTE, comm ); }
		else                        { MPI_Neighbor_alltoallv (in, lengthByteIn, displByteIn, MPI_BYTE, out, lengthByteOut, displByteOut, MPI_BYTE, comm ); }
	}




private:

	template<class T>
	void print_array(T * arr, int length, std::string name) const
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
