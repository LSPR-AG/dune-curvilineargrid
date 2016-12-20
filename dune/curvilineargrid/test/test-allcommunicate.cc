// Autotool generated header.
#include <config.h>

// Stl headers.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>


#include <dune/common/parallel/mpihelper.hh>

#include <dune/curvilineargrid/utility/allcommunication.hh>



/*************************************
 *  [TODO] Write assert-type tests to neighbour communication the same as it is done for all-to-all and gatherv
 *
 *
 */

using namespace Dune;

using namespace Dune::CurvGrid;



struct teststruct
{
	int rank_s;
	int rank_r;
	int data;

	teststruct() {}
	teststruct(int s, int r, int d) : rank_s(s), rank_r(r), data(d) {}
};




void test_gatherv_arr(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    AllCommunication allcomm(mpihelper);

    int lin = (rank + 1) * 2;
    int linTot = size * (size + 1);

    int lout[size];
    teststruct dataIn[lin];
    teststruct dataOut[linTot];
	for (int j = 0; j < lin; j++) { dataIn[j] = teststruct(rank, MPI_MASTER_RANK, j); }

    allcomm.gatherv(dataIn, lin, dataOut, lout);

    if (rank == MPI_MASTER_RANK) {
        int count = 0;
        for (int iProc = 0; iProc < size; iProc++) {
        	std::cout << "for proc " << iProc << " lout=" << lout[iProc] << std::endl;
            for (int i = 0; i < lout[iProc]; i++) {
            	std::cout
					<< "received from=" << dataOut[count].rank_s << " expected=" << iProc
            		<< "; sentto=" << dataOut[count].rank_r << " expected=" << MPI_MASTER_RANK
					<< "; data=" << dataOut[count].data << " expected=" << i
					<< std::endl;

            	assert(dataOut[count].rank_s == iProc);
            	assert(dataOut[count].rank_r == MPI_MASTER_RANK);
				assert(dataOut[count].data == i);
            	count++;
            }
        }
    }
}


void test_gatherv_vec(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    AllCommunication allcomm(mpihelper);

    int lin = (rank + 1) * 2;

    std::vector<int> lout;
    std::vector<teststruct> dataIn;
    std::vector<teststruct> dataOut;

    for (int j = 0; j < lin; j++) { dataIn.push_back(teststruct(rank, MPI_MASTER_RANK, j)); }

    allcomm.gatherv(dataIn, dataOut, lout);

    if (rank == MPI_MASTER_RANK) {
        int count = 0;
        for (int iProc = 0; iProc < size; iProc++) {
        	std::cout << "for proc " << iProc << " lout=" << lout[iProc] << std::endl;
            for (int i = 0; i < lout[iProc]; i++) {
            	std::cout
					<< "received from=" << dataOut[count].rank_s << " expected=" << iProc
            		<< "; sentto=" << dataOut[count].rank_r << " expected=" << MPI_MASTER_RANK
					<< "; data=" << dataOut[count].data << " expected=" << i
					<< std::endl;

            	assert(dataOut[count].rank_s == iProc);
            	assert(dataOut[count].rank_r == MPI_MASTER_RANK);
				assert(dataOut[count].data == i);
            	count++;
            }
        }
    }
}


void test_all2all_arr(MPIHelper & mpihelper) {
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    AllCommunication allcomm(mpihelper);

    std::cout << "process_" << rank << " ::: Testing pointer interface" << std::endl;

    int tmpsize = 5;

    int lin[size];
    int lout[size];

    teststruct dataIn[size * tmpsize];
    teststruct dataOut[size * tmpsize];

    for (int i = 0; i < size; i++) {
    	lin[i] = tmpsize;
    	for (int j = 0; j < tmpsize; j++) { dataIn[i * tmpsize + j] = teststruct(rank, i, j); }
    }

    allcomm.all2all(dataIn, lin, dataOut, lout);

    for (int iProc = 0; iProc < size; iProc++) {
    	for (int i = 0; i < tmpsize; i++) {
    		int iData = iProc * tmpsize + i;

        	std::cout
    			<< " from=" << dataOut[iData].rank_s << " expected=" << iProc
    			<< "; to=" << dataOut[iData].rank_r << " expected=" << rank
    			<< "; data=" << dataOut[iData].data << " expected=" << i << std::endl;

        	assert(dataOut[iData].rank_s == iProc);
        	assert(dataOut[iData].rank_r == rank);
			assert(dataOut[iData].data == i);
    	}
    }
}


void test_all2all_vec(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    AllCommunication allcomm(mpihelper);

    std::cout << "process_" << rank << " ::: Testing vector interface" << std::endl;

    int tmpsize = 5;

    std::vector<int> lin;
    std::vector<int> lout;

    std::vector<teststruct> dataIn(size * tmpsize);
    std::vector<teststruct> dataOut;

    for (int i = 0; i < size; i++) {
    	lin.push_back(tmpsize);
    	for (int j = 0; j < tmpsize; j++)  { dataIn[i * tmpsize + j] = teststruct(rank, i, j);  }
    }

    allcomm.all2all(dataIn, lin, dataOut, lout);


    for (int iProc = 0; iProc < size; iProc++) {
    	for (int i = 0; i < tmpsize; i++) {
    		int iData = iProc * tmpsize + i;

        	std::cout
    			<< " from=" << dataOut[iData].rank_s << " expected=" << iProc
    			<< "; to=" << dataOut[iData].rank_r << " expected=" << rank
    			<< "; data=" << dataOut[iData].data << " expected=" << i << std::endl;

        	assert(dataOut[iData].rank_s == iProc);
        	assert(dataOut[iData].rank_r == rank);
			assert(dataOut[iData].data == i);
    	}
    }
}


void test_all2allnb_arr(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    AllCommunication allcomm(mpihelper);

    std::cout << "process_" << rank << " ::: Testing neighbour pointer interface" << std::endl;

    // Communication sparsity pattern - communicate to next and previous (cyclic)
	int ranksIn[2];
	int ranksOut[2];
	ranksIn[0] = (rank == 0)      ? size - 1 : rank - 1;
	ranksIn[1] = (rank == size-1) ? 0        : rank + 1;

	int nNeighborIn = 2;
	int nNeighborOut;
    int tmpsize = 5;

    int lin[nNeighborIn];
    int lout[nNeighborIn];

    teststruct dataIn[nNeighborIn * tmpsize];
    teststruct dataOut[nNeighborIn * tmpsize];

    for (int i = 0; i < nNeighborIn; i++)
    {
    	lin[i] = tmpsize;

    	for (int j = 0; j < tmpsize; j++)
    	{
    		int l = i * tmpsize + j;
    		dataIn[l].rank_s = rank;
    		dataIn[l].rank_r = ranksIn[i];
    		dataIn[l].data = j;
    	}
    }

    //std::cout << "process_" << rank << " start comm of array of sizeof(T)=" << sizeof(teststruct) << " sizeof(arr<T>)=" << sizeof(dataIn) << std::endl;

    allcomm.all2allnb(dataIn, nNeighborIn, ranksIn, lin, dataOut, nNeighborOut, ranksOut, lout);
    assert(nNeighborOut == 2);

    int iData = 0;
    for (int i = 0; i < nNeighborOut; i++)
    {
    	for (int j = 0; j < tmpsize; j++)
    	{
    		teststruct thisData = dataOut[iData++];
    		std::cout << "process_" << rank << " received from=" << ranksOut[i];
    		std::cout << " data={ from=" << thisData.rank_s << " to=" << thisData.rank_r << " data=" << thisData.data << "}" << std::endl;
    	}
    }
}


void test_all2allnb_vec(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    AllCommunication allcomm(mpihelper);

    std::cout << "process_" << rank << " ::: Testing neighbour vector interface" << std::endl;

    // Communication sparsity pattern - communicate to next and previous (cyclic)
	std::vector<int> ranksIn(2);
	std::vector<int> ranksOut;
	ranksIn[0] = (rank == 0)      ? size - 1 : rank - 1;
	ranksIn[1] = (rank == size-1) ? 0        : rank + 1;

	int nNeighborIn = 2;
    int tmpsize = 5;

    std::vector<int> lin;
    std::vector<int> lout;

    std::vector<teststruct> dataIn(nNeighborIn * tmpsize);
    std::vector<teststruct> dataOut;

    for (int i = 0; i < nNeighborIn; i++)
    {
    	lin.push_back(tmpsize);

    	for (int j = 0; j < tmpsize; j++)
    	{
    		int l = i * tmpsize + j;
    		dataIn[l].rank_s = rank;
    		dataIn[l].rank_r = ranksIn[i];
    		dataIn[l].data = j;
    	}
    }

    allcomm.all2allnb(dataIn, ranksIn, lin, dataOut, ranksOut, lout);

    int iData = 0;
    for (int i = 0; i < lout.size(); i++)
    {
    	for (int j = 0; j < lout[i]; j++)
    	{
    		teststruct thisData = dataOut[iData++];
    		std::cout << "process_" << rank << " received from=" << ranksOut[i];
    		std::cout << " data={ from=" << thisData.rank_s << " to=" << thisData.rank_r << " data=" << thisData.data << "}" << std::endl;
    	}
    }
}



int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);

    if (mpihelper.size() > 1)
    {
        test_all2all_arr(mpihelper);
        test_all2all_vec(mpihelper);

        test_all2allnb_arr(mpihelper);
        test_all2allnb_vec(mpihelper);

        test_gatherv_arr(mpihelper);
    } else
    {
    	std::cout << "Skipping Allcommunicate test, as it is designed for parallel case" << std::endl;
    }





    /** \brief leave program peacefully */
    return(0);
}
