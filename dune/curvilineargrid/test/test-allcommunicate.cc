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



using namespace Dune;

using namespace CurvGrid;



struct teststruct
{
	int rank_s;
	int rank_r;
	int data;
};



void test_pointerinterface(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    Dune::CurvGrid::AllCommunication allcomm(mpihelper);

    std::cout << "process_" << rank << " ::: Testing pointer interface" << std::endl;

    int tmpsize = 5;

    int lin[size];
    int lout[size];

    teststruct dataIn[size * tmpsize];
    teststruct dataOut[size * tmpsize];

    for (int i = 0; i < size; i++)
    {
    	lin[i] = tmpsize;

    	for (int j = 0; j < tmpsize; j++)
    	{
    		int l = i * tmpsize + j;
    		dataIn[l].rank_s = rank;
    		dataIn[l].rank_r = i;
    		dataIn[l].data = j;
    	}
    }

    //std::cout << "process_" << rank << " start comm of array of sizeof(T)=" << sizeof(teststruct) << " sizeof(arr<T>)=" << sizeof(dataIn) << std::endl;

    allcomm.communicate(dataIn, lin, dataOut, lout);

    for (int i = 0; i < size * tmpsize; i++)
    {
    	std::cout << "process_" << rank << " received-struct from=" << dataOut[i].rank_s << " to=" << dataOut[i].rank_r << " data=" << dataOut[i].data << std::endl;
    }
}


void test_vectorinterface(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    Dune::CurvGrid::AllCommunication allcomm(mpihelper);

    std::cout << "process_" << rank << " ::: Testing vector interface" << std::endl;

    int tmpsize = 5;

    std::vector<int> lin;
    std::vector<int> lout;

    std::vector<teststruct> dataIn(size * tmpsize);
    std::vector<teststruct> dataOut;

    for (int i = 0; i < size; i++)
    {
    	lin.push_back(tmpsize);

    	for (int j = 0; j < tmpsize; j++)
    	{
    		int l = i * tmpsize + j;
    		dataIn[l].rank_s = rank;
    		dataIn[l].rank_r = i;
    		dataIn[l].data = j;
    	}
    }

    allcomm.communicate(dataIn, lin, dataOut, lout);


    for (int i = 0; i < dataOut.size(); i++)
    {
    	std::cout << "process_" << rank << " received-struct from=" << dataOut[i].rank_s << " to=" << dataOut[i].rank_r << " data=" << dataOut[i].data << std::endl;
    }
}


void test_neighbor_pointerinterface(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    Dune::CurvGrid::AllCommunication allcomm(mpihelper);

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

    allcomm.communicate_neighbors(dataIn, nNeighborIn, ranksIn, lin, dataOut, nNeighborOut, ranksOut, lout);
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


void test_neighbor_vectorinterface(MPIHelper & mpihelper)
{
    int rank = mpihelper.rank();
    int size = mpihelper.size();
    Dune::CurvGrid::AllCommunication allcomm(mpihelper);

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

    allcomm.communicate_neighbors(dataIn, ranksIn, lin, dataOut, ranksOut, lout);

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
        //test_pointerinterface(mpihelper);
        //test_vectorinterface(mpihelper);

        //test_neighbor_pointerinterface(mpihelper);
        test_neighbor_vectorinterface(mpihelper);
    } else
    {
    	std::cout << "Skipping Allcommunicate test, as it is designed for parallel case" << std::endl;
    }





    /** \brief leave program peacefully */
    return(0);
}
