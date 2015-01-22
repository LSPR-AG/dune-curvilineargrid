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
    	std::cout << "process_" << rank << " received-struct " << dataOut[i].rank_s << " " << dataOut[i].rank_r << " " << dataOut[i].data << std::endl;
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

    	for (int j = 0; j < 5; j++)
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
    	std::cout << "process_" << rank << " received-struct " << dataOut[i].rank_s << " " << dataOut[i].rank_r << " " << dataOut[i].data << std::endl;
    }
}


int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);


    //test_pointerinterface(mpihelper);
    test_vectorinterface(mpihelper);




    /** \brief leave program peacefully */
    return(0);
}
