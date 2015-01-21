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

/** *\brief include functionality that encapsulates MPI */
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/curvilineargrid/utility/allcommunication.hh>



using namespace Dune;

using namespace CurvGrid;



struct teststruct
{
	int rank_s;
	int rank_r;
	int data;
};

template <class T>
std::string vector2string(const T & V)
{
    std::stringstream tmp_stream;

    int nEntry = V.size();
    if (nEntry == 0)  { tmp_stream << "Null"; }
    for (int i = 0; i < nEntry; i++) {
    	tmp_stream << V[i];
    	if (i != nEntry - 1) { tmp_stream << " "; }
    }
    return tmp_stream.str();
}



int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);

    Dune::CurvGrid::AllCommunication allcomm(mpihelper);


    std::vector<int> lin;
    std::vector<int> lout;

    std::vector<teststruct> dataIn;
    std::vector<teststruct> dataOut;


    for (int i = 0; i < mpihelper.size(); i++)
    {
    	lin.push_back(5);

    	for (int j = 0; j < 5; j++)
    	{
    		teststruct tmp;
    		tmp.rank_s = mpihelper.rank();
    		tmp.rank_r = i;
    		tmp.data = j;
    		dataIn.push_back(tmp);
    	}
    }

    allcomm.communicate(dataIn, lin, dataOut, lout);

    std::cout << "process_" << mpihelper.rank() << " datasize-send " << vector2string(lin) << std::endl;
    std::cout << "process_" << mpihelper.rank() << " datasize-recv " << vector2string(lout) << std::endl;

    for (int i = 0; i < dataOut.size(); i++)
    {
    	std::cout << "process_" << mpihelper.rank() << " received-struct " << dataOut[i].rank_s << " " << dataOut[i].rank_r << " " << dataOut[i].data << std::endl;
    }




    /** \brief leave program peacefully */
    return(0);
}
