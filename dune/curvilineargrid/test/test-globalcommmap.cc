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

#include <dune/curvilineargrid/utility/globalcommmap.hh>



using namespace Dune;

using namespace Dune::CurvGrid;



int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    static MPIHelper &mpihelper=Dune::MPIHelper::instance(argc,argv);

    typedef std::pair<int, int> Key;
    typedef std::pair<Key, int> DataPair;
    typedef GlobalCommMap<Key, int> GlobalMap;
    typedef typename GlobalMap::DataMapConstIter  MapConstIter;

    std::vector<DataPair> vec;
    for (int i = 0; i < 5; i++) {
    	vec.push_back(DataPair(Key(mpihelper.rank(), i), i*i + mpihelper.rank()));
    }

    GlobalMap globalmap;
    globalmap.init(mpihelper, vec);

    for (MapConstIter iter = globalmap.map().begin(); iter != globalmap.map().end(); iter++) {
    	std::cout << "On rank " << mpihelper.rank() << " Key " << iter->first.first << "," << iter->first.second << " Data " << iter->second << std::endl;
    }

    /** \brief leave program peacefully */
    return(0);
}
