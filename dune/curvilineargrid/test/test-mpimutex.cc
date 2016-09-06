#include <iostream>
#include <fstream>
#include <dune/curvilineargrid/common/mpimutex.hh>


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    {
        MPIMutex m(rank, size, 0);
        m.lock();
        std::fstream outfile;
        outfile.open("testfile.txt", std::fstream::out | std::fstream::app);
        
        for (int i = 0; i < 100000; i++)  {
            outfile << rank << " ";
        }
        
        outfile << std::endl;
        outfile.close();
        m.unlock();
    }
    
    MPI_Finalize();
    
    return 0;
}
