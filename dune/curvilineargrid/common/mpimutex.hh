#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>


#ifndef DUNE_MPI_MUTEX_HH_
#define DUNE_MPI_MUTEX_HH_

#define MPI_MUTEX_MSG_TAG_BASE 1023

namespace Dune {

class MPIMutex {
    
public:

    MPIMutex(int rank, int size, int home) :
        rank_(rank),
        size_(size),
        home_(home),
        comm_(MPI_COMM_WORLD)
    {
        //std::cout << "Started Constructor" << std::endl;
        
        // What the hell does this do?
        static int tag = MPI_MUTEX_MSG_TAG_BASE;    

        tag_ = tag++;

        if (rank == home) {
            // Allocate and expose waitlist
            MPI_Alloc_mem(size_, MPI_INFO_NULL, &waitlist_);
            if (!waitlist_) {
                fprintf(stderr, "Warning: MPI_Alloc_mem failed on worker %2d\n", rank_);
                exit(1);
            }
            memset(waitlist_, 0, size_);
            MPI_Win_create(waitlist_, size_, 1, MPI_INFO_NULL, comm_, &win_);
        } else {
            // Don't expose anything
            waitlist_ = NULL;
            MPI_Win_create(waitlist_, 0, 1, MPI_INFO_NULL, comm_, &win_);
        }

        // It is essential that all processes finish the constructor before the mutex is used
        MPI_Barrier(MPI_COMM_WORLD);
        
        //std::cout << "Finished Constructor" << std::endl;
    }
    
    
    ~MPIMutex() {
        //std::cout << "Started Destructor" << std::endl;
        // Is it necessary to call barrier on destruction?
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank_ == home_) {
            // Free waitlist
            MPI_Win_free(&win_);
            MPI_Free_mem(waitlist_);
        } else {
            MPI_Win_free(&win_);
            assert(waitlist_ == NULL);
        }
        
        //std::cout << "Finished Destructor" << std::endl;
    }
    
    
    // Assign lock to this process
    // If another process has the lock, wait
    int lock() {
        //std::cout << "Started Lock" << std::endl;
        
        unsigned char waitlist[size_]; 
        unsigned char lockVar = 1;
        int i;
        
        // Try to acquire lock in one access epoch
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, home_, 0, win_);
        MPI_Put(&lockVar, 1, MPI_CHAR, home_, rank_ /* &win_[rank_] */, 1, MPI_CHAR, win_);
        MPI_Get(waitlist, size_, MPI_CHAR, home_, 0, size_, MPI_CHAR, win_);
        MPI_Win_unlock(home_, win_);

        assert(waitlist[rank_] == 1);

        // Count the 1's
        for (i = 0; i < size_; i++) {
            if (waitlist[i] == 1 && i != rank_) {
                // We have to wait for the lock
                // Dummy receive, no payload
                //printf("Worker %d waits for lock\n", rank_);
                MPI_Recv(&lockVar, 0, MPI_CHAR, MPI_ANY_SOURCE, tag_, comm_, MPI_STATUS_IGNORE);
                break;
            }
        }
        //printf("Worker %d has the lock\n", rank_);
        
        //std::cout << "Finished Lock" << std::endl;
        
        return 0;
    }
    
    
    // Check if anybody has the lock at the moment
    int trylock() {
        //std::cout << "Started TryLock" << std::endl;
        
        unsigned char waitlist[size_]; 
        unsigned char lockVar = 1;
        int i;

        // Try to acquire lock in one access epoch
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, home_, 0, win_);
        MPI_Put(&lockVar, 1, MPI_CHAR, home_, rank_ /* &win[rank_] */, 1, MPI_CHAR, win_);
        MPI_Get(waitlist, size_, MPI_CHAR, home_, 0, size_, MPI_CHAR, win_);
        MPI_Win_unlock(home_, win_);

        assert(waitlist[rank_] == 1);

        // Count the 1's
        for (i = 0; i < size_; i++) {
            if (waitlist[i] == 1 && i != rank_) {
                //Lock is already held, return immediately
                return 1;
            }
        }
        //printf("Worker %d has the lock\n", rank_);
        
        //std::cout << "Finished TryLock" << std::endl;
        return 0;
    }
    
    
    // Remove lock from this process. After, any process can lock
    int unlock() {
        //std::cout << "Started Unlock" << std::endl;
        
        unsigned char waitlist[size_]; 
        unsigned char lockVar = 0;
        int i, next;

        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, home_, 0, win_);
        MPI_Put(&lockVar, 1, MPI_CHAR, home_, rank_ /* &win[rank_] */, 1, MPI_CHAR, win_);
        MPI_Get(waitlist, size_, MPI_CHAR, home_, 0, size_, MPI_CHAR, win_);
        MPI_Win_unlock(home_, win_);
        
        assert(waitlist[rank_] == 0);

        // If there are other processes waiting for the lock, transfer ownership
        next = (rank_ + 1 + size_) % size_;
        for (i = 0; i < size_; i++, next = (next + 1) % size_) {
            if (waitlist[next] == 1) {
                // Dummy send, no payload
                //printf("Worker %d transfers lock ownership to worker %d\n", rank_, i);
                MPI_Send(&lockVar, 0, MPI_CHAR, next, tag_, comm_);
                break;
            }
        }
        
        //std::cout << "Finished Unlock" << std::endl;
        return 0;
    }
    
    
    int home() { return home_; }
    
    int tag() { return tag_; }
    
private:

    int rank_;
    int size_;
    int home_;
    int tag_;
    MPI_Comm comm_;
    MPI_Win win_;
    unsigned char *waitlist_;
};


}


#endif //DUNE_MPI_MUTEX_HH_
