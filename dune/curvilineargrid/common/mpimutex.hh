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

namespace CurvGrid {


// [FIXME] Currently uses CHAR for waitlist. ATM, this code only scales to 256 processors
class MPIMutex {

	static const char MPI_MUTEX_UNLOCKED = 0;
	static const char MPI_MUTEX_LOCKED = 1;
    
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
            // Allocate and expose MPI Memory. This is standard CHAR* memory, optimized for Remote Memory Access (RMA)
            MPI_Alloc_mem(size_, MPI_INFO_NULL, &waitlist_);
            if (!waitlist_) { fprintf(stderr, "Warning: MPI_Alloc_mem failed on worker %2d\n", rank_); exit(1); }

            // Initialize all entries of waitlist as unlocked
            memset(waitlist_, MPI_MUTEX_UNLOCKED, size_);

            // Create the Remote Memory Access window for waitlist
            MPI_Win_create(waitlist_, size_, 1, MPI_INFO_NULL, comm_, &win_);
        } else {
            // Don't expose anything - the waitlist is stored only on one process
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
    /******
     * Algorithm: data
     * 1) Exists shared array waitlist_[size_], which has 1's for every process currently having lock or waiting for lock, and 0's otherwise
     *
     * Algorithm: lock()
     * 1) When a process locks, it writes its intention to lock into waitlist_, then gets copy of waitlist
     * 2) The received local copy of waitlist has 1's for all processes that attempted to lock before this process
     * 3) This process waits to receive an (empty) note from each process that is before it. Then it is considered to have the lock
     *
     * Algorithm: unlock()
     * 1) When a process unlocks(), it removes itself from waitlist_ by PUT 0, then gets copy of waitlist
     * 2) The received local copy of waitlist has 1's for all processes that attempted to lock after this process, so far
     * 3) This process waits to send an (empty) note to each process that is after it at this moment. Then it is considered to have unlocked
     *
     *
     *
     */
    int lock() {
        //std::cout << "Started Lock" << std::endl;
        
        unsigned char waitlist[size_]; 
        unsigned char lockVar = MPI_MUTEX_LOCKED;
        int i;

        // Lock the waitlist_ while doing operations on it, so there are no collisions
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, home_, 0, win_);								// Try to acquire lock in one access epoch
        MPI_Put(&lockVar, 1, MPI_CHAR, home_, rank_, 1, MPI_CHAR, win_);		// Attempt to insert value of lockVar into window at displacement rank_
        MPI_Get(waitlist, size_, MPI_CHAR, home_, 0, size_, MPI_CHAR, win_);		// Attempt to get a copy of the waitlist_ array from the home_ rank
        MPI_Win_unlock(home_, win_);																	// Unlock the window

        assert(waitlist[rank_] == MPI_MUTEX_LOCKED);										// Verify that this process has managed to PUT its lock intention into waitlist_

        // For each process that has the lock before this process, wait to RECV an empty note that it has finished
        for (i = 0; i < size_; i++) {
            if (waitlist[i] == MPI_MUTEX_LOCKED && i != rank_) {
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
        unsigned char lockVar = MPI_MUTEX_LOCKED;
        int i;

        // Try to acquire lock in one access epoch
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, home_, 0, win_);
        MPI_Put(&lockVar, 1, MPI_CHAR, home_, rank_ /* &win[rank_] */, 1, MPI_CHAR, win_);
        MPI_Get(waitlist, size_, MPI_CHAR, home_, 0, size_, MPI_CHAR, win_);
        MPI_Win_unlock(home_, win_);

        assert(waitlist[rank_] == MPI_MUTEX_LOCKED);

        // Count the 1's
        for (i = 0; i < size_; i++) {
            if (waitlist[i] == MPI_MUTEX_LOCKED && i != rank_) {
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
        unsigned char lockVar = MPI_MUTEX_UNLOCKED;
        int i, next;

        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, home_, 0, win_);
        MPI_Put(&lockVar, 1, MPI_CHAR, home_, rank_, 1, MPI_CHAR, win_);
        MPI_Get(waitlist, size_, MPI_CHAR, home_, 0, size_, MPI_CHAR, win_);
        MPI_Win_unlock(home_, win_);
        
        assert(waitlist[rank_] == MPI_MUTEX_UNLOCKED);

        // If there are other processes waiting for the lock, transfer ownership
        next = (rank_ + 1 + size_) % size_;
        for (i = 0; i < size_; i++, next = (next + 1) % size_) {
            if (waitlist[next] == MPI_MUTEX_LOCKED) {
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

} // namespace CurvGrid

} // Dune


#endif //DUNE_MPI_MUTEX_HH_
