#ifndef CURVGRID_REAL_TIME_LOG_HH_
#define CURVGRID_REAL_TIME_LOG_HH_

/** Include headers. */
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>

// Threads
#include <thread>
#include <mutex>
#include <future>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/common/mpimutex.hh>


namespace Dune
{

namespace CurvGrid {



#ifdef HAVE_CURVGRID_LINUX_MEMORY_LOG


/** Memory Statistics - Logs memory used by each process individually into a file
 *
 * Algorithm:
 * 1) Communicate PiD of all processes to the master process
 * 2) Start thread on the master process
 * 2.1) Thread on master process reads memory data for all processes every couple of seconds
 * 2.2) It writes info to the file for every 100 sample times or so
 * 3) At the end of code, finish thread
 * 3.1) Call thread.join()
 * 3.2) Set boolean variable on=false
 * 3.3) Thread finishes last memory read, writes remainder to the file, and exits
 *
 *
 * */
class CurvGridRealTimeLog {

	typedef CurvGridRealTimeLog   This;
	typedef std::pair<time_t, std::string> NoteType;

	CurvGridRealTimeLog() {}
	CurvGridRealTimeLog(This const&)  = delete;          // Disallow copy-constructors
    void operator=(This const&)  = delete;

public:

	/** \brief Container for memory consumption data. */
	struct memoryType {
		int size_;          // total program size, in kb
		int res_;           // size of memory portions, in kb
		int shared_;        // memory of pages that are shared, in kb
		int text_;          // memory of pages that are code, in kb
		int sharedLibs_;    // memory of pages of data/stack, in kb
		int stack_;         // memory of library pages, in kb
		int dirtyPages_;    // memory of dirty pages, in kb
	};


public:

    static This & getInstance()
    {
        static This instance;
        return instance;
    }


    // Initialize the memory logging mechanism. This starts a memory logging thread on the master process
    static void init(Dune::MPIHelper & mpihelper, std::string filenameBase, int pauseIntervalSeconds = 1, int packageWriteSize = 5) {
    	getInstance().initImpl(mpihelper, filenameBase, pauseIntervalSeconds, packageWriteSize);
    }

    // Note a certain point in time, to later associate the memory consumption
    // This function can be called multiple times, on any process
    static void note(std::string text) { getInstance().noteImpl(text); }

    // Logs the value of a variable to the dynamic output file
    static void logvar(std::string name, std::string val) {  getInstance().logvarImpl(name, val); }

    // Finish the memory logging. Stops the memory logging thread, and writes all unwritten data to a log file
    static void stop() { getInstance().stopImpl(); }


protected:

    static time_t timeNow() { return std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); }

    // Convert time to string and remove the unnecessary end of line at the end of the timestring
    static std::string time2str (time_t t) {
    	std::string strtime = asctime(localtime (&t));
    	strtime.resize(strtime.size() - 1);
    	return strtime;
    }

	void initImpl(Dune::MPIHelper & mpihelper, std::string filenameBase, int pauseIntervalSeconds, int packageWriteSize)
	{
		rank_ = mpihelper.rank();
		size_ = mpihelper.size();
		pauseIntervalSeconds_ = pauseIntervalSeconds;
		packageWriteSize_ = packageWriteSize;
		filename_ = filenameBase + ".mem";
		mpimutex_ = new MPIMutex(rank_, size_, MPI_MASTER_RANK);

		int pid = getpid();
		std::vector<int> pidGlobal (mpihelper.size());

		MPI_Comm comm = Dune::MPIHelper::getCommunicator();
		MPI_Gather(&pid, 1, MPI_INT, static_cast<int *>(pidGlobal.data()), 1, MPI_INT, MPI_MASTER_RANK, comm);

		if (rank_ == MPI_MASTER_RANK) {
			// Create the memory file
			std::fstream memFile;
			memFile.open(filename_, std::ios::out);
			memFile.close();

			on_ = true;
			fi_ = std::async(std::launch::async, &CurvGridRealTimeLog::memLogLoop, this, pidGlobal, std::ref(on_));
		}
	}


	// Write a note of arbitrary text to the log file. Useful to determine where exactly the code is at this moment of time
	void noteImpl(std::string text) {
		std::stringstream sstr;
		sstr << time2str(timeNow())  << "; " << rank_ << "; " << "NOTE" << "; " << text;
		writeFileSecure(std::vector<std::string>(1, sstr.str()));
	}

	// Save the name and value of a variable at a given point in time
	void logvarImpl(std::string name, std::string val) {
		std::stringstream sstr;
		sstr << time2str(timeNow())  << "; " << rank_ << "; " << "VAR" << "; " << name << "; " << val;
		writeFileSecure(std::vector<std::string>(1, sstr.str()));
	}


	void stopImpl() {
		//std::cout << "[[[Closing Real time log...." << std::endl;

		// Finish the thread and write all remaining memory data to the file
		if (rank_ == MPI_MASTER_RANK) {
			on_ = false;
			fi_.get();
			//t.join();
		}

		// Delete the mutex
		assert(mpimutex_);  // If init was called, the mpimutex should be initialized at this stage
		delete mpimutex_;

		//std::cout << "[[[Closed Real time log...." << std::endl;
	}



	// DO NOT MAKE THIS VECTOR A REFERENCE!!! THE REFERENCE OWNER WILL BE DESTROYED BEFORE THIS FINISHES
	int memLogLoop(std::vector<int> pidGlobal, std::atomic<bool> & on)
	{
		std::cout << "Starting memory log loop for processes: ";
		for (int i = 0; i < pidGlobal.size(); i++) { std::cout << pidGlobal[i] << " "; }
		std::cout << std::endl;

		std::vector<time_t> memTimes;
		std::vector<std::vector<int>> memData;

		unsigned int pageSize = getpagesize() / 1024;

		while (on) {
			std::vector<int> thisMem;
			for (int i = 0; i < pidGlobal.size(); i++)  { thisMem.push_back(getMemoryTotal(pidGlobal[i]) * pageSize); }
			memData.push_back(thisMem);
			memTimes.push_back(timeNow());
			if (memData.size() >= packageWriteSize_) { writeMem(memData, memTimes); }

			std::this_thread::sleep_for(std::chrono::seconds(pauseIntervalSeconds_));
		}

		// Write last stored memory entries before exit
		writeMem(memData, memTimes);
		return 0;
	}


	/** \brief Get all memory of a process with given pid. */
	static int getMemoryTotal(int pid)
	{
	  std::ifstream memfile;
	  memfile.open("/proc/" + std::to_string(pid) + "/statm");

	  memoryType mem;
	  if (memfile.is_open()) {
	      memfile
	        >> mem.size_
	        >> mem.res_
	        >> mem.shared_
	        >> mem.text_
	        >> mem.sharedLibs_
	        >> mem.stack_
	        >> mem.dirtyPages_;

	      memfile.close();

	      return mem.res_;
	  } else {
		  std::cout << "/proc/" + std::to_string(pid) + "/statm" << std::endl;
		  DUNE_THROW(Dune::IOError, "Linux process file not find");
	  }
	}


	// Write memory data to a file and clear the memory data
	void writeMem(std::vector<std::vector<int>> & memData, std::vector<time_t> & memTimes) {
		std::vector<std::string> datastr(memData.size(), "");
		assert(memTimes.size() == memData.size());
		for (int i = 0; i < memData.size(); i++) {
			std::stringstream str;
			str << time2str(memTimes[i])  << "; " << rank_ << "; " << "MEM" << "; ";
			for (int j = 0; j < memData[i].size() - 1; j++) { str << memData[i][j] << "; "; }
			str << memData[i][memData[i].size() - 1];  // write the last one without separator at the end
			datastr[i] = str.str();
		}

		writeFileSecure(datastr);

		memData.clear();
		memTimes.clear();
	}


	// Write parallel data to file, using mutexes to avoid clashes
	// Can write several instances of one command to save on file writing
	void writeFileSecure(const std::vector<std::string> & output) {
		// If this rank has threads writing to the memory file, lock file output for threads
		std::lock_guard<std::mutex> * lock;
		if (rank_ == MPI_MASTER_RANK)  { lock = new std::lock_guard<std::mutex>(threadmutex_); }

		// Lock file output for all ranks except this one, so other ranks can't write to file
		mpimutex_->lock();

		std::fstream outfile;
		outfile.open(filename_, std::ios::out | std::ios::app);
		for (int i = 0; i < output.size(); i++) {
			outfile << output[i] << std::endl;
		}
		outfile.close();

		// Unlock MPI Lock
		mpimutex_->unlock();

		// Remove thread-lock
		if (rank_ == MPI_MASTER_RANK)  { delete lock; }
	}



private:


	// MPI
	int rank_;
	int size_;
	MPIMutex * mpimutex_;

	// Threads
	std::atomic<bool> on_;										// Determines if the memory log loop is on or off
	std::future<int> fi_;											// The asynchronous thread for the memory loop
	std::mutex threadmutex_;									// The mutex for protecting output of threads on the same rank

	// Memory data
	std::string filename_;
	int pauseIntervalSeconds_;
	int packageWriteSize_;

};



#else


// Fake implementation. Can use as skeleton for impl on other operating systems, e.g. MAC
class CurvGridRealTimeLog {
	typedef CurvGridRealTimeLog   This;

	CurvGridRealTimeLog() {}
	CurvGridRealTimeLog(This const&)  = delete;          // Disallow copy-constructors
    void operator=(This const&)  = delete;
public:

    static This & getInstance()
    {
        static This instance;
        return instance;
    }

	static void init(Dune::MPIHelper & mpihelper, std::string filename, int pauseIntervalSeconds = 1, int packageWriteSize = 5) { }

	static void note(std::string text) { }

	static void logvar(std::string name, std::string val) {  }

	static void stop() { }
};



#endif



} // namespace CurvGrid

}  // End of namespace Dune

#endif // CURVGRID_REAL_TIME_LOG_HH_

