#ifndef DUNE_LOGGING_TIMER_HH
#define DUNE_LOGGING_TIMER_HH

#include <config.h>

/** Include headers. */
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <ratio>
#include <chrono>
#include <map>


#include <dune/common/parallel/mpihelper.hh>

using namespace std::chrono;

namespace Dune
{

template<typename LoggingMessage>
class LoggingTimer
{
	typedef LoggingTimer<LoggingMessage>  This;

	LoggingTimer()  {  }
	LoggingTimer(This const&)  = delete;
    void operator=(This const&)  = delete;

    typedef high_resolution_clock::time_point  TimePoint;
    typedef duration<double>                   TimeInterval;

    typedef std::map<std::string, unsigned int>   ActionNameMap;
    typedef typename ActionNameMap::iterator      ActionNameIter;
    typedef std::pair<std::string, unsigned int>  ActionNameKey;

    static const unsigned int  DEFAULT_N_CHAR_ITER = 6;   // Default number of characters used for number of iterations of a timed operation
    static const unsigned int  DEFAULT_N_CHAR_TIME = 13;  // Default number of characters used for time output as a string

    const bool DEFAULT_REALVERBOSE = false;  // Do not write text every time something is timed, unless user really wants to see a lot of text output

    /***********************************************************/
    /* Auxiliary structures                                    */
    /***********************************************************/

    // This structure stores active timed actions
    struct TimeIntervalStorage
    {
    	bool         on_;
    	TimePoint    first_;
    	TimePoint    last_;
    	TimeInterval duration_;
    	unsigned int count_;
    };


    struct TimeStatisticsStorate
    {
    	TimePoint    start_;
    	TimeInterval minDuration_;
    	TimeInterval maxDuration_;
    	unsigned int minCount_;
    	unsigned int maxCount_;
    	std::string  actionname_;

    	bool operator< (const TimeStatisticsStorate & other) const { return start_ < other.start_; }
    };

    typedef typename std::pair<std::string, TimeIntervalStorage>  NameIntervalPair;

    static bool comparefirst (const NameIntervalPair & i, const NameIntervalPair & j)  { return (i.first < j.first); }

public:

    /***********************************************************/
    /* Interface                                               */
    /***********************************************************/

    static This & getInstance()
    {
        static This instance;
        return instance;
    }

    /** \brief Static singleton initialization */
    static void init(Dune::MPIHelper & mpihelper) { getInstance().initImpl(mpihelper); }

    static void setRealVerbose(bool realverbose)  { getInstance().setRealVerboseImpl(realverbose); }


    /** \brief Times an operation. This operation has to be called an even number of times for each keyword
     * First call records the starting time of an operation, 2nd call records the elapsed time during the operation
     * The consequent call pairs calculate the number of counts of each keyword and sum up the total time */
    static void time(std::string actionName)  { getInstance().timeImpl(actionName); }


    /** \brief Checks if a given action has ever been timed  **/
    static bool isTimed(std::string actionName) {
    	return getInstance().isTimedImpl(actionName);
    }


    /** \brief Reports all operations timed on this process. This is a serial operation: if called in parallel, each process will report */
    static void report() { getInstance().reportImpl(); }


    /** Communicates all finished timers to master process.
    /* Outputs minimal and maximal duration for each action over all processes on master process
    /* NOTE: this procedure requires that EXACTLY THE SAME action names are timed on all processes
    /* NOTE: The order of output is sorted by the starting time of action on master process
    /* NOTE: Strings are not communicated. Only the total number of timers is compared. Communication is based on sorting strings on each process */
    static void reportParallel() { getInstance().reportParallelImpl(); }


protected:

    /***********************************************************/
    /* Static initialization                                   */
    /***********************************************************/

    void initImpl(Dune::MPIHelper & mpihelper) {
    	mpihelper_ = &mpihelper;
    	realverbose_ = DEFAULT_REALVERBOSE;
    }

    void setRealVerboseImpl(bool realverbose)  { realverbose_ = realverbose; }



    /***********************************************************/
    /* Auxiliary Methods                                       */
    /***********************************************************/

    static std::string removeNewline(const std::string & s)
    {
    	std::string newS(s);
    	newS.erase(std::remove(newS.begin(), newS.end(), '\n'), newS.end());
    	return newS;
    }


    template<typename T>
    static std::string niceObjStr(T obj, int newSize)
    {
    	std::stringstream s;
    	s << obj;
    	std::string str = s.str();

    	if (str.size() < newSize)  { str.insert(str.size(), "                ", newSize - str.size()); }

    	return str;
    }


    /***********************************************************/
    /* Implementation                                          */
    /***********************************************************/

    bool isTimedImpl(std::string actionName) {
    	return namemap_.find(actionName) != namemap_.end();
    }


    // Starts timer by mapping the action name to current time
    // If action already mapped, notes duration of the action
    // NOTE: names of all actions MUST BE unique
    void timeImpl(std::string actionName)
    {
    	TimePoint  timeNow = high_resolution_clock::now();

    	ActionNameIter iter = namemap_.find(actionName);
    	if (iter == namemap_.end())  // If this is the first time this keyword was used
    	{
    		TimeIntervalStorage interv;
    		interv.on_    = true;      // Currently being timed
    		interv.first_ = timeNow;   // Init first time this operation has been recorded
    		interv.last_  = timeNow;   // Init last time this operation has been recorded
    		interv.count_ = 0;         // This operation has not been timed a single time yet
    		namemap_.insert(ActionNameKey(actionName, timestorage_.size()));   // Map this action name to its timing info index
    		timestorage_.push_back(NameIntervalPair(actionName, interv));      // Store this timing info

    		if (realverbose_)  {
    			std::stringstream logstr;
    			std::time_t tt = system_clock::to_time_t(timeNow);
    			logstr << "-----Started Action: (" << actionName << ") at " << ctime(&tt) << std::endl;
    			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, logstr.str());
    		}
    	} else  // If this keyword was used before
    	{
    		TimeIntervalStorage & interv = timestorage_[(*iter).second].second;

    		if (interv.on_) {  // If this operation is currently being timed
    			interv.on_ = false;      // Switch off the timer
    			interv.count_++;         // Increase the number of times this action has been timed

    			TimeInterval  thisInterv = duration_cast<TimeInterval>(timeNow - interv.last_);   // Calculate the time taken during last interval

    			interv.last_     = timeNow;  // Record the last time the timer has been called (for statistics purposes)

    			if (interv.count_ == 1)  { interv.duration_ = thisInterv; }    // Initialize duration
    			else                     { interv.duration_ += thisInterv; }   // Increase duration

        		if (realverbose_)  {
        			std::stringstream logstr;
        			std::time_t tt = system_clock::to_time_t(timeNow);
        			logstr << "-----Finished Action: (" << actionName << ") at " << ctime(&tt) << std::endl;
        			LoggingMessage::template write<CurvGrid::LOG_MSG_DVERB>(__FILE__, __LINE__, logstr.str());
        		}
    		} else             // If this operation was timed before, but is not currently being timed
    		{
    			interv.on_ = true;       // Switch on the timer
    			interv.last_ = timeNow;  // Record the starting time
    		}
    	}
    }


    // Writes currently active and finished timers locally for each processor
    void reportImpl()
    {
    	LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, "---- Starting current timer status ------");

    	LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, "------ Currently active timers: ---------");
    	for (ActionNameIter iter = namemap_.begin(); iter != namemap_.end(); iter++)
    	{
    		TimeIntervalStorage & interv = timestorage_[(*iter).second].second;

    		if (interv.on_)
    		{
        		std::stringstream logstr;
        		logstr << "Action: " << (*iter).first << " has been called " << interv.count_ + 1 << " times and is still active";
        		LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, logstr.str());
    		}
    	}

    	LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, "------ Currently finished timers: ---------");
    	for (ActionNameIter iter = namemap_.begin(); iter != namemap_.end(); iter++)
    	{
    		TimeIntervalStorage & interv = timestorage_[(*iter).second].second;

    		if (!interv.on_)
    		{
        		std::stringstream logstr;

        		//logstr << "Started at " << ctime(&interv.first_) << " ||| ";
        		//logstr << " finished at " << ctime(&interv.last_) << " ||| ";
        		logstr << "Called "         << niceObjStr(interv.count_,            DEFAULT_N_CHAR_ITER) << " times ||| ";
        		logstr << "Total duration " << niceObjStr(interv.duration_.count(), DEFAULT_N_CHAR_TIME) << " ||| ";
        		logstr << " action: " << (*iter).first;

    			LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, removeNewline(logstr.str()));
    		}
    	}
    	LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, "---- Finishing current timer status ------");
    }


    // Communicates all finished timers to master process.
    void reportParallelImpl()
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_->getCollectiveCommunication();
    	int rank = mpihelper_->rank();
    	int size = mpihelper_->size();

    	//for (int i = 0; i < timestorage_.size(); i++)  { std::cout << "process [" << rank << " operation " << timestorage_[i].first << std::endl;  }

    	// -------------------Self-Test-------------------------------------------
    	// Communicate total number of finished processes to master. Assert that this number is the same on all processes
    	// [TODO] Calculate number of started and finished processes separately
    	{
    		int thisNTimer = timestorage_.size();
    		std::vector<int> nTimer(size);
    		collective_comm.gather(&thisNTimer, nTimer.data(), 1, CurvGrid::MPI_MASTER_RANK);

    		if (rank == CurvGrid::MPI_MASTER_RANK) {
    			for (int i = 0; i < size; i++)  { assert(thisNTimer == nTimer[i]); }
    		}
    	}
    	// -------------------Self-Test Over-------------------------------------------

    	// Copy the time-storage, as it will be sorted in a moment
    	std::vector<NameIntervalPair> tmpstorage = timestorage_;

    	// Sort all timelog_ entries by action name
    	std::sort(tmpstorage.begin(), tmpstorage.end(), comparefirst);

    	// Communicate start and finish time to master
    	std::vector<TimeStatisticsStorate> timeIntervalStorage;

    	for (int i = 0; i < tmpstorage.size(); i++)
    	{
    		TimeIntervalStorage & interv = tmpstorage[i].second;

    		std::vector<TimePoint>     startTimeSet;
    		std::vector<TimeInterval>  durationTimeSet;
    		std::vector<unsigned int>  countTimeSet;

    		if (rank == CurvGrid::MPI_MASTER_RANK)  {
    			startTimeSet.resize(size);
    			durationTimeSet.resize(size);
    			countTimeSet.resize(size);
    		}

    		collective_comm.gather(& interv.first_,    startTimeSet.data(),    1, CurvGrid::MPI_MASTER_RANK);
    		collective_comm.gather(& interv.duration_, durationTimeSet.data(), 1, CurvGrid::MPI_MASTER_RANK);
    		collective_comm.gather(& interv.count_,    countTimeSet.data(),    1, CurvGrid::MPI_MASTER_RANK);


    		if (rank == CurvGrid::MPI_MASTER_RANK) {
    			// Initialize statistics
    			TimeStatisticsStorate  thisStatStorage;
    			thisStatStorage.actionname_  = tmpstorage[i].first;
    			thisStatStorage.start_       = interv.first_;
    			thisStatStorage.minDuration_ = interv.duration_;
    			thisStatStorage.maxDuration_ = interv.duration_;
    			thisStatStorage.minCount_    = interv.count_;
    			thisStatStorage.maxCount_    = interv.count_;

    			// Compute min and max over all processes
    			for (int i = 0; i < size; i++)
    			{
        			thisStatStorage.minDuration_ = std::min(thisStatStorage.minDuration_, durationTimeSet[i]);
        			thisStatStorage.maxDuration_ = std::max(thisStatStorage.maxDuration_, durationTimeSet[i]);
        			thisStatStorage.minCount_    = std::min(thisStatStorage.minCount_, countTimeSet[i]);
        			thisStatStorage.maxCount_    = std::max(thisStatStorage.maxCount_, countTimeSet[i]);
    			}
    			timeIntervalStorage.push_back(thisStatStorage);
    		}
    	}

		// Sort all received timelog according to starting time on master
    	std::sort(timeIntervalStorage.begin(), timeIntervalStorage.end());

    	// Write all output
		if (rank == CurvGrid::MPI_MASTER_RANK) {
			LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, "--------- Parallel Timer Statistics ------------------");

			for (int i = 0; i < timeIntervalStorage.size(); i++)
    		{
				std::stringstream logstr;
				logstr << " min_time: "  << niceObjStr(timeIntervalStorage[i].minDuration_.count(), DEFAULT_N_CHAR_TIME);
				logstr << " max_time: "  << niceObjStr(timeIntervalStorage[i].maxDuration_.count(), DEFAULT_N_CHAR_TIME);
				logstr << " min_count: " << niceObjStr(timeIntervalStorage[i].minCount_, DEFAULT_N_CHAR_ITER);
				logstr << " max_count: " << niceObjStr(timeIntervalStorage[i].maxCount_, DEFAULT_N_CHAR_ITER);
				logstr << " Action: " << timeIntervalStorage[i].actionname_;
				LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, logstr.str());
    		}
			LoggingMessage::template write<CurvGrid::LOG_MSG_PRODUCTION>(__FILE__, __LINE__, "--------- Parallel Timer Statistics - Finished -------");
		}
    }


private:

    bool realverbose_;
    ActionNameMap namemap_;
    std::vector<NameIntervalPair> timestorage_;

    Dune::MPIHelper * mpihelper_;

};

} // End namespace Dune.

#endif // DUNE_LOGGING_TIMER_HH
