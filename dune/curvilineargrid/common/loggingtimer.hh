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


namespace Dune
{

  using namespace std::chrono;

class LoggingTimer
{
public:
    static const int MPI_MASTER_RANK = 0;

    typedef high_resolution_clock::time_point  TimePoint;
    typedef duration<double>                   TimeInterval;

    typedef std::map<std::string, TimePoint>   TimedMap;
    typedef std::pair<std::string, TimePoint>  TimedMapKey;
    typedef std::pair<TimePoint, TimePoint>    TimedPair;

    typedef typename TimedMap::iterator        TimedMapIter;
    typedef typename std::pair<std::string, TimedPair>  TimedData;


    // Logging Message Routines
    static const unsigned int LOG_PHASE_DEV = Dune::LoggingMessage::Phase::DEVELOPMENT_PHASE;
    static const unsigned int LOG_CATEGORY_DEBUG = Dune::LoggingMessage::Category::DEBUG;


    LoggingTimer(bool realverbose, MPIHelper & mpihelper, LoggingMessage & loggingmessage )
      : realverbose_(realverbose),
        mpihelper_(mpihelper),
        loggingmessage_(loggingmessage)
    {
    	rank_ = mpihelper.rank();
    	size_ = mpihelper.size();
    }


    // Starts timer by mapping the action name to current time
    // If action already mapped, notes duration of the action
    // NOTE: names of all actions MUST BE unique
    void time(std::string actionName)
    {
    	TimePoint  tmpPoint = high_resolution_clock::now();

    	TimedMapIter iter = timedmap_.find(actionName);
    	if (iter == timedmap_.end())
    	{
    		if (realverbose_)  {
    			std::stringstream logstr;
    			std::time_t tt = system_clock::to_time_t(tmpPoint);
    			logstr << "-----Started Action: (" << actionName << ") at " << ctime(&tt) << std::endl;
    			loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    		}

    		timedmap_.insert(TimedMapKey(actionName, tmpPoint));
    	} else
    	{
    		if (realverbose_)  {
    			std::stringstream logstr;
    			std::time_t tt = system_clock::to_time_t(tmpPoint);
    			logstr << "-----Finished Action: (" << actionName << ") at " << ctime(&tt) << std::endl;
    			loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    		}

    		timelog_.push_back(TimedData(actionName, TimedPair((*iter).second, tmpPoint)));
    		timedmap_.erase(iter);
    	}
    }


    // Writes currently active and finished timers locally for each processor
    void report()
    {
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "---- Starting current timer status ------");
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "------ Currently active timers: ---------");

    	for (TimedMapIter iter = timedmap_.begin(); iter != timedmap_.end(); iter++)
    	{
    		std::time_t tt = system_clock::to_time_t((*iter).second);

    		std::stringstream logstr;
    		logstr << "Action: " << (*iter).first << " time " << ctime(&tt);
    		loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    	}

    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "------ Currently finished timers: ---------");
    	for (int i = 0; i < timelog_.size(); i++)
    	{
    		std::time_t      tt_start = system_clock::to_time_t(timelog_[i].second.first);
    		std::time_t      tt_end = system_clock::to_time_t(timelog_[i].second.second);
    		duration<double> thisInterv = duration_cast<duration<double>>(timelog_[i].second.second - timelog_[i].second.first);

    		std::stringstream logstr;

    		logstr << "Started at " << ctime(&tt_start) << " ||| ";
    		logstr << " finished at " << ctime(&tt_end) << " ||| ";
    		logstr << " duration is " << niceObjStr(thisInterv.count(), 13) << " ||| ";
    		logstr << " action: " << timelog_[i].first << " ||| ";

			loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, removeNewline(logstr.str()));
    	}
    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "---- Finishing current timer status ------");
    }


    // Communicates all finished timers to master process.
    // Outputs minimal and maximal duration for each action over all processes on master process
    // NOTE: this procedure requires that EXACTLY THE SAME action names are timed on all processes
    // NOTE: The order of output is sorted by the starting time of action on master process
    // NOTE: Strings are not communicated. Only the total number of timers is compared. Communication is based on sorting strings on each process
    void reportParallel()
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper_.getCollectiveCommunication();

    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__,  "Started Self-test" );

    	// -------------------Self-Test-------------------------------------------
    	// Communicate total number of finished processes to master. Assert that this number is the same on all processes
    	{
        	int thisNTimer = timelog_.size();
        	std::vector<int> nTimer(size_);
        	collective_comm.gather(&thisNTimer, nTimer.data(), 1, MPI_MASTER_RANK);

        	if (rank_ == MPI_MASTER_RANK) {
        		for (int i = 0; i < size_; i++)  { assert(thisNTimer == nTimer[i]); }
        	}
    	}
    	// -------------------Self-Test Over-------------------------------------------


    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "Started Sort" );

    	// Sort all timelog_ entries by action name
    	std::sort(timelog_.begin(), timelog_.end(), comparefirst);

    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__,  "Started MinMax gather" );

    	// Communicate start and finish time to master
    	std::vector<TimeIntervalStorage> timeIntervalStorage;
    	for (int i = 0; i < timelog_.size(); i++)
    	{
    		std::vector<TimePoint> startVec, endVec;
    		collective_comm.gather(& timelog_[i].second.first, startVec.data(), 1, MPI_MASTER_RANK);
    		collective_comm.gather(& timelog_[i].second.second, endVec.data(), 1, MPI_MASTER_RANK);

    		std::cout << "--finished gather" << std::endl;

    		// Compute min and max over all processes
    		if (rank_ == MPI_MASTER_RANK) {

    			std::cout << "aaa" << std::endl;

        		TimeIntervalStorage  thisIntervStorage;
        		TimeInterval tmpInterv = duration_cast<duration<double>>(timelog_[i].second.second - timelog_[i].second.first);

        		std::cout << "bbbb" << std::endl;

        		thisIntervStorage.actionname_ = timelog_[i].first;
        		thisIntervStorage.start_ = timelog_[i].second.first;
        		thisIntervStorage.min_ = tmpInterv;
        		thisIntervStorage.max_ = tmpInterv;

        		std::cout << "cccc" << std::endl;

        		for (int i = 0; i < size_; i++)
        		{
        			TimeInterval thisInterv = duration_cast<duration<double>>(endVec[i] - startVec[i]);

        			std::cout << "start compare" << std::endl;

        			thisIntervStorage.min_ = std::min(thisIntervStorage.min_, thisInterv);
        			thisIntervStorage.max_ = std::max(thisIntervStorage.max_, thisInterv);

        			std::cout << "finishe compare" << std::endl;
        		}
        		timeIntervalStorage.push_back(thisIntervStorage);
    		}
    	}


    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__,  "Started 2nd sort" );

		// Sort all received timelog according to starting time on master
    	std::sort(timeIntervalStorage.begin(), timeIntervalStorage.end());

    	loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__,  "Started output" );

    	// Write all output
		if (rank_ == MPI_MASTER_RANK) {
			loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "--------- Parallel Timer Statistics ------------------");

			for (int i = 0; i < size_; i++)
    		{
				std::stringstream logstr;
				logstr << " Action: " << timeIntervalStorage[i].actionname_;
				logstr << " min_time: " << timeIntervalStorage[i].min_.count();
				logstr << " max_time: " << timeIntervalStorage[i].max_.count();
				loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    		}
			loggingmessage_.write<LOG_PHASE_DEV, LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "--------- Parallel Timer Statistics - Finished -------");
		}
    }


protected:

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

    	if (str.size() < newSize)  { str.insert(str.size(), " ", newSize - str.size()); }

    	return str;
    }

    static bool comparefirst (TimedData i,TimedData j)  { return (i.first < j.first); }

    struct TimeIntervalStorage
    {
    	TimePoint    start_;
    	TimeInterval min_;
    	TimeInterval max_;
    	std::string  actionname_;

    	bool operator< (const TimeIntervalStorage & other) const { return start_ < other.start_; }
    };


private:

    MPIHelper & mpihelper_;
    LoggingMessage & loggingmessage_;

    int rank_;
    int size_;

    bool realverbose_;
    TimedMap timedmap_;
    std::vector<TimedData> timelog_;

};

} // End namespace Dune.

#endif // DUNE_LOGGING_TIMER_HH
