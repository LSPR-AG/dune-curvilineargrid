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

template<typename LoggingMessage>
class LoggingTimer
{
	typedef LoggingTimer<LoggingMessage>  This;

	LoggingTimer()  {  }
	LoggingTimer(This const&)  = delete;
    void operator=(This const&)  = delete;



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
    static const unsigned int LOG_CATEGORY_DEBUG = LoggingMessage::Category::DEBUG;


    static This & getInstance()
    {
        static This instance;
        return instance;
    }

    void init(bool realverbose)  { realverbose_ = realverbose; }


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
    			LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    		}

    		timedmap_.insert(TimedMapKey(actionName, tmpPoint));
    	} else
    	{
    		if (realverbose_)  {
    			std::stringstream logstr;
    			std::time_t tt = system_clock::to_time_t(tmpPoint);
    			logstr << "-----Finished Action: (" << actionName << ") at " << ctime(&tt) << std::endl;
    			LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    		}

    		timelog_.push_back(TimedData(actionName, TimedPair((*iter).second, tmpPoint)));
    		timedmap_.erase(iter);
    	}
    }


    // Writes currently active and finished timers locally for each processor
    void report()
    {
    	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "---- Starting current timer status ------");
    	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "------ Currently active timers: ---------");

    	for (TimedMapIter iter = timedmap_.begin(); iter != timedmap_.end(); iter++)
    	{
    		std::time_t tt = system_clock::to_time_t((*iter).second);

    		std::stringstream logstr;
    		logstr << "Action: " << (*iter).first << " time " << ctime(&tt);
    		LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    	}

    	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "------ Currently finished timers: ---------");
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

			LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, removeNewline(logstr.str()));
    	}
    	LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "---- Finishing current timer status ------");
    }


    // Communicates all finished timers to master process.
    // Outputs minimal and maximal duration for each action over all processes on master process
    // NOTE: this procedure requires that EXACTLY THE SAME action names are timed on all processes
    // NOTE: The order of output is sorted by the starting time of action on master process
    // NOTE: Strings are not communicated. Only the total number of timers is compared. Communication is based on sorting strings on each process
    void reportParallel(MPIHelper & mpihelper)
    {
    	Dune::CollectiveCommunication<MPI_Comm> collective_comm = mpihelper.getCollectiveCommunication();
    	int rank = mpihelper.rank();
    	int size = mpihelper.size();

    	// -------------------Self-Test-------------------------------------------
    	// Communicate total number of finished processes to master. Assert that this number is the same on all processes
    	{
        	int thisNTimer = timelog_.size();
        	std::vector<int> nTimer(size);
        	collective_comm.gather(&thisNTimer, nTimer.data(), 1, MPI_MASTER_RANK);

        	if (rank == MPI_MASTER_RANK) {
        		for (int i = 0; i < size; i++)  { assert(thisNTimer == nTimer[i]); }
        	}
    	}
    	// -------------------Self-Test Over-------------------------------------------

    	// Sort all timelog_ entries by action name
    	std::sort(timelog_.begin(), timelog_.end(), comparefirst);

    	// Communicate start and finish time to master
    	std::vector<TimeIntervalStorage> timeIntervalStorage;

    	std::cout << "TmpLogTime: " << timelog_.size() << std::endl;

    	for (int i = 0; i < timelog_.size(); i++)
    	{
    		std::vector<TimePoint> startVec, endVec;
    		if (rank == MPI_MASTER_RANK)  {
    			startVec.resize(size);
    			endVec.resize(size);
    		}

    		collective_comm.gather(& timelog_[i].second.first, startVec.data(), 1, MPI_MASTER_RANK);
    		collective_comm.gather(& timelog_[i].second.second, endVec.data(), 1, MPI_MASTER_RANK);

    		// Compute min and max over all processes
    		if (rank == MPI_MASTER_RANK) {

        		TimeIntervalStorage  thisIntervStorage;
        		TimeInterval tmpInterv = duration_cast<duration<double>>(timelog_[i].second.second - timelog_[i].second.first);

        		thisIntervStorage.actionname_ = timelog_[i].first;
        		thisIntervStorage.start_ = timelog_[i].second.first;
        		thisIntervStorage.min_ = tmpInterv;
        		thisIntervStorage.max_ = tmpInterv;

        		for (int i = 0; i < size; i++)
        		{
        			TimeInterval thisInterv = duration_cast<duration<double>>(endVec[i] - startVec[i]);
        			thisIntervStorage.min_ = std::min(thisIntervStorage.min_, thisInterv);
        			thisIntervStorage.max_ = std::max(thisIntervStorage.max_, thisInterv);
        		}
        		timeIntervalStorage.push_back(thisIntervStorage);
    		}
    	}

		// Sort all received timelog according to starting time on master
    	std::sort(timeIntervalStorage.begin(), timeIntervalStorage.end());

    	// Write all output
		if (rank == MPI_MASTER_RANK) {
			LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "--------- Parallel Timer Statistics ------------------");

			for (int i = 0; i < timeIntervalStorage.size(); i++)
    		{
				std::stringstream logstr;
				logstr << " min_time: " << niceObjStr(timeIntervalStorage[i].min_.count(), 13);
				logstr << " max_time: " << niceObjStr(timeIntervalStorage[i].max_.count(), 13);
				logstr << " Action: " << timeIntervalStorage[i].actionname_;
				LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, logstr.str());
    		}
			LoggingMessage::getInstance().template write<LOG_CATEGORY_DEBUG>(__FILE__, __LINE__, "--------- Parallel Timer Statistics - Finished -------");
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

    	if (str.size() < newSize)  { str.insert(str.size(), "                ", newSize - str.size()); }

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

    bool realverbose_;
    TimedMap timedmap_;
    std::vector<TimedData> timelog_;

};

} // End namespace Dune.

#endif // DUNE_LOGGING_TIMER_HH
