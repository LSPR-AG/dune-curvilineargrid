#ifndef DUNE_LOGGING_MESSAGE_X_HH
#define DUNE_LOGGING_MESSAGE_X_HH

#include <config.h>

/** Include headers. */
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <time.h>


#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvilineargrid/common/constant.hh>
#include <dune/curvilineargrid/common/colorcoding.hh>

//#include <boost/date_time/posix_time/posix_time.hpp>

namespace Dune
{

namespace CurvGrid {


// A specialized set of templated functions that determine whether to write output based on compile parameters
namespace LoggingMessageCompile
{
	template <unsigned int messageCat>
    bool writeMsg(bool writeThisProcess)  { return false; }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_PERSISTENT>(bool writeThisProcess)  { return true; }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_PRODUCTION>(bool writeThisProcess)  { return writeThisProcess; }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_DVERB>(bool writeThisProcess)  {
#if HAVE_LOG_MSG_DVERB || HAVE_LOG_MSG_DDVERB || HAVE_LOG_MSG_DDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDDDVERB
    	return writeThisProcess;
#endif
    	return false;
    }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_DDVERB>(bool writeThisProcess)  {
#if HAVE_LOG_MSG_DDVERB || HAVE_LOG_MSG_DDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDDDVERB
    	return writeThisProcess;
#endif
    	return false;
    }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_DDDVERB>(bool writeThisProcess)  {
#if HAVE_LOG_MSG_DDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDDDVERB
    	return writeThisProcess;
#endif
    	return false;
    }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_DDDDVERB>(bool writeThisProcess)  {
#if HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDDDVERB
    	return writeThisProcess;
#endif
    	return false;
    }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_DDDDDVERB>(bool writeThisProcess)  {
#if HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDDDVERB
    	return writeThisProcess;
#endif
    	return false;
    }

    template <>
    bool writeMsg<CurvGrid::LOG_MSG_DDDDDDVERB>(bool writeThisProcess)  {
#if HAVE_LOG_MSG_DDDDDDVERB
    	return writeThisProcess;
#endif
    	return false;
    }

}










class LoggingMessage
{
	const bool DEFAULT_VERBOSE  = true;
	const bool DEFAULT_PVERBOSE = false;

private:
	typedef LoggingMessage   This;

	LoggingMessage()  {  }
	LoggingMessage(This const&)  = delete;          // Disallow copy-constructors
    void operator=(This const&)  = delete;


public:

    static This & getInstance()
    {
        static This instance;
        return instance;
    }


    /** \brief Initialises the logging message. This is compulsory before calling write routines
     * \param[in]  mpihelper       Parallel MPIHelper class (or its fake) provided by Dune
     *  */
    static void init(MPIHelper & mpihelper)  { getInstance().initImpl(mpihelper); }

    /** \brief Allows user to dynamically change verbosity and parallel verbosity of the code */
    static void setVerbose(bool verbose)    { getInstance().setVerboseImpl(verbose); }
    static void setPVerbose(bool pverbose)  { getInstance().setPVerboseImpl(pverbose); }


    /** \brief Logging message writer. This is a singleton routine intended to be run like LoggingMessage::write(...);
     *
     * \param[in]  filename        Name of the logged file (via C macro command __FILE__)
     * \param[in]  linenumber      Number of the line the writer is called from (via C macro command __LINE__)
     * \param[in]  message         User defined text output message
     *
     *  */
    template <unsigned int messageCat>
    static void write(std::string filename, unsigned int linenumber, std::string message)
    {
    	getInstance().template writeImpl<messageCat>(filename, linenumber, message);
    }

    template <class T>
    static void assertWrite(std::string filename, unsigned int linenumber, T a, T b, std::string message) {
    	if (a != b) {
    		std::stringstream messageEff;
    		messageEff << message << "; expected=" << a << ", received=" << b;
    		getInstance().writeImpl(filename, linenumber, messageEff.str(), true, true);
    		DUNE_THROW(Dune::IOError, message);
    	}
    }


    /** \brief Logging message patience writer. This is a singleton routine intended to be run like LoggingMessage::write(...);
     *
     * Only one process should writePatience at the same time. Normally it is the Master process.
     * If there is need for another process (e.g. rank 3) to write patience data, then use allowParallel=true
     * NOTE: if allowParallel=true, user must ensure that only one process is writing at the same time
     * */
    static void writePatience(std::string message, int iter, int total, bool allowParallel = false ) {
#if HAVE_LOG_MSG_DVERB || HAVE_LOG_MSG_DDVERB || HAVE_LOG_MSG_DDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDVERB || HAVE_LOG_MSG_DDDDDDVERB
    // Do nothing. For debugging purposes patience writer is harmful, as inside log file it overwrites debug output.
    // It is only useful for production/optimization time, where clean report of performance of a fully working code is required
#else
    getInstance().writePatienceImpl(message, iter, total, allowParallel);
#endif
    }


protected:

    //*******************************************************/
    // Initialization routines                             **/
    //*******************************************************/

    void initImpl(MPIHelper & mpihelper)
    {
    	verbose_        = DEFAULT_VERBOSE;   // To log or not to log
    	processVerbose_ = DEFAULT_PVERBOSE;  // If parallel=false, only master process logs, otherwise every process logs

    	rank_ = mpihelper.rank();
    	size_ = mpihelper.size();
    }

    void setVerboseImpl(bool verbose)    { verbose_ = verbose; }
    void setPVerboseImpl(bool pverbose)  { processVerbose_ = pverbose; }

    //*******************************************************/
    // Auxiliary routines                                  **/
    //*******************************************************/

    // Checks if we intend to write anything on this process.
    bool writeThisProcess() const  { return processVerbose_ || (verbose_ && (rank_ == CurvGrid::MPI_MASTER_RANK)); }


    // Creates a string with current time using C-style <time.h>
    std::string time2string() const
    {
    	time_t rawtime;
    	time ( &rawtime );

    	std::string strtime = asctime(localtime ( &rawtime ));
    	strtime.resize(strtime.size() - 1);                      // Remove the unnecessary end of line at the end of the timestring

    	return strtime;
    }


    // Increases the length of a string with white space until it is of desired length
    std::string extendString(std::string s, unsigned int n)  {
    	s.resize(n, ' ');
    	return s;
    }


    //*******************************************************/
    // Write process implementation                        **/
    //*******************************************************/

    void writePatienceImpl(std::string message, int iter, int total, bool allowParallel) {
    	if ((rank_ == CurvGrid::MPI_MASTER_RANK) || allowParallel)
    	{
    		int totalInclusive = (total > 1) ? total - 1 : 1;
    		int percent    = std::round((100.0 * iter     )  / totalInclusive);
    		int percentOld = std::round((100.0 * (iter - 1)) / totalInclusive);

    		if ((iter == 0) || (percent != percentOld))
    		{
    			std::cout << "\r[" << ColorCoding::FG_BLUE << extendString(std::to_string(percent) + "%", 4) << ColorCoding::FG_DEFAULT << "] " << message;
    			if (percent == 100)  { std::cout << std::endl; }
    			else                 { std::cout.flush(); }
    		}
    	}
    }


    template <unsigned int messageCat>
    void writeImpl(std::string filename, unsigned int linenumber, std::string message)
    {
    	if (LoggingMessageCompile::writeMsg<messageCat>(writeThisProcess()))
    	{
    		bool withPath = (messageCat > 1);
    		bool withRank = (messageCat > 1);
    		writeImpl(filename, linenumber, message, withPath, withRank);
    	}
    }


    void writeImpl(const std::string & filename, unsigned int linenumber, const std::string & message, bool withPath, bool withRank) const
    {
        /** Set the stream to create for the message. */
        std::stringstream printedMessage;

        printedMessage << time2string();
        if (withRank)  { printedMessage << "| process[" << rank_ << "] | "; }
        if (withPath)  { printedMessage << std::setw(20) << filename << ":" << std::setw(6) << linenumber; }
        printedMessage << " ::: " << message << std::endl;

        std::cout << printedMessage.str();
    };


private:
    int rank_;
    int size_;

    bool verbose_;
    bool processVerbose_;
};


} // namespace CurvGrid

} // End namespace Dune.

#endif // DUNE_LOGGING_MESSAGE_X_HH
