#ifndef DUNE_LOGGING_MESSAGE_HH
#define DUNE_LOGGING_MESSAGE_HH

#include <config.h>

/** Include headers. */
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


#include <dune/common/parallel/mpihelper.hh>

//#include <boost/date_time/posix_time/posix_time.hpp>

namespace Dune
{








struct LoggingMessageHelper
{
    /** \brief Determine if the logging message is only used in the development
     *  phase or also during the production runs */
    struct Phase
    {
        enum
        {
          DEVELOPMENT_PHASE = 0,
          PRODUCTION_PHASE  = 1
        };
    };
};


/** \brief Provide a logging message.
 * The output behaviour can be controlled with some compilation flags:
 *  HAVE_DEBUG:          We get all output including line number and file
 *                       name; only the root cpu says something.
 *  HAVE_DEBUG_PARALLEL: Instead to the option above you get output from all
 *                       cpus!
 *  If neither HAVE_DEBUG nor HAVE_DEBUG_PARALLEL is defined we get only:
 *  -poduction messages (phase == PRODUCTION_PHASE) without any additional
 *   information like line number and file (except the time, if selected).
 */
/*-----------------------------------------------------------------------
| loggingmessage                                                        |
-----------------------------------------------------------------------*/

template<unsigned int phase>
class LoggingMessage
{
private:
	typedef LoggingMessage<phase> This;
	typedef LoggingMessageHelper  Helper;

	LoggingMessage()  {  }
	LoggingMessage(This const&)  = delete;
    void operator=(This const&)  = delete;


public:
    static const int MPI_MASTER_RANK = 0;

    /** \brief Define the type of logging message */
    struct Category
    {
        enum
        {
          DEBUG    = 0,
          MESSAGE  = 1,
          WARNING  = 2,
          ERROR    = 3,
          TIME     = 4,
          MEMORY   = 5,
          PARALLEL = 6,
          CLEAN    = 7
        };
    };

    static const std::vector<std::string> PHASE_TEXT;
    static const std::vector<std::string> CATEGORY_TEXT;

    static This & getInstance()
    {
        static This instance;
        return instance;
    }


    // Initialises the logging message
    void init(MPIHelper & mpihelper, bool verbose, bool processVerbose)
    {
    	verbose_ = verbose;
    	processVerbose_ = processVerbose;

    	rank_ = mpihelper.rank();
    	size_ = mpihelper.size();
    }


    template <unsigned int messageCat>
    static void writeStatic(std::string filename, unsigned int linenumber, std::string message)
    {
    	getInstance().template write<messageCat>(filename, linenumber, message);
    }


    /** \brief Logging message writer
     *
     * \param[in]  mpihelper       Parallel MPIHelper class (or its fake) provided by Dune
     * \param[in]  verbose         To log or not to log
     * \param[in]  parallel        If parallel=false, only master process logs, otherwise every process logs
     * \param[in]  filename        Name of the logged file (via C macro command __FILE__)
     * \param[in]  linenumber      Number of the line the writer is called from (via C macro command __LINE__)
     * \param[in]  message         User defined text output message
     *
     * [FIXME] - Include boost in CMAKE
     * [FIXME] - How exactly is the memory calculated
     *
     *  */
    template <unsigned int messageCat>
    void write(std::string filename, unsigned int linenumber, std::string message)
    {
        if (processVerbose_ || (verbose_ && (rank_ == MPI_MASTER_RANK)))
        {

            /** Set the stream to create for the message. */
            std::stringstream printedMessage;
            printedMessage.clear();

            /** Set the phase string text: */
            std::string phasestring = PHASE_TEXT[phase];
            std::string categorystring = std::string("");


            /** Set the category string text: */
            if ((messageCat != Category::TIME) && (messageCat != Category::MEMORY))
            {
                categorystring = CATEGORY_TEXT[messageCat];
            } else
            {
                // Get the time ...
                //boost::posix_time::ptime time = boost::posix_time::microsec_clock::local_time();
                // ... and print it out.
                //printedMessage << time << " ::: ";
            }

            if (messageCat != Category::CLEAN)
            {
                printedMessage << std::setw(20) << filename << ":" << std::setw(6) << linenumber;
                printedMessage << " ::: process[" << rank_ << "]";

                printedMessage << " ::: " << phasestring << " ";
                printedMessage << categorystring;
            }

            if (messageCat == Category::MEMORY)
            {
                unsigned int memTotal  = 0; //getMemoryTotal();
                unsigned int memShared = 0; //getMemoryShared();
                printedMessage << " ::: total memory = " << std::setw(10) << memTotal << " [kb]"
                               << " ::: shared memory = " << std::setw(10) << memShared << " [kb]";
                printedMessage << " ::: ";
            }

            printedMessage << message;


            // Only print, if there is something written.
            if (printedMessage.str().length() > 0 )  { std::cout << printedMessage.str() << std::endl; }
            else                                     { std::cout << std::endl; }
        }

    };


private:
    int rank_;
    int size_;

    bool verbose_;
    bool processVerbose_;
};





/** Define the text output for LoggingMessagePhase. */
template<unsigned int phase>
const std::vector<std::string> LoggingMessage<phase>::PHASE_TEXT
(
	{
		"DEVELOPMENT",
		"PRODUCTION "
	}
);


/** Define the text output for LoggingMessageCategory. */
template<unsigned int phase>
const std::vector<std::string> LoggingMessage<phase>::CATEGORY_TEXT
(
	{
    	"DEBUG   ",
    	"MESSAGE ",
    	"WARNING ",
    	"ERROR   ",
    	"TIME    ",
    	"MEMORY  ",
    	"PARALLEL",
    	"CLEAN   "
	}
);

} // End namespace Dune.

#endif // DUNE_LOGGING_MESSAGE_HH
