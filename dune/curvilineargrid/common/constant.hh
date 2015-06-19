#ifndef DUNE_CURVILINEARGRID_CONSTANT_HH
#define DUNE_CURVILINEARGRID_CONSTANT_HH

namespace Dune
{

namespace CurvGrid
{


// ******************************************************/
// MPI Constants                                       **/
// ******************************************************/
const int MPI_MASTER_RANK = 0;



// ******************************************************/
// Curvilinear Grid Constants                          **/
// ******************************************************/
const int BOUNDARY_SEGMENT_PARTITION_TYPE = 500;  // Additional partition type for boundary segments to distinguish them from interior faces







// ******************************************************/
// Constants associated with the LoggingMessage        **/
// ******************************************************/
const bool LOG_MSG_PATH_ON  = true;  const bool LOG_MSG_PATH_OFF = false;
const bool LOG_MSG_RANK_ON  = true;  const bool LOG_MSG_RANK_OFF = false;

enum
{
	LOG_MSG_PERSISTENT = 0,     // This output is always written regardless of verbosity options
	LOG_MSG_PRODUCTION = 1,     // This output is written if LoggingMessage is initialized with verbosity

	// Debugging verbosity is activated if LoggingMessage is initialized with verbosity, and if the corresponding flag has been compiled in
	// The intended use is to increase the verbosity level within each loop, thus the amount of output will decrease exponentially with each order
	LOG_MSG_DVERB = 2,          // Debugging Verbosity lvl=0
	LOG_MSG_DDVERB = 3,          // Debugging Verbosity lvl=1
	LOG_MSG_DDDVERB = 4,          // Debugging Verbosity lvl=2
	LOG_MSG_DDDDVERB = 5,          // Debugging Verbosity lvl=3
	LOG_MSG_DDDDDVERB = 6,          // Debugging Verbosity lvl=4
	LOG_MSG_DDDDDDVERB = 7,          // Debugging Verbosity lvl=5
};


} // namespace CurvGrid

} // namespace Dune






#endif
