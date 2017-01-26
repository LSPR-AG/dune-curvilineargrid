#ifndef DUNE_CURVILINEARGRID_CONSTANT_HH
#define DUNE_CURVILINEARGRID_CONSTANT_HH

#include <stdint.h>

namespace Dune
{

namespace CurvGrid
{


// ******************************************************/
// Arithmetic Constants                                       **/
// ******************************************************/
const double NUMERICS_RELATIVE_TOLERANCE = 1.0e-10;  // Maximal error allowed to be accumulated by double Mantissa through basic arithmetic


// ******************************************************/
// MPI Constants                                       **/
// ******************************************************/
const int MPI_MASTER_RANK = 0;

// ******************************************************/
// Curvilinear Grid Constants                          **/
// ******************************************************/
const int PERIODIC_GHOST_BIND_EDGE_TYPE = 1001;  // Additional partition type for periodic ghost binding edges, used for VTK visualization
const int PERIODIC_GHOST_PARTITION_TYPE = 1000;  // Additional partition type for periodic ghosts
const int PERIODIC_BOUNDARY_PARTITION_TYPE = 600;  // Additional partition type for periodic boundary faces
const int BOUNDARY_SEGMENT_PARTITION_TYPE = 500;  // Additional partition type for boundary segments to distinguish them from interior faces
const int INTERIOR_BOUNDARY_SEGMENT_PARTITION_TYPE = 400;  // Additional partition type for interior boundary surfaces

// ******************************************************/
// Constants associated with the LoggingMessage        **/
// ******************************************************/
const bool LOG_MSG_PATH_ON  = true;  const bool LOG_MSG_PATH_OFF = false;
const bool LOG_MSG_RANK_ON  = true;  const bool LOG_MSG_RANK_OFF = false;

// ******************************************************/
// Curvilinear VTK Writer Constants        **/
// ******************************************************/
const std::string VTK_XML_VERSION = "1.0";
const std::string VTK_GRID_TYPE = "UnstructuredGrid";
const std::string VTK_VTU_VERSION = "0.1";
const std::string VTK_BYTE_ORDER = "LittleEndian";

const std::string VTK_DATA_FORMAT_ASCII = "ASCII";
const std::string VTK_DATA_FORMAT_BINARY = "BINARY";
const std::string VTU_DATA_FORMAT_ASCII = "ascii";
const std::string VTU_DATA_FORMAT_BINARY = "binary";

const std::string VTU_FLOAT_TYPE = "Float32";
const std::string VTU_UINT_8_TYPE = "UInt8";
const std::string VTU_INT_32_TYPE = "Int32";

typedef float  VTKFloat;
typedef uint8_t VTKUInt8;
typedef int32_t VTKInt32;



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
