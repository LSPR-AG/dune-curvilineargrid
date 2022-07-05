# File for module specific CMake tests.

#########################################
# Find necessary packages               #
#########################################


#########################################
# Find ParMETIS                         #
#########################################
#find_package(ParMETIS)


#########################################
# Find Boost                            #  # nota bene: Boost is no longer necessary, since c++ <chrono> has all necessary functionality
#########################################
# find_package(Boost)              # find_package(Boost 1.55.0 COMPONENTS system)
# if(Boost_FOUND)
#   include_directories(${Boost_INCLUDE_DIRS})
#endif()
