# File for module specific CMake tests.
find_package(ParMETIS)
# find_package(Boost 1.55.0 COMPONENTS system)
find_package(Boost)

include(AddParMETISFlags)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
endif()
