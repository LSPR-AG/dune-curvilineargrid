# Define variable with path to the grids to be used in the code
add_definitions("-DDUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"")

# Logging for debug purposes
# add_definitions("-DHAVE_LOG_MSG_DVERB")

# Define tests to be run
set(TESTS
  curvilineargmshtest)

add_executable(curvilineargmshtest curvilineargmshtest.cc)
target_link_libraries(curvilineargmshtest PUBLIC dunecommon dunegeometry dunegrid)  # dunecurvilineargeometry


# We do not want want to build the tests during make all,
# but just build them on demand
#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS})

#foreach(_test ${TESTS})
#  add_test(${_test} ${_test})
#endforeach(_test ${TESTS})


# Add parmetis
if(HAVE_MPI)
    add_dune_parmetis_flags(curvilineargmshtest)
endif()

# Add boost
if(Boost_FOUND)
   target_link_libraries(curvilineargmshtest ${Boost_LIBRARIES})
endif()


