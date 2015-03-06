# Define variable with path to the grids to be used in the code
add_definitions("-DDUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"")

# Define tests to be run
set(TESTS
	tutorial1gettingstarted
	tutorial2traversal
	tutorial3visualisation
)

add_executable(tutorial1gettingstarted tutorial1gettingstarted.cc)
add_executable(tutorial2traversal tutorial2traversal.cc)
add_executable(tutorial3visualisation tutorial3visualisation.cc)
target_link_libraries(tutorial1gettingstarted dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial2traversal dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial3visualisation dunecommon dunegeometry dunegrid) # dunecurvilineargeometry


# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})


# Add parmetis
if(HAVE_MPI)
    add_dune_mpi_flags(tutorial1gettingstarted)
    add_dune_mpi_flags(tutorial2traversal)
    add_dune_mpi_flags(tutorial3visualisation)
    add_dune_parmetis_flags(tutorial1gettingstarted)
    add_dune_parmetis_flags(tutorial2traversal)
    add_dune_parmetis_flags(tutorial3visualisation)
endif(HAVE_MPI)

# Add boost
if(Boost_FOUND)
   target_link_libraries(tutorial1gettingstarted ${Boost_LIBRARIES})
   target_link_libraries(tutorial2traversal ${Boost_LIBRARIES})
   target_link_libraries(tutorial3visualisation ${Boost_LIBRARIES})
endif()
