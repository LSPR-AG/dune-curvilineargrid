# Define variable with path to the grids to be used in the code
add_definitions("-DDUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"")
#add_definitions("-DHAVE_LOG_MSG_DVERB")

# Define tests to be run
set(TESTS
  test-diagnostics
  test-diagnostics-base
  test-allcommunicate
  test-vectorhelper
  test-curvgrid
)

add_executable(test-diagnostics test-diagnostics.cc)
add_executable(test-diagnostics-base test-diagnostics-base.cc)
add_executable(test-allcommunicate test-allcommunicate.cc)
add_executable(test-vectorhelper test-vectorhelper.cc)
add_executable(test-curvgrid test-curvgrid.cc)

target_link_libraries(test-diagnostics dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(test-diagnostics-base dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(test-allcommunicate dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(test-vectorhelper dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(test-curvgrid dunecommon dunegeometry dunegrid) # dunecurvilineargeometry


# We do not want want to build the tests during make all,
# but just build them on demand
#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS})

#foreach(_test ${TESTS})
#  add_test(${_test} ${_test})
#endforeach(_test ${TESTS})


# Add parmetis
if(HAVE_MPI)
    add_dune_mpi_flags(test-diagnostics)
    add_dune_mpi_flags(test-diagnostics-base)
    add_dune_mpi_flags(test-allcommunicate)
    add_dune_parmetis_flags(test-diagnostics)
    add_dune_parmetis_flags(test-diagnostics-base)
    add_dune_parmetis_flags(test-curvgrid)
endif(HAVE_MPI)

# Add boost
if(Boost_FOUND)
   target_link_libraries(test-diagnostics ${Boost_LIBRARIES})
   target_link_libraries(test-diagnostics-base ${Boost_LIBRARIES})
   target_link_libraries(test-curvgrid ${Boost_LIBRARIES})
endif()


