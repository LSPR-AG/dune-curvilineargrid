# Define variable with path to the grids to be used in the code
add_definitions("-DDUNE_CURVILINEARGRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"")
# add_definitions("-DHAVE_LOG_MSG_DVERB")

# Define tests to be run
set(TESTS
	tutorial1gettingstarted
	tutorial2traversal
	tutorial3visualisation
	tutorial3visualisation
	tutorial4IntegrationGauss
	tutorial4IntegrationNormal
	tutorial5CommunicationConstant
	tutorial5CommunicationGlobalIndex
	tutorial6ParallelData
	tutorial7GlobalBoundary
	tutorial8InteriorBoundary
	tutorial9PeriodicBoundary
)

add_executable(tutorial1gettingstarted tutorial1gettingstarted.cc)
add_executable(tutorial2traversal tutorial2traversal.cc)
add_executable(tutorial3visualisation tutorial3visualisation.cc)
add_executable(tutorial4IntegrationGauss tutorial4IntegrationGauss.cc)
add_executable(tutorial4IntegrationNormal tutorial4IntegrationNormal.cc)
add_executable(tutorial5CommunicationConstant tutorial5CommunicationConstant.cc)
add_executable(tutorial5CommunicationGlobalIndex tutorial5CommunicationGlobalIndex.cc)
add_executable(tutorial6ParallelData tutorial6ParallelData.cc)
add_executable(tutorial7GlobalBoundary tutorial7GlobalBoundary.cc)
add_executable(tutorial8InteriorBoundary tutorial8InteriorBoundary.cc)
add_executable(tutorial9PeriodicBoundary tutorial9PeriodicBoundary.cc)
target_link_libraries(tutorial1gettingstarted dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial2traversal dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial3visualisation dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial4IntegrationGauss dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial4IntegrationNormal dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial5CommunicationConstant dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial5CommunicationGlobalIndex dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial6ParallelData dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial7GlobalBoundary dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial8InteriorBoundary dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial9PeriodicBoundary dunecommon dunegeometry dunegrid) # dunecurvilineargeometry


# We do not want want to build the tests during make all,
# but just build them on demand
#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS})

#foreach(_test ${TESTS})
#  add_test(${_test} ${_test})
#endforeach(_test ${TESTS})


# Add parmetis
if(HAVE_MPI)
    add_dune_mpi_flags(tutorial1gettingstarted)
    add_dune_mpi_flags(tutorial2traversal)
    add_dune_mpi_flags(tutorial3visualisation)
    add_dune_mpi_flags(tutorial4IntegrationGauss)
    add_dune_mpi_flags(tutorial4IntegrationNormal)
    add_dune_mpi_flags(tutorial5CommunicationConstant)
    add_dune_mpi_flags(tutorial5CommunicationGlobalIndex)
    add_dune_mpi_flags(tutorial6ParallelData)
    add_dune_mpi_flags(tutorial7GlobalBoundary)
    add_dune_mpi_flags(tutorial8InteriorBoundary)
    add_dune_mpi_flags(tutorial9PeriodicBoundary)
    add_dune_parmetis_flags(tutorial1gettingstarted)
    add_dune_parmetis_flags(tutorial2traversal)
    add_dune_parmetis_flags(tutorial3visualisation)
    add_dune_parmetis_flags(tutorial4IntegrationGauss)
    add_dune_parmetis_flags(tutorial4IntegrationNormal)
    add_dune_parmetis_flags(tutorial5CommunicationConstant)
    add_dune_parmetis_flags(tutorial5CommunicationGlobalIndex)
    add_dune_parmetis_flags(tutorial6ParallelData)
    add_dune_parmetis_flags(tutorial7GlobalBoundary)
    add_dune_parmetis_flags(tutorial8InteriorBoundary)
    add_dune_parmetis_flags(tutorial9PeriodicBoundary)
endif(HAVE_MPI)

# Add boost
#if(Boost_FOUND)
#   target_link_libraries(tutorial1gettingstarted ${Boost_LIBRARIES})
#   target_link_libraries(tutorial2traversal ${Boost_LIBRARIES})
#   target_link_libraries(tutorial3visualisation ${Boost_LIBRARIES})
#   target_link_libraries(tutorial4IntegrationGauss ${Boost_LIBRARIES})
#   target_link_libraries(tutorial4IntegrationNormal ${Boost_LIBRARIES})
#   target_link_libraries(tutorial5CommunicationConstant ${Boost_LIBRARIES})
#   target_link_libraries(tutorial5CommunicationGlobalIndex ${Boost_LIBRARIES})
#   target_link_libraries(tutorial6ParallelData ${Boost_LIBRARIES})
#   target_link_libraries(tutorial7GlobalBoundary ${Boost_LIBRARIES})
#   target_link_libraries(tutorial8InteriorBoundary ${Boost_LIBRARIES})
#endif()

# set_property(TARGET tutorial1gettingstarted APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial2traversal APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial3visualisation APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial4IntegrationGauss APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial4IntegrationNormal APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial5CommunicationConstant APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial5CommunicationGlobalIndex APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial6ParallelData APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial7GlobalBoundary APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})
# set_property(TARGET tutorial8InteriorBoundary APPEND_STRING PROPERTY COMPILE_FLAGS ${CURVGRID_FLAGS})


