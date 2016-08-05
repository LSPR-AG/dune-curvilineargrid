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
	tutorial6CommunicationConstant
	tutorial6CommunicationGlobalIndex
	tutorial7ParallelData
	tutorial8GlobalBoundary
	tutorial9InteriorBoundary
)

add_executable(tutorial1gettingstarted tutorial1gettingstarted.cc)
add_executable(tutorial2traversal tutorial2traversal.cc)
add_executable(tutorial3visualisation tutorial3visualisation.cc)
add_executable(tutorial4IntegrationGauss tutorial4IntegrationGauss.cc)
add_executable(tutorial4IntegrationNormal tutorial4IntegrationNormal.cc)
add_executable(tutorial6CommunicationConstant tutorial6CommunicationConstant.cc)
add_executable(tutorial6CommunicationGlobalIndex tutorial6CommunicationGlobalIndex.cc)
add_executable(tutorial7ParallelData tutorial7ParallelData.cc)
add_executable(tutorial8GlobalBoundary tutorial8GlobalBoundary.cc)
add_executable(tutorial9InteriorBoundary tutorial9InteriorBoundary.cc)
target_link_libraries(tutorial1gettingstarted dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial2traversal dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial3visualisation dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial4IntegrationGauss dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial4IntegrationNormal dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial6CommunicationConstant dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial6CommunicationGlobalIndex dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial7ParallelData dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial8GlobalBoundary dunecommon dunegeometry dunegrid) # dunecurvilineargeometry
target_link_libraries(tutorial9InteriorBoundary dunecommon dunegeometry dunegrid) # dunecurvilineargeometry


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
    add_dune_mpi_flags(tutorial6CommunicationConstant)
    add_dune_mpi_flags(tutorial6CommunicationGlobalIndex)
    add_dune_mpi_flags(tutorial7ParallelData)
    add_dune_mpi_flags(tutorial8GlobalBoundary)
    add_dune_mpi_flags(tutorial9InteriorBoundary)
    add_dune_parmetis_flags(tutorial1gettingstarted)
    add_dune_parmetis_flags(tutorial2traversal)
    add_dune_parmetis_flags(tutorial3visualisation)
    add_dune_parmetis_flags(tutorial4IntegrationGauss)
    add_dune_parmetis_flags(tutorial4IntegrationNormal)
    add_dune_parmetis_flags(tutorial6CommunicationConstant)
    add_dune_parmetis_flags(tutorial6CommunicationGlobalIndex)
    add_dune_parmetis_flags(tutorial7ParallelData)
    add_dune_parmetis_flags(tutorial8GlobalBoundary)
    add_dune_parmetis_flags(tutorial9InteriorBoundary)
endif(HAVE_MPI)

# Add boost
if(Boost_FOUND)
   target_link_libraries(tutorial1gettingstarted ${Boost_LIBRARIES})
   target_link_libraries(tutorial2traversal ${Boost_LIBRARIES})
   target_link_libraries(tutorial3visualisation ${Boost_LIBRARIES})
   target_link_libraries(tutorial4IntegrationGauss ${Boost_LIBRARIES})
   target_link_libraries(tutorial4IntegrationNormal ${Boost_LIBRARIES})
   target_link_libraries(tutorial6CommunicationConstant ${Boost_LIBRARIES})
   target_link_libraries(tutorial6CommunicationGlobalIndex ${Boost_LIBRARIES})
   target_link_libraries(tutorial7ParallelData ${Boost_LIBRARIES})
   target_link_libraries(tutorial8GlobalBoundary ${Boost_LIBRARIES})
   target_link_libraries(tutorial9InteriorBoundary ${Boost_LIBRARIES})
endif()
