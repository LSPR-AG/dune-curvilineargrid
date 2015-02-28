add_subdirectory(curvilineargridbase)
add_subdirectory(curvilineargrid)
add_subdirectory(curvilineargridhowto)
add_subdirectory(io)
add_subdirectory(test)
add_subdirectory(common)
add_subdirectory(utility)

#install headers
install(FILES curvilineargrid.hh DESTINATION include/dune/curvilineargrid)

