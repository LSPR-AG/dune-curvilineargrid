add_subdirectory(curvilineargridbase)
add_subdirectory(curvilineargrid)
add_subdirectory(io)
add_subdirectory(test)
add_subdirectory(feedback)

#install headers
install(FILES curvilineargrid.hh DESTINATION include/dune/curvilineargrid)

