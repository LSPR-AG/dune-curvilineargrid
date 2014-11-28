#add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(test)

set(HEADERS
  curvilineargmshreader.hh
  curvilinearvtkwriter.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/io/file)
