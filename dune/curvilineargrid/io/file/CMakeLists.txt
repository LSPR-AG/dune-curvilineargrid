add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(vtk)
add_subdirectory(gnuplot)

set(HEADERS
  gmshreader.hh
  gmshwriter.hh
  curvilineargmshreader.hh
  curvilinearvtkreader.hh
  gnuplot.hh
  vtk.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file)
