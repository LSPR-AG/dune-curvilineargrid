#add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(test)

install(FILES
   curvilineargmshreader.hh
   curvilinearvtkentitysubset.hh
   curvilinearvtkformat.hh
   curvilinearvtkgridwriter.hh
   curvilinearvtkwriter.hh
   gmsh2dunemapper.hh
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/io/file)
