set(HEADERS
  curvilineargridstorage.hh
  curvilineargridconstructor.hh
  curvilineargridbase.hh
  curvilinearlooseoctree.hh
  curvilinearoctant.hh
  curvilinearoctreenode.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/curvilineargridbase)
