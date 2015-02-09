set(HEADERS
  curvilineargridstorage.hh
  curvilineargridconstructor.hh
  curvilineargridbase.hh
  curvilinearghostconstructor.hh
  curvilinearpostconstructor.hh
  curvilinearlooseoctree.hh
  curvilinearoctant.hh
  curvilinearoctreenode.hh
  curvilineargridbasefactory.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/curvilineargridbase)
