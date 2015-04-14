set(HEADERS
  griddiagnostic.hh
  gridbasediagnostic.hh
  allcommunication.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/utility)
