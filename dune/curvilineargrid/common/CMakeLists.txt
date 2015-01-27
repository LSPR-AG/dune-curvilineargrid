set(HEADERS
  loggingmessage.hh
  vectorhelper.hh
  )

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/common)
