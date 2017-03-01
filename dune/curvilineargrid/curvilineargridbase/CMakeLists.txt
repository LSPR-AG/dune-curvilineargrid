
install(FILES
   curvilineargridbase.hh
   curvilineargridbasefactory.hh
   methods/curvilineargridbasecommunication.hh
   methods/curvilineargridbaseedge.hh
   methods/curvilineargridbaseentity.hh
   methods/curvilineargridbaseintersection.hh
   methods/curvilineargridbaseproperty.hh
   methods/curvilinearoctreewrapper.hh
   impl/curvilinearghostconstructor.hh
   impl/curvilinearglobalindexconstructor.hh
   impl/curvilineargridconstructor.hh
   impl/curvilineargridstorage.hh
   impl/curvilinearlooseoctree.hh
   impl/curvilinearoctant.hh
   impl/curvilinearoctreenode.hh
   impl/curvilinearperiodicconstructor.hh
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/curvilineargrid/curvilineargridbase)
