set(modules
DuneCurvilineargrid.cmake
DuneCurvilineargridMacros.cmake
AddMETISFlags.cmake
AddParMETISFlags.cmake
)

install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
