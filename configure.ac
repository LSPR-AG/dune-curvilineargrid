# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.62])
DUNE_AC_INIT # gets module version from dune.module file
AC_CONFIG_SRCDIR([src/dune_curvilineargrid.cc])
AC_CONFIG_HEADERS([config.h])

# --enable-extra-utils
DUNE_EXTRA_UTILS

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-geometry dune-grid
# this implies checking for [dune-common], [dune-geometry], [dune-grid]
DUNE_CHECK_ALL

# set up flags for the automated test system
DUNE_AUTOBUILD_FLAGS

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

DUNE_SUMMARY_ALL

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  cmake/Makefile
  cmake/modules/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/curvilineargrid/Makefile
  m4/Makefile
  dune/curvilineargrid/test/Makefile
  dune/curvilineargrid/io/file/test/Makefile
  dune/curvilineargrid/io/file/gnuplot/Makefile
  dune/curvilineargrid/io/file/vtk/Makefile
  dune/curvilineargrid/io/file/Makefile
  dune/curvilineargrid/io/Makefile  
  dune-curvilineargrid.pc
])
AC_CONFIG_FILES([dune/curvilineargrid/io/file/test/mpivtktest],
    [chmod +x dune/curvilineargrid/io/file/test/mpivtktest])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
