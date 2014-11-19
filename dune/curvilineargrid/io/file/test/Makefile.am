ALLTESTS = curvilineargmshtest

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"


# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)

# list of tests to run
TESTS = $(ALLTESTS)

curvilineargmshtest_SOURCES  = curvilineargmshtest.cc
curvilineargmshtest_CPPFLAGS = $(AM_CPPFLAGS)  $(ALL_PKG_CPPFLAGS) -static
curvilineargmshtest_LDFLAGS  = $(AM_LDFLAGS)   $(ALL_PKG_LDFLAGS)
curvilineargmshtest_LDADD    = $(ALL_PKG_LIBS) $(LDADD)


include $(top_srcdir)/am/global-rules

CLEANFILES = *.vtu *.vtp *.data *.pvtu *.pvtp *.pvd

EXTRA_DIST = CMakeLists.txt
