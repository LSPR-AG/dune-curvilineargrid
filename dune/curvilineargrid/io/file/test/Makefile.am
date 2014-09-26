ALLTESTS = vtktest gnuplottest vtksequencetest subsamplingvtktest gmshtest curvilineargmshtest

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"


# programs just to build when "make check" is used
check_PROGRAMS = $(ALLTESTS)

# list of tests to run
TESTS = $(ALLTESTS) mpivtktest

ALLTESTS += conformvolumevtktest
conformvolumevtktest_SOURCES = conformvolumevtktest.cc

ALLTESTS += nonconformboundaryvtktest
nonconformboundaryvtktest_SOURCES = nonconformboundaryvtktest.cc

vtktest_SOURCES = vtktest.cc
vtktest_CPPFLAGS = $(AM_CPPFLAGS)		\
	$(DUNEMPICPPFLAGS)
vtktest_LDFLAGS = $(AM_LDFLAGS)			\
	$(DUNEMPILDFLAGS)
vtktest_LDADD =					\
	$(DUNEMPILIBS)				\
	$(LDADD)

vtksequencetest_SOURCES = vtksequencetest.cc

gnuplottest_SOURCES = gnuplottest.cc

subsamplingvtktest_SOURCES = subsamplingvtktest.cc test-linking.cc
subsamplingvtktest_CPPFLAGS = $(AM_CPPFLAGS)	\
	$(DUNEMPICPPFLAGS)
subsamplingvtktest_LDFLAGS = $(AM_LDFLAGS)	\
	$(DUNEMPILDFLAGS)
subsamplingvtktest_LDADD =			\
	$(DUNEMPILIBS)				\
	$(LDADD)

gmshtest_SOURCES  = gmshtest.cc
gmshtest_CPPFLAGS = $(AM_CPPFLAGS)  $(ALL_PKG_CPPFLAGS)
gmshtest_LDFLAGS  = $(AM_LDFLAGS)   $(ALL_PKG_LDFLAGS)
gmshtest_LDADD    = $(ALL_PKG_LIBS) $(LDADD)

gmshtest_alberta2d_SOURCES  = gmshtest.cc
gmshtest_alberta2d_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS) $(GRAPE_CPPFLAGS) -DGRIDDIM=$(ALBERTA_DIM)
gmshtest_alberta2d_LDFLAGS  = $(AM_LDFLAGS)  $(ALBERTA_LDFLAGS)  $(GRAPE_LDFLAGS)
gmshtest_alberta2d_LDADD    =	$(GRAPE_LIBS)	 $(ALBERTA_LIBS)		 $(LDADD)

gmshtest_alberta3d_SOURCES  = gmshtest.cc
gmshtest_alberta3d_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA3D_CPPFLAGS) $(GRAPE_CPPFLAGS)
gmshtest_alberta3d_LDFLAGS  = $(AM_LDFLAGS)  $(ALBERTA3D_LDFLAGS)  $(GRAPE_LDFLAGS)
gmshtest_alberta3d_LDADD    =	$(GRAPE_LIBS)	 $(ALBERTA3D_LIBS)		 $(LDADD)

gmshtest_alugrid_SOURCES  = gmshtest.cc
gmshtest_alugrid_CPPFLAGS = $(AM_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(GRAPE_CPPFLAGS) $(GRIDDIM_CPPFLAGS)
gmshtest_alugrid_LDFLAGS  = $(AM_LDFLAGS)	 $(ALUGRID_LDFLAGS)	 $(GRAPE_LDFLAGS)
gmshtest_alugrid_LDADD    =	$(GRAPE_LIBS)	 $(ALUGRID_LIBS)		 $(LDADD)

curvilineargmshtest_SOURCES  = curvilineargmshtest.cc
curvilineargmshtest_CPPFLAGS = $(AM_CPPFLAGS)  $(ALL_PKG_CPPFLAGS) -static
curvilineargmshtest_LDFLAGS  = $(AM_LDFLAGS)   $(ALL_PKG_LDFLAGS)
curvilineargmshtest_LDADD    = $(ALL_PKG_LIBS) $(LDADD)


include $(top_srcdir)/am/global-rules

CLEANFILES = *.vtu *.vtp *.data sgrid*.am *.pvtu *.pvtp *.pvd

EXTRA_DIST = CMakeLists.txt
