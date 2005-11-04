#
# Include Trilinos header that defines pertinent compilation information
# Substitute these for the appropriate locations
#
TRILINOS_HOME     = /Users/cgbaker/Trilinos
TRILINOS_BUILD    = /Users/cgbaker/Trilinos/BUILD_SERIAL

# Also, we need the libModeLaplacde.a library, in the anasazi/util/ModeLaplace directory
# This should be built if configured with --enable-anasazi-tests or --enable-tests
# If these options were not configured, then the library can be built manually:
#    (after building anasazi)
#    cd ${BUILD_DIR}
#    cd packages/anasazi/util/ModeLaplace
#    make

include Makefile.export.anasazi.macros # This will be built automatically!
include ${TRILINOS_BUILD}/packages/anasazi/Makefile.export.anasazi # Use the code in the build directory!

#
# Linear algebra object libraries
#
OLIBS = -L$(TRILINOS_BUILD)/packages/anasazi/util/ModeLaplace -lModeLaplace $(ANASAZI_LIBS)

#
# Compilation rules
#
.c.o:
	$(anasazi_CC) $(anasazi_CFLAGS) $(ANASAZI_INCLUDES) -c $(anasazi_DEFS) $<

.cpp.o:
	$(anasazi_CXX) $(anasazi_CPPFLAGS) $(anasazi_CXXFLAGS) $(ANASAZI_INCLUDES) -c $(anasazi_DEFS) \
  -I$(TRILINOS_HOME)/packages/anasazi/util/ModeLaplace -I$(anasazi_prefix)/include $<

all: BlockDavidson1DQ1_e.exe BlockDavidson1DQ1_t.exe BlockDavidson1DQ1_te.exe  \
	LOBPCG1DQ1_e.exe LOBPCG1DQ1_t.exe LOBPCG1DQ1_te.exe 

LOBPCG1DQ1_e.exe:  LOBPCG1DQ1_e.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

LOBPCG1DQ1_t.exe:  LOBPCG1DQ1_t.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

LOBPCG1DQ1_te.exe:  LOBPCG1DQ1_te.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

BlockDavidson1DQ1_e.exe:  BlockDavidson1DQ1_e.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

BlockDavidson1DQ1_t.exe:  BlockDavidson1DQ1_t.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

BlockDavidson1DQ1_te.exe:  BlockDavidson1DQ1_te.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

clean: 
	/bin/rm -f *.o core

cleanall:
	/bin/rm -f *.o *.exe core libmodal.a

top_srcdir = $(TRILINOS_HOME)/packages/anasazi
top_builddir = $(TRILINOS_BUILD)/packages/anasazi

# This like was taken directly out of anasazi/Makefile.am!
Makefile.export.anasazi.macros:
	perl $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile anasazi > Makefile.export.anasazi.macros
