HEADER    = ./
#
# Include Trilinos header that defines pertinent compilation information
# Substitute these for the appropriate locations
#   TRILINOS_HOME = /home/u/cgbaker/Trilinos
#   TRILINOS_ARCH = JETER_MPI
#

TRILINOS_HOME = /Users/cgbaker/Trilinos
# TRILINOS_ARCH = MPI
TRILINOS_ARCH = SERIAL
include ${TRILINOS_HOME}/${TRILINOS_ARCH}/include/Makefile.export.anasazi.macros
#
# Linear algebra object libraries
#
OLIBS = -L$(TRILINOS_HOME)/BUILD_${TRILINOS_ARCH}/packages/anasazi/util/ModeLaplace -L$(anasazi_prefix)/lib -L. -lModeLaplace -lanasazi -lepetra -lepetrathyra -lthyra -lteuchos $(anasazi_LDFLAGS) $(anasazi_BLAS_LIBS) $(anasazi_LAPACK_LIBS) $(anasazi_LIBS)

#
# Compilation rules
#
.c.o:
	$(anasazi_CC) $(anasazi_CFLAGS) -c $(anasazi_DEFS) $<

.cpp.o:
	$(anasazi_CXX) $(anasazi_CPPFLAGS) $(anasazi_CXXFLAGS) -c $(anasazi_DEFS) -I$(TRILINOS_HOME)/packages/anasazi/util/ModeLaplace -I$(anasazi_prefix)/include $<

all: BlockDavidson1DQ1_e.exe BlockDavidson1DQ1_t.exe BlockDavidson1DQ1_te.exe \
	LOBPCG1DQ1_e.exe LOBPCG1DQ1_t.exe LOBPCG1DQ1_te.exe

BlockDavidson1DQ1_e.exe:  BlockDavidson1DQ1_e.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

BlockDavidson1DQ1_t.exe:  BlockDavidson1DQ1_t.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

BlockDavidson1DQ1_te.exe:  BlockDavidson1DQ1_te.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

LOBPCG1DQ1_e.exe:  LOBPCG1DQ1_e.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

LOBPCG1DQ1_t.exe:  LOBPCG1DQ1_t.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

LOBPCG1DQ1_te.exe:  LOBPCG1DQ1_te.o
	$(anasazi_CXX) $(anasazi_CXXFLAGS) -o $@ $^ $(OLIBS) $(anasazi_FLIBS)

clean: 
	/bin/rm -f *.o core

cleanall:
	/bin/rm -f *.o *.exe core libmodal.a
