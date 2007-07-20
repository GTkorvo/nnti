HEADER    = ./
#
# Include Trilinos header that defines pertinent compilation information
# Substitute these for the appropriate locations
#   TRILINOS_HOME = /home/u/cgbaker/Trilinos
#   TRILINOS_ARCH = JETER_MPI
#

TRILINOS_HOME = /scratch/Rel7
# TRILINOS_ARCH = MPI
TRILINOS_ARCH = TOMS
include ${TRILINOS_HOME}/${TRILINOS_ARCH}/include/Makefile.export.anasazi.macros
#
# Linear algebra object libraries
#
OLIBS = -L$(ANASAZI_prefix)/lib -L. -lanasazi -lepetra -lteuchos $(ANASAZI_LDFLAGS) $(ANASAZI_BLAS_LIBS) $(ANASAZI_LAPACK_LIBS) $(ANASAZI_LIBS)

# FNIT = -fno-implicit-templates
FNIT = 

# EO = explicits.o
EO = 

#
# Compilation rules
#
.c.o:
	$(ANASAZI_CC) $(ANASAZI_CFLAGS) -c $(ANASAZI_DEFS) $<

.cpp.o:
	$(ANASAZI_CXX) $(ANASAZI_CPPFLAGS) $(ANASAZI_CXXFLAGS) -c $(ANASAZI_DEFS) -I. -I$(ANASAZI_prefix)/include ${FNIT} $<

all:    BKS_2DQ1.exe

explicits.o: explicits.cpp
	$(ANASAZI_CXX) $(ANASAZI_CPPFLAGS) $(ANASAZI_CXXFLAGS) -c $(ANASAZI_DEFS) -I$(ANASAZI_prefix)/include $<

BKS_2DQ1.exe:  BKS_2DQ1.o $(EO)
	$(ANASAZI_CXX) $(ANASAZI_CXXFLAGS) -o $@ $^ $(OLIBS) $(ANASAZI_FLIBS) $(EO)

clean: 
	/bin/rm -f *.o core

cleanall:
	/bin/rm -f *.o *.exe core libmodal.a
