
#
# The default target is help
#
help:
	@echo "Usage: make target, where target is"
	@echo "   help       -- print this message"
	@echo "   configure  -- configure package (??)"
	@echo "   clean      -- clean package"
	@echo "   all        -- compile package"
	@echo "   install    -- install package"
	@echo " remaining targets are for developers only"
	@echo "   dox        -- generate online documentation"
	@echo "   autoconf   -- regenerate configure script"
	@echo "   uninstall  -- remove installation"
	@echo "   test       -- build tests"
	@echo "   example    -- build examples"
	@echo "   dist       -- build distribution"

#
# Include (if we can) the current configuration.
# If not available we won't fail but some targets
# will not be available.
#
-include config.dir/package.mk

#
# Run the configure script. You can also do this by
# hand, with many options. See also the scripts in
# config.dir for various systems. Run configure --help
# for info.
#
configure:
	./config.dir/configure

#
# Run autoconf to generate the configure script. This should normally
# NOT be necessary.
#
autoconf:
	(cd config.dir; autoconf;)

#
# Build documentation
#   Use version 1.3.7 of doxgyen or better, earlier has problems
#   with formulas (will function but may produce erroneous output).
#
dox:
	@echo "Build online documentation"
	doxygen config.dir/config.dox

#
# Clean
#
# (do not remove Makefiles or .mk files)
clean:
	@echo "Clean cmmlib"
	(cd src/base; make clean)
	(cd src/rand; make clean)
	(cd src/geom3; make clean)
	(cd src/blas; make clean)
	(cd src/matrix; make clean)
	(cd lib; rm *.a);
	(cd include; rm *.h);
	(cd config.dir; rm -rf autom4te.cache config.log config.status)

#
# Scrub. After this you must run autoconf to
# regenerate configure and .mk files
scrub: clean
	@echo "Scrub cmmlib"
	(cd config.dir; rm -f configure *.mk)

#*********************************************************
#
# All of the following targets will fail if the package
# is not configured!
#
#*********************************************************

#
# Compile the package
#
all: config.dir/package.mk
	@echo "Build cmmlib"
	(cd src/base; make all;)
	(cd src/rand; make all;)
	(cd src/geom3; make all;)
	(cd src/blas; make all;)
	(cd src/matrix; make all;)

#
# Install the package
#
install: config.dir/package.mk
	@echo "Install cmmlib"
	mkdir -p ${includedir}
	mkdir -p ${libdir}
	mkdir -p ${docdir}
	install -d ${includedir}/cmmlib
	install include/*.h config.dir/cmmlib.mk ${includedir}/cmmlib
	install lib/* ${libdir}
	install -d  ${docdir}/cmmlib
	install doc/dox/html/*  ${docdir}/cmmlib


#
# Uninstall the package
#
uninstall: config.dir/package.mk
	@echo "Uninstall cmmlib"
	rm -rf ${includedir}/cmmlib
	rm -f ${libdir}/libcmm.a
	rm -rf ${docdir}/cmmlib

#
# Build test programs
#
test: config.dir/package.mk
	@echo "Build cmmlib test programs"

#
# Build examples
#
examples: config.dir/package.mk
	@echo "Build cmmlib examples"

#
# Build a distribution. This works by building a copy of the source tree in a temp directory
# then tar-gzip it.
#
DISTTEMP=..

VERSION=$(PACKAGE_VERSION)
DISTDIR=$(PACKAGE_NAME)-$(VERSION)
DISTTAR=$(DISTDIR).tar
DISTTARGZ=$(DISTDIR).tar.gz

SRCDIR=$(abs_srcdir)/..

dist: config.dir/package.mk
	@echo "Build cmmlib distribution $(DISTDIR) -> $(DISTTAR).gz"
	(j=`pwd`; cd $(DISTTEMP); rm -rf $(DISTDIR) $(DISTTARGZ); mkdir $(DISTDIR); cp -r $(SRCDIR)/* $(DISTDIR); \
            tar cf $(DISTTAR) $(DISTDIR); gzip $(DISTTAR); rm -rf $(DISTDIR); )
	ls -al ../$(DISTTARGZ)

PHONY: all test examples clean scrub install uninstall dist doc

################ Debugging targets ###########

test1: test1.o
	$(CC) $(CCFLAGS) -o test1 test1.o $(LIBS)

#
#
#

#
# Debug config
#
echo:
	@echo "srcdir $(srcdir)"
	@echo "topdir $(top_srcdir)"
	@echo "CC     $(CC) $(CCFLAGS)"
	@echo "CXX    $(CXX) $(CXXFLAGS)"
	@echo "AR     $(AR)"
	@echo "RANLIB $(RANLIB)"
	@echo "MPI_CC_EXISTS $(MPI_CC_EXISTS)"
	@echo "MPI_CC     $(MPI_CC) $(CCFLAGS)"
	@echo "MPI_CXX    $(MPI_CXX) $(CXXFLAGS)"

