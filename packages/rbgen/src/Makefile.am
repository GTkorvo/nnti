# @HEADER
# ************************************************************************
# 
#                           RBGen Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package rbgen subdirectory src built successfully."
	@echo ""

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside RBGen. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.
CORE_H = \
        $(srcdir)/RBGen_ConfigDefs.h \
        $(srcdir)/RBGen_FileIOFactory.hpp \
        $(srcdir)/RBGen_FileIOHandler.hpp \
        $(srcdir)/RBGen_Filter.hpp \
        $(srcdir)/RBGen_MethodFactory.hpp \
        $(srcdir)/RBGen_Method.hpp \
        $(srcdir)/RBGen_NoPreprocessor.hpp \
        $(srcdir)/RBGen_Params.h \
        $(srcdir)/RBGen_PODMethod.hpp \
        $(srcdir)/RBGen_PreprocessorFactory.hpp \
        $(srcdir)/RBGen_Preprocessor.hpp


CORE = \
        $(srcdir)/RBGen_Params.cpp

#-----------------------------------------------------

EXTRA_EPETRA_H = \
        $(srcdir)/RBGen_BurkardtFileIOHandler.h \
        $(srcdir)/RBGen_EpetraMVFileIOFactory.h \
        $(srcdir)/RBGen_EpetraMVMethodFactory.h \
        $(srcdir)/RBGen_EpetraMVPreprocessorFactory.h \
        $(srcdir)/RBGen_LapackPOD.h \
        $(srcdir)/RBGen_MSPreprocessor.h \
        $(srcdir)/RBGen_Utils.h

EXTRA_EPETRA = \
        $(srcdir)/RBGen_BurkardtFileIOHandler.cpp \
        $(srcdir)/RBGen_EpetraMVFileIOFactory.cpp \
        $(srcdir)/RBGen_EpetraMVMethodFactory.cpp \
        $(srcdir)/RBGen_EpetraMVPreprocessorFactory.cpp \
        $(srcdir)/RBGen_LapackPOD.cpp \
        $(srcdir)/RBGen_MSPreprocessor.cpp \
        $(srcdir)/RBGen_Utils.cpp
	
if BUILD_RBGEN_EPETRA

EPETRA_H = $(EXTRA_EXTENDED_H)
EPETRA = $(EXTRA_EPETRA)

else 

EPETRA_H =
EPETRA =

endif

#-----------------------------------------------------

EXTRA_EPETRAEXT_H = \
        $(srcdir)/RBGen_EpetraCrsMatrixFileIOHandler.h \
        $(srcdir)/RBGen_MatrixMarketFileIOHandler.h

EXTRA_EPETRAEXT = \
        $(srcdir)/RBGen_EpetraCrsMatrixFileIOHandler.cpp \
        $(srcdir)/RBGen_MatrixMarketFileIOHandler.cpp

if BUILD_RBGEN_EPETRAEXT

EPETRAEXT_H = $(EXTRA_EPETRAEXT_H)
EPETRAEXT = $(EXTRA_EPETRAEXT)

else

EPETRAEXT_H =
EPETRAEXT = 

endif

#-----------------------------------------------------

EXTRA_ANASAZI_EPETRA_H = \
        $(srcdir)/RBGen_AnasaziPOD.h \
        $(srcdir)/RBGen_IncSVDPOD.h \
		$(srcdir)/RBGen_ISVDSingle.h \
		$(srcdir)/RBGen_ISVDMultiCD.h \
		$(srcdir)/RBGen_ISVDUDV.h \
		$(srcdir)/RBGen_ISVD_SingleUDV.h 

EXTRA_ANASAZI_EPETRA = \
        $(srcdir)/RBGen_AnasaziPOD.cpp \
        $(srcdir)/RBGen_IncSVDPOD.cpp \
		$(srcdir)/RBGen_ISVDSingle.cpp \
		$(srcdir)/RBGen_ISVDMultiCD.cpp \
		$(srcdir)/RBGen_ISVDUDV.cpp \
		$(srcdir)/RBGen_ISVD_SingleUDV.cpp 

if BUILD_RBGEN_ANASAZI


ANASAZI_EPETRA_H = $(EXTRA_ANASAZI_EPETRA_H)
ANASAZI_EPETRA = $(EXTRA_ANASAZI_EPETRA)

else

ANASAZI_EPETRA_H = 
ANASAZI_EPETRA = 

endif

#-----------------------------------------------------

EXTRA_NETCDF_H = \
        $(srcdir)/RBGen_NetCDFFileIOHandler.h

EXTRA_NETCDF = \
        $(srcdir)/RBGen_NetCDFFileIOHandler.cpp

if BUILD_RBGEN_NETCDF
if BUILD_RBGEN_EPETRA

NETCDF_EPETRA_H = $(EXTRA_NETCDF_H)
NETCDF_EPETRA = $(EXTRA_NETCDF)

else

NETCDF_EPETRA_H = 
NETCDF_EPETRA = 

endif

else

NETCDF_EPETRA_H = 
NETCDF_EPETRA = 

endif

#-----------------------------------------------------

UTIL_H = \
	$(srcdir)/RBGen_Version.h

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


#np# List locations of include files within Trilinos but outside of this 
#np# package
#  
include $(top_builddir)/Makefile.export.rbgen

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(RBGEN_INCLUDES))
else
EXPORT_INCLUDES = $(RBGEN_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace rbgen with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# RBGen library specifications
# ------------------------------------------------------------------------
#np# replace rbgen with the name of the package being autotool'ed here
RBGEN_LIB = librbgen.a

#np# replace rbgen with the name of the package being autotool'ed here
RBGEN_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(EPETRA_H) \
	$(EPETRAEXT_H) \
	$(ANASAZI_EPETRA_H) \
	$(NETCDF_EPETRA_H)	

#np# replace rbgen with the name of the package being autotool'ed here
librbgen_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(EPETRA) \
	$(EPETRAEXT) \
	$(ANASAZI_EPETRA) \
	$(NETCDF_EPETRA)	

#np# replace rbgen with the name of the package being autotool'ed here
EXTRA_librbgen_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace rbgen with the name of the package being autotool'ed here
librbgen_a_AR = $(ALTERNATE_AR)
else

#np# replace rbgen with the name of the package being autotool'ed here
librbgen_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace rbgen with the name of the package being autotool'ed here
librbgen_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace rbgen with the name of the package being autotool'ed here
lib_LIBRARIES = $(RBGEN_LIB)

# ------------------------------------------------------------------------
# List of all drivers to build
# ------------------------------------------------------------------------

if BUILD_RBGEN_DRIVERS
if BUILD_RBGEN_EPETRA
bin_PROGRAMS = RBGenDriver_EpetraMV
RBGenDriver_EpetraMV_SOURCES = $(srcdir)/RBGenDriver_EpetraMV.cpp
RBGenDriver_EpetraMV_LDADD = $(RBGEN_LIB) $(RBGEN_LIBS)
RBGenDriver_EpetraMV_LDFLAGS = 
endif
endif

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace rbgen with the name of the package being autotool'ed here
include_HEADERS = $(RBGEN_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace rbgen with the name of the package being autotool'ed here
# SPECIAL NOTE: RBGen_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace rbgen with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = RBGen_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






