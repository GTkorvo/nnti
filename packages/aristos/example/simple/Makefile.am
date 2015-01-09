# @HEADER
# ************************************************************************
#
#                     Aristos Optimization Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Denis Ridzal (dridzal@sandia.gov)
#
# ************************************************************************
# @HEADER

include $(top_builddir)/Makefile.export.aristos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ARISTOS_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ARISTOS_TEST_INCLUDES))
else
EXPORT_LIBS = $(ARISTOS_TEST_LIBS)
EXPORT_INC_PATH = $(ARISTOS_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)  

EXEEXT = .exe

noinst_PROGRAMS = testSimple

testSimple_SOURCES = \
			$(srcdir)/testSimple.cpp \
			$(srcdir)/SimpleApp_SimpleEpetraLagMult.cpp \
			$(srcdir)/SimpleApp_SimpleEpetraFeasStep.cpp \
			$(srcdir)/SimpleApp_SimpleEpetraConstraints.cpp \
			$(srcdir)/SimpleApp_SimpleEpetraObjective.cpp \
			$(srcdir)/SimpleApp_SimpleEpetraHessVec.cpp \
			$(srcdir)/SimpleApp_SimpleEpetraDataPool.cpp

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
testSimple_DEPENDENCIES = \
	$(top_builddir)/src/libaristos.a \
	$(top_builddir)/../amesos/src/libamesos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

testSimple_CXXFLAGS = -I$(top_srcdir)/src -I$(srcdir)

#np# List libraries that np is dependent on.
testSimple_LDADD = $(EXPORT_LIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 

# Dummy for new export makefiles
LIBS = 
