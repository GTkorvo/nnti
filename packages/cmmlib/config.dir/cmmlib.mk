# To generate/regenerate cmmlib.mk, run ./configure

#
# This file is a make fragment exported by cmmlib for
# use by dependent packages.
#

#
# Installation directories
#

# If you override bindir, etc this will not work properly,
# you will have to go back and munge the installed cmmlib.mk
#

cmmlib_prefix        = /home/mpsears/Packages/install
cmmlib_exec_prefix   = ${cmmlib_prefix}
cmmlib_bindir        = ${cmmlib_exec_prefix}/bin
cmmlib_libdir        = ${cmmlib_exec_prefix}/lib
cmmlib_libexecdir    = ${cmmlib_exec_prefix}/libexec
cmmlib_includedir    = ${cmmlib_prefix}/include
cmmlib_datadir       = ${cmmlib_prefix}/share

# HTML goes here
cmmlib_docdir        = ${cmmlib_prefix}/doc

# dont use these
cmmlib_infodir       = ${cmmlib_prefix}/info
cmmlib_mandir        = ${cmmlib_prefix}/man

#
# Source directories
#
cmmlib_abs_srcdir    = /home/mpsears/Projects/cmmlib/config.dir
cmmlib_abs_topsrcdir = /home/mpsears/Projects/cmmlib/config.dir
cmmlib_srcdir        = .
cmmlib_top_srcdir    = .
cmmlib_VPATH         = .

cmmlib_host          = i686-pc-linux-gnu
cmmlib_host_alias    = 
cmmlib_host_cpu      = i686
cmmlib_host_os       = linux-gnu
cmmlib_host_vendor   = pc

cmmlib_target        = i686-pc-linux-gnu
cmmlib_target_alias  = 
cmmlib_target_cpu    = i686
cmmlib_target_os     = linux-gnu
cmmlib_target_vendor = pc

cmmlib_build         = i686-pc-linux-gnu
cmmlib_build_alias   = i686-pc-linux
cmmlib_build_cpu     = i686
cmmlib_build_os      = linux-gnu
cmmlib_build_vendor  = pc


#
# Define compilers, etc from autoconf/configure
#
PATH_SEPARATOR = :
RANLIB = ranlib
SHELL = /bin/sh

CC=cc
CFLAGS=-g -O2
CPPFLAGS=

CXX = g++
CXXCPP = g++ -E
CXXFLAGS = -g -O2

MPI_CC_EXISTS = 
MPI_CC = @MPI_CC@

MPI_CXX_EXISTS = 
MPI_CXX = 

LIBS = 
FLIBS = -L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -lg2c

#
# Should I put BLAS_LIBS =  -L/home/mpsears/cplane/misc/ATLAS/lib/Linux_PIIISSE1 -lcblas -lf77blas -latlas  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -lg2c  ??
#
BLAS_LIBS = -L/home/mpsears/cplane/misc/ATLAS/lib/Linux_PIIISSE1 -lcblas -lf77blas -latlas
BLAS_INCLUDE = -I/home/mpsears/cplane/misc/ATLAS/include -DUSE_CBLAS

MPI_LIBS = @MPI_LIBS@
MPI_INCLUDE = @MPI_INCLUDE@

CMMLIB_LIBS = -L${cmmlib_libdir} -lcmm
CMMLIB_INCLUDE = -I${cmmlib_includedir}/cmmlib

#
# Try to figure out Fortran name mangling
#

sgemm = sgemm_

ifeq ("$(sgemm)", "sgemm")
 FORTRAN_MANGLE=-DFORTRAN_MANGLE_NONE
endif

ifeq ("$(sgemm)", "sgemm_")
 FORTRAN_MANGLE=-DFORTRAN_MANGLE_
endif

ifeq ("$(sgemm)", "SGEMM")
 FORTRAN_MANGLE=-DFORTRAN_MANGLE_UPPER
endif

ifeq ("$(sgemm)", "SGEMM_")
 FORTRAN_MANGLE=-DFORTRAN_MANGLE_UPPER_
endif
