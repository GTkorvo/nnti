# To generate/regenerate package.mk, run ./configure
#
# This file is the internal make fragment included by all the
# subdirectory Makefiles.
#

PACKAGE_BUGREPORT = mpsears@sandia.gov
PACKAGE_NAME = cmmlib
PACKAGE_STRING = cmmlib 1.0
PACKAGE_TARNAME = cmmlib
PACKAGE_VERSION = 1.0

#
# Source directories
#
abs_srcdir    = /home/mpsears/Projects/cmmlib/config.dir
abs_topsrcdir = /home/mpsears/Projects/cmmlib/config.dir
srcdir        = .
top_srcdir    = .


host          = i686-pc-linux-gnu
host_alias    = 
host_cpu      = i686
host_os       = linux-gnu
host_vendor   = pc

target        = i686-pc-linux-gnu
target_alias  = 
target_cpu    = i686
target_os     = linux-gnu
target_vendor = pc

build         = i686-pc-linux-gnu
build_alias   = i686-pc-linux
build_cpu     = i686
build_os      = linux-gnu
build_vendor  = pc

#
# Installation directories
#
prefix        = /home/mpsears/Packages/install
exec_prefix   = ${prefix}
bindir        = ${exec_prefix}/bin
libdir        = ${exec_prefix}/lib
libexecdir    = ${exec_prefix}/libexec
includedir    = ${prefix}/include
datadir       = ${prefix}/share

# install html documentation here
docdir        = /home/mpsears/Packages/install/doc

# i dont use these
infodir       = ${prefix}/info
mandir        = ${prefix}/man

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
