#!/bin/sh

# script to run configure for ICC

#
# mpiCC, etc compilers for HP cluster. These are wrappers for
# Intel X86 compilers. HP also defines MPI.
#
# BLAS is MKL library.
#

#
#
#
./configure                 \
   CC=mpicc \
   CXX=mpiCC \
   F77=mpif77 \
   F90=mpif90 \
   --with-blasinclude="-I/opt/intel/mkl/include -DUSE_MKL_CBLAS " \
   --with-blas=" -L/opt/intel/mkl/lib/32 -lguide -lmkl_ia32 " \
   --build=i686-pc-linux    \
   --prefix=${HOME}/Packages/install \
   --disable-default-packages \
   --enable-cmmlib


