#!/bin/sh

# script to run configure
#
# Use Intel compilers
# Use ATLAS BLAS (compiled with gcc/g77)
#    This is why need FLIBS!
# Use my mpich library.
#

GCCLIB="-L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -lg2c"

#
# This works
#
./configure                 \
   --build=i686-pc-linux    \
   CC=/usr/local/intel/cc   \
   CFLAGS=" -O3"            \
   CXX=/usr/local/intel/c++ \
   CXXFLAGS=" -03"          \
   F77=/usr/local/intel/f77 \
   F77FLAGS=" -03"          \
   F90=/usr/local/intel/f90 \
   F90FLAGS=" -03"          \
   FLIBS="-L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -lg2c" \
   --prefix=${HOME}/Packages/install \
   --with-blas="-L${HOME}/cplane/misc/ATLAS/lib/Linux_PIIISSE1 -lcblas -lf77blas -latlas" \
   --enable-mpi \
   --with-mpi-incdir=${HOME}/include/mpi \
   --with-mpi-libdir="${HOME}/lib/mpi -lmpich" \
   --disable-default-packages \
   --enable-cmmlib
