#!/bin/sh

# script to run configure
#
# Generic gcc compilers
# Use ATLAS BLAS
#

#
#
#
./configure                 \
   FLIBS="-L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -lg2c" \
   --with-blasinclude="-I${HOME}/cplane/misc/ATLAS/include -DUSE_CBLAS" \
   --with-blas="-L${HOME}/cplane/misc/ATLAS/lib/Linux_PIIISSE1 -lcblas -lf77blas -latlas" \
   --build=i686-pc-linux    \
   --prefix=${HOME}/Packages/install \
   --disable-default-packages \
   --enable-cmmlib
