#!/bin/sh

# script to run configure
#
# Generic gcc compilers
# Use the internal reference blas
#

#
#
#
./configure                 \
   --build=i686-pc-linux    \
   --with-blasinclude="-DUSE_MYZGEMM" \
   --prefix=${HOME}/Packages/install \
   --disable-default-packages \
   --enable-cmmlib
