#!/bin/sh

# script to run configure on standard linux system
#
# Use generic gcc compilers
# Lets see if it can find a blas on its own.
#

#
# Note the use of CFLAGS=... 
# Without this flag, the configure script will
# specify -g -O2.
#
# --with-blasinclude=... specifies that we
# want to compile with the Fortran blas found
# by the configure script.
#
# You can prepend any flags you want to CFLAGS
# with the --with-cflags=... argument.
#

./configure                 \
   CFLAGS="-O3"             \
   --with-blasinclude="-DUSE_FBLAS" \
   --prefix=${HOME}/Packages/install \
   --disable-default-packages \
   --enable-cmmlib

