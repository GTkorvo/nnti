#!/bin/sh

# script to run configure
#
# Use Cplant (Compaq) compilers and MPI
# Use GOTO blas
#

#
# This works BUT you must go into cmmlib.mk and package.mk files
# afterwards and delete the -std1 flag!! NOT MY FAULT! I think this
# is due to autoconf looking at an older version of the
# Compaq compiler.
#
./configure                 \
   CPPFLAGS=-D__USE_STD_IOSTREAM \
   CXX="/usr/local/cplant/ross/current/bin/c++" \
   CC="/usr/local/cplant/ross/current/bin/cc" \
   F77="/usr/local/cplant/ross/current/bin/f77" \
   F90="/usr/local/cplant/ross/current/bin/f90" \
   --host=alpha-unknown-linux \
   --with-blasinclude="-DUSE_FBLAS" \
   --with-blas=" -L/home/mpsears/packages/build/goto-blas-alpha/blas-alpha -lstatabo_ev6 -lfor -lFutil" \
   --prefix=${HOME}/Packages/install \
   --disable-default-packages \
   --enable-cmmlib


