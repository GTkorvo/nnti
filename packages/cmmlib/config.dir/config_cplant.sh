#!/bin/sh

# script to run configure
#
# Use Cplant (OSF) compilers and MPI
# Use CXML blas
#

#
# This works BUT you must go into cmmlib.mk afterwards
# and delete the -std1 flag!!
#
./configure                 \
   CPPFLAGS=-D__USE_STD_IOSTREAM \
   CXX="/usr/local/cplant/ross/current/bin/c++" \
   CC="/usr/local/cplant/ross/current/bin/cc" \
   F77="/usr/local/cplant/ross/current/bin/f77" \
   F90="/usr/local/cplant/ross/current/bin/f90" \
   --host=alpha-unknown-linux \
   --with-blasinclude="-DUSE_FBLAS" \
   --with-blas=" -lcxml_ev6" \
   --prefix=${HOME}/Packages/install \
   --disable-default-packages \
   --enable-cmmlib

