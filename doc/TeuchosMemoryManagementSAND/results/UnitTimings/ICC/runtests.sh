#!/bin/bash

# This is a Release build with MS Visual C++ 2008 with:
#
#  -DTrilinos_ENABLE_BOOST
#  -DCMAKE_CXX_FLAG:STRING="... /D_SECURE_SCL=0 /DBOOST_SP_DISABLE_THREADS /Ox"
#

# to turn off the checked STL.

time ./Teuchos_RCP_PerformanceTests.exe \
  --details=ALL \
  --rel-cpu-speed=5e+3 \
  --max-array-size=100000 \
  2>&1 | tee  Teuchos_RCP_PerformanceTest.out

time ./Teuchos_Array_PerformanceTests.exe \
  --details=ALL \
  --rel-cpu-speed=5e+3 \
  --max-array-size=10000 \
  2>&1 | tee  Teuchos_Array_PerformanceTest.out
