#!/bin/bash

#
# Checkin-test script for Trios in the preCopyrightTrilinos repository. 
#
# For macos with OpenMPI installed in /opt/local.
#
# We need to make sure Trios configures and builds with and without
# an RDMA transport library.   On some machines, this won't be 
# available. 
#
# The Trios checkin test only evaluates the MPI_DEBUG build.  Do we 
# need to have a SERIAL_RELEASE build?
# 

# Use this to pass in addtitional args (e.g., --push)
EXTRA_ARGS=$@

# Set this point to your source
TRILINOS_SRC_DIR=${HOME}/research/workspace/Trilinos

echo "
-D BUILD_SHARED:BOOL=ON
-D TPL_ENABLE_Portals:BOOL=ON
-D Portals_INCLUDE_DIRS:PATH=${HOME}/research/support/lib/portals/include
-D Portals_LIBRARY_DIRS:PATH=${HOME}/research/support/lib/portals/lib
" > COMMON.config

echo "
-D TPL_ENABLE_MPI:BOOL=ON
-D MPI_BASE_DIR:PATH=/opt/local/lib/openmpi
" > MPI_DEBUG.config


# Test the MPI_DEBUG.  This script will do pull, configure, build, test. 
${TRILINOS_SRC_DIR}/checkin-test.py \
   --no-eg-git-version-check \
   --trilinos-src-dir=${TRILINOS_SRC_DIR} \
   --enable-all-packages=off \
   --no-enable-fwd-packages \
   --without-serial-release \
   --extra-repos=preCopyrightTrilinos \
   --enable-packages=Trios \
   --extra-build=MPI_DEBUG \
   --pull \
   --configure \
   --build \
   --test \
   $EXTRA_ARGS

