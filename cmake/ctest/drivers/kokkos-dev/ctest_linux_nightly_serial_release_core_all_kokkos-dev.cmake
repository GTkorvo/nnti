# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.kokkos-dev.gcc.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_RELEASE_DEV_KokkosCore_all)
SET(CTEST_PARALLEL_LEVEL 4)
SET(CTEST_TEST_TYPE Experimental)
SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_PACKAGES KokkosCore KokkosLinAlg KokkosContainers KokkosAlgorithms KokkosExample)
SET(EXTRA_EXCLUDE_PACKAGES KokkosTSQR KokkosClassic)


SET(EXTRA_CONFIGURE_OPTIONS
  "-DTPL_ENABLE_HWLOC:BOOL=ON"
  "-DTPL_ENABLE_CUDA:BOOL=ON"
  "-DTPL_ENABLE_CUSPARSE:BOOL=ON"
  "-DTPL_ENABLE_Pthread:BOOL=ON"
  "-DKokkos_ENABLE_Pthread:BOOL=ON"
  "-DTrilinos_ENABLE_OpenMP:BOOL=ON"
  "-DKokkos_ENABLE_OpenMP:BOOL=ON"

#Turning qthreads off until where/how it will be installed can be determined
  "TPL_ENABLE_QTHREAD:BOOL=OFF"
  "QTHREAD_INCLUDE_DIRS:FILEPATH=path_to_qthread/include"
  "QTHREAD_LIBRARY_DIRS:FILEPATH=path_to_qthread/lib"

  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
    
  "-DTrilinos_ENABLE_KokkosCore=ON"
  "-DTrilinos_ENABLE_KokkosContainers=ON"
  "-DTrilinos_ENABLE_KokkosLinAlg=ON"
  "-DTrilinos_ENABLE_KokkosAlgorithms=ON"
  "-DTrilinos_ENABLE_KokkosExample=ON"
  "-DTrilinos_ENABLE_KokkosTSQR=OFF"
  "-DTrilinos_ENABLE_KokkosClassic=OFF"

  "-DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF"

)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
