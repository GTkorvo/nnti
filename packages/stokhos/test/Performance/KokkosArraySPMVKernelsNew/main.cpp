// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include <string>
#include <iostream>
#include <cstdlib>

#include "Kokkos_Threads.hpp"
#include "Kokkos_hwloc.hpp"

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_OPENMP) && defined(HAVE_STOKHOS_MKL)
#include <omp.h>
#endif

template <typename scalar>
int mainHost(bool test_flat, bool test_orig, bool test_deg, bool test_lin,
             bool test_block, bool symmetric, bool mkl);

template <typename scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_lin,
             bool test_block, bool symmetric, int device_id);

int main(int argc, char *argv[])
{
  // Defaults
  bool test_host = true;
  bool test_cuda = true;
  bool test_block = true;
  bool test_flat = true;
  bool test_orig = true;
  bool test_deg = false;
  bool test_lin = false;
  bool symmetric = true;
  bool single = false;
  bool mkl = false;
  int device = 0;

  // Parse command line arguments
  bool print_usage = false;
  int i=1;
  while (i<argc) {
    std::string s(argv[i]);
    if (s == "cuda")
      test_cuda = true;
    else if (s == "no-cuda")
      test_cuda = false;
    else if (s == "host")
      test_host = true;
    else if (s == "no-host")
      test_host = false;
    else if (s == "block")
      test_block = true;
    else if (s == "no-block")
      test_block = false;
    else if (s == "flat")
      test_flat = true;
    else if (s == "no-flat")
      test_flat = false;
    else if (s == "orig")
      test_orig = true;
    else if (s == "no-orig")
      test_orig = false;
    else if (s == "deg")
      test_deg = true;
    else if (s == "no-deg")
      test_deg = false;
    else if (s == "linear")
      test_lin = true;
    else if (s == "no-linear")
      test_lin = false;
    else if (s == "symmetric")
      symmetric = true;
    else if (s == "no-symmetric")
      symmetric = false;
    else if (s == "mkl")
      mkl = true;
    else if (s == "no-mkl")
      mkl = false;
    else if (s == "single")
      single = true;
    else if (s == "double")
      single = false;
    else if (s == "device") {
      ++i;
      device = std::atoi(argv[i]);
    }
    else if (s == "-h" || s == "--help")
      print_usage = true;
    else {
      std::cout << "Invalid argument:  " << s << std::endl;
      print_usage = true;
    }
    ++i;
  }
  if (print_usage) {
    std::cout << "Usage:" << std::endl
              << "\t" << argv[0]
              << " [no-][cuda|host|block|flat|orig|deg|linear|symmetric] [single|double] [device device_id]"
              << std::endl << "Defaults are all enabled." << std::endl;
    return -1;
  }

  // Always initialize threads (for Cuda::device_mirror)
  const size_t team_count =
    Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  const size_t threads_per_team =
    Kokkos::hwloc::get_available_threads_per_core();

#if defined(HAVE_STOKHOS_OPENMP) && defined(HAVE_STOKHOS_MKL)
  // Call a little OpenMP parallel region so that MKL will get the right
  // number of threads.  This isn't perfect in that the thread binding
  // doesn't seem right, and only works at all when using GNU threads with MKL.
#pragma omp parallel
  {
    int numThreads = omp_get_num_threads();
#pragma omp single
    std::cout << " num_omp_threads = " << numThreads << std::endl;
  }
#endif

  Kokkos::Threads::initialize( team_count , threads_per_team );

  if (test_host) {
    if (single)
      mainHost<float>(test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
    else
      mainHost<double>(test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
  }
  if (test_cuda) {
    if (single)
      mainCuda<float>(test_flat, test_orig, test_lin, test_block, symmetric, device);
    else
      mainCuda<double>(test_flat, test_orig, test_lin, test_block, symmetric, device);
  }

  Kokkos::Threads::finalize();

  return 0 ;
}
