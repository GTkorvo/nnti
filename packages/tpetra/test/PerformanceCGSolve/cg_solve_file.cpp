/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "cg-solve_file.hpp"
#ifndef USE_KOKKOS_CLASSIC
  #ifdef KOKKOS_HAVE_CUDA
    int run_cudawrapper(int& argc, char* argv[]);
  #endif
  #ifdef KOKKOS_HAVE_OPENMP
    int run_openmpwrapper(int& argc, char* argv[]);
  #endif
  #ifdef KOKKOS_HAVE_PTHREAD
    int run_threadswrapper(int& argc, char* argv[]);
  #endif
#else
  int run_tpi(int argc, char* argv[]);
#endif

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  //
  // Specify types used in this example
  //

#ifndef USE_KOKKOS_CLASSIC
  #ifdef KOKKOS_HAVE_CUDA
    run_cudawrapper(argc,argv);
  #endif
  #ifdef KOKKOS_HAVE_OPENMP
    run_openmpwrapper(argc,argv);
  #endif
  #ifdef KOKKOS_HAVE_PTHREAD
    run_threadswrapper(argc,argv);
  #endif
#else
    run_tpi(argc,argv);
#endif
  return 0;
}
