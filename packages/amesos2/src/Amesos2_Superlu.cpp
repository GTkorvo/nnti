// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

#include "Amesos2_config.h"
#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#  include "Amesos2_Superlu_decl.hpp"
#  include "Amesos2_Superlu_def.hpp"
#  include "Amesos2_ExplicitInstantiationHelpers.hpp"

namespace Amesos2 {
#ifdef HAVE_AMESOS2_EPETRA
  AMESOS2_SOLVER_EPETRA_INST(Superlu);
#endif

#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superlu,float,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superlu,float,int,unsigned int);
#endif
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superlu,double,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superlu,double,int,unsigned int);
#endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<float>,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<float>,int,unsigned int);
#endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<double>,int,int);
#ifdef HAVE_TPETRA_INST_INT_UNSIGNED
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<double>,int,unsigned int);
#endif
#endif

}

#ifdef HAVE_TPETRA_INST_INT_LONG
namespace Amesos2 {
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superlu,float,int,long);
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superlu,double,int,long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<float>,int,long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<double>,int,long);
#endif
}
#endif

#ifdef HAVE_TPETRA_INST_INT_LONG_LONG
namespace Amesos2 {
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superlu,float,int,long long);
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superlu,double,int,long long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<float>,int,long long);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(Superlu,std::complex<double>,int,long long);
#endif
}
#endif
//
// mfh 16 Jan 2014: Hack to make explicit instantiation for
// Ifpack2::Details::Amesos2Wrapper work.
//
#include "Kokkos_DefaultNode.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#define AMESOS2_SUPERLU_LOCAL_INSTANT(S,LO,GO,N)                        \
  template class Amesos2::Superlu<Tpetra::CrsMatrix<S, LO, GO, N>,      \
                                  Tpetra::MultiVector<S, LO, GO,  N> >;

TPETRA_ETI_MANGLING_TYPEDEFS()

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && defined(HAVE_TPETRA_INST_DOUBLE) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)
  AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, int, KokkosClassic_TPINode)
#endif

#if defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_TPETRA_INST_DOUBLE) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THRUSTGPUNODE)
  AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, int, KokkosClassic_ThrustGPUNode)
#endif

#if defined(HAVE_TPETRA_INST_SERIAL) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#define NODETYPE Kokkos_Compat_KokkosSerialWrapperNode
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#if defined(HAVE_TPETRA_INST_PTHREAD) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#define NODETYPE Kokkos_Compat_KokkosThreadsWrapperNode
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#if defined(HAVE_TPETRA_INST_OPENMP) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#define NODETYPE Kokkos_Compat_KokkosOpenMPWrapperNode
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#if defined(HAVE_TPETRA_INST_CUDA) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#define NODETYPE Kokkos_Compat_KokkosCudaWrapperNode
#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(float, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(double, int, unsigned int, NODETYPE)
    #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, int, NODETYPE)
  #ifdef HAVE_TPETRA_INST_INT_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, long long, NODETYPE)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<float>, int, unsigned int, NODETYPE)
  #endif
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, int, NODETYPE)
    #ifdef HAVE_TPETRA_INST_INT_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, long long, NODETYPE)
    #endif
    #ifdef HAVE_TPETRA_INST_INT_UNSIGNED
      AMESOS2_SUPERLU_LOCAL_INSTANT(std::complex<double>, int, unsigned int, NODETYPE)
    #endif
#endif
#undef NODETYPE
#endif

#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
