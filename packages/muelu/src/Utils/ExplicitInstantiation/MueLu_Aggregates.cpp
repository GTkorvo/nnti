// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_Aggregates_def.hpp"

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
template class MueLu::Aggregates<int, int, KokkosClassic::DefaultNode::DefaultNodeType, KokkosClassic::DefaultKernels<void, int, KokkosClassic::DefaultNode::DefaultNodeType>::SparseOps>;
#endif

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::Aggregates<int, long long int, KokkosClassic::DefaultNode::DefaultNodeType, KokkosClassic::DefaultKernels<void, int, KokkosClassic::DefaultNode::DefaultNodeType>::SparseOps>;
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

#if defined(HAVE_KOKKOSCLASSIC_THRUST) && defined(HAVE_KOKKOSCLASSIC_CUDA_DOUBLE) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
template class MueLu::Aggregates<int, int, KokkosClassic::ThrustGPUNode, KokkosClassic::DefaultKernels<void, int, KokkosClassic::ThrustGPUNode>::SparseOps>;
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
template class MueLu::Aggregates<int, int, KokkosClassic::TPINode, KokkosClassic::DefaultKernels<void, int, KokkosClassic::TPINode>::SparseOps>;
#endif


