/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#include "Kokkos_DummySparseKernelClass.hpp"
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>

/** \file DummySparseKernelDriver.cpp
    \brief A file testing the build of the DummySparseKernel class and illustrating its usage.
 */

int main() {

  typedef Kokkos::DefaultNode::DefaultNodeType                          Node;
  typedef KokkosExamples::DummySparseKernel<Node>                  SparseOps;
  typedef typename SparseOps::graph<int,Node>::graph_type              Graph;
  typedef typename SparseOps::matrix<double,int,Node>::matrix_type   DMatrix;
  typedef typename SparseOps::matrix< float,int,Node>::matrix_type   FMatrix;
  typedef Kokkos::MultiVector<double,Node>                         DoubleVec;
  typedef Kokkos::MultiVector<float,Node>                           FloatVec;

  std::cout << "Note, this class doesn't actually do anything. We are only testing that it compiles." << std::endl;

  // get a pointer to the default node
  Teuchos::RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

  // create the graph G
  const size_t numRows = 5;
  Teuchos::RCP<Graph> G = Teuchos::rcp(new Graph(numRows,node));

  // create a double-valued matrix dM using the graph G
  Teuchos::RCP<DMatrix> dM = Teuchos::rcp(new DMatrix(G.getConst()));
  // create a double-valued sparse kernel using the rebind functionality
  SparseOps::rebind<double>::other doubleKernel(node);
  // initialize it with G and dM
  doubleKernel.setGraphAndMatrix(G,dM);
  // create double-valued vectors and initialize them
  DoubleVec dx(node), dy(node);
  // test the sparse kernel operator interfaces
  doubleKernel.multiply( Teuchos::NO_TRANS, 1.0, dx, dy);
  doubleKernel.multiply( Teuchos::NO_TRANS, 1.0, dx, 1.0, dy);
  doubleKernel.solve( Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::UNIT_DIAG, dy, dx);

  // create a float-valued matrix fM using the graph G
  Teuchos::RCP<FMatrix> fM = Teuchos::rcp(new FMatrix(G));
  // create a double-valued sparse kernel using the rebind functionality
  SparseOps::rebind<float>::other floatKernel(node);
  // initialize it with G and fM
  floatKernel.setGraphAndMatrix(G,fM);
  // create float-valued vectors and initialize them
  FloatVec fx(node), fy(node);
  // test the sparse kernel operator interfaces
  floatKernel.multiply( Teuchos::NO_TRANS, 1.0f, fx, fy);
  floatKernel.multiply( Teuchos::NO_TRANS, 1.0f, fx, 1.0f, fy);
  floatKernel.solve( Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::UNIT_DIAG, fy, fx);

  std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
