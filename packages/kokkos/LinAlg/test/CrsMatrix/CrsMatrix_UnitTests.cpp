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

#include <Teuchos_UnitTestHarness.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::arcp;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  int N = 100;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  // intialize using static graph
  TEUCHOS_UNIT_TEST( CrsGraph, RuntimeExceptions )
  {
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps SparseOps;
    typedef SparseOps::graph<int,Node>::graph_type             Graph;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    RCP<Graph> G = rcp(new Graph(N,node));
    {
      ArrayRCP<size_t>  ptrs_tooSmall(N), ptrs_tooBig(N+2);
      ArrayRCP<int>     inds;
      TEST_THROW( G->setStructure(ptrs_tooSmall, inds), std::runtime_error );
      TEST_THROW( G->setStructure(ptrs_tooBig,   inds), std::runtime_error );
    }
    {
      ArrayRCP<size_t> ptrs(N+1);
      for (size_t i=0; i<=N; ++i) ptrs[i] = i;
      ArrayRCP<int> tooFewInds(N-1);
      TEST_THROW( G->setStructure(ptrs, tooFewInds), std::runtime_error );
    }
  }

  TEUCHOS_UNIT_TEST( CrsGraph, StatusTests )
  {
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps SparseOps;
    typedef SparseOps::graph<int,Node>::graph_type             Graph;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    RCP<Graph> G = rcp(new Graph(N,node));
    G->finalize(null);
    TEST_EQUALITY_CONST( G->isEmpty(), true );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, StaticGraph, Scalar, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<Scalar,Ordinal,Node>::SparseOps SparseOps;
    typedef SparseOps::graph<Ordinal,Node>::graph_type                      Graph;
    typedef SparseOps::matrix<Scalar,int,Node>::matrix_type                 Matrix;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    // build a non-empty graph, tridiagonal
    const size_t testNumEntries = 3*N-2;
    ArrayRCP<Ordinal> inds(testNumEntries);
    ArrayRCP<Scalar > vals(testNumEntries);
    ArrayRCP<size_t> ptrs(N+1);
    {
      std::fill( inds.begin(), inds.end(), 0 );
      std::fill( vals.begin(), vals.end(), 0 );
      size_t curoffset = 0;
      for (size_t r=0; r < N; ++r) {
        ptrs[r] = curoffset;
        if (r > 0 && r < N-1) curoffset += 3;
        else                  curoffset += 2;
      }
      ptrs[N] = curoffset;
    }
    RCP<Graph> G = rcp(new Graph(N,node));
    TEST_EQUALITY( G->getNumRows(), N );
    TEST_EQUALITY( G->getNode(), node );
    G->setStructure(ptrs, inds);
    ArrayRCP<Ordinal> chkInds;
    ArrayRCP<size_t> chkPtrs;
    chkInds = G->getIndices();
    chkPtrs = G->getPointers();
    TEST_EQUALITY( inds, chkInds );
    TEST_EQUALITY( ptrs, chkPtrs );
    TEST_EQUALITY_CONST( G->isFinalized(), false );
    G->finalize(null);
    TEST_EQUALITY_CONST( G->isFinalized(), true );
    TEST_EQUALITY_CONST( G->isEmpty(), false );
    TEST_EQUALITY( staticG->getNumRows(), N );
    TEST_EQUALITY( staticG->getNode(), node );
    Matrix M(G);
    TEST_EQUALITY( M.getNumRows(),    G->getNumRows() );
    M.setValues(vals);
    TEST_EQUALITY_CONST( M.isEmpty(), false );        // can only query this pre-finalize for a static graph
    TEST_EQUALITY_CONST( M.isFinalized(), false );
    M.finalize(null);
      // HERE
      }
      TEST_EQUALITY_CONST( M.isEmpty(), false );
      TEST_EQUALITY_CONST( M.isFinalized(), true );
      TEST_EQUALITY_CONST( M.isOptimized(), G->isOptimized() );
      {
        // test graph clear
        G->clear();
        TEST_EQUALITY_CONST( G->isFinalized(), false );
        TEST_EQUALITY( G->getNumRows(), N );
        TEST_EQUALITY( G->getNumEntries(), 0 );
        TEST_EQUALITY_CONST( G->is1DStructure(), false );
        TEST_EQUALITY_CONST( G->is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G->get1DStructure(chkInds1D, chkBegs, chkEnds);
        G->get2DStructure(chkInds2D, chkNumPer);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
    }
  }

  // intialize using static graph
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, DynamicGraph, Scalar, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<Scalar,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                   Graph;
    typedef CrsMatrix<Scalar,Ordinal,Node,SparseOps>                           Matrix;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a non-empty graph, tridiagonal
    const size_t testNumEntries = 3*N-2;
    ArrayRCP<Ordinal> inds1D(testNumEntries);
    ArrayRCP<Scalar > vals1D(testNumEntries);
    ArrayRCP<ArrayRCP<Ordinal> > inds2D(N);
    ArrayRCP<ArrayRCP<Scalar > > vals2D(N);
    ArrayRCP<size_t> begs(N+1), ends(N), numper(N);
    {
      std::fill( inds1D.begin(), inds1D.end(), 0 );
      std::fill( vals1D.begin(), vals1D.end(), 0 );
      size_t curoffset = 0;
      for (size_t r=0; r < N; ++r) {
        if (r > 0 && r < N-1) numper[r] = 3;
        else numper[r] = 2;
        begs[r] = curoffset;
        ends[r] = begs[r] + numper[r];
        inds2D[r] = inds1D.persistingView(begs[r],numper[r]);
        vals2D[r] = vals1D.persistingView(begs[r],numper[r]);
        curoffset += numper[r];
      }
      begs[N] = curoffset;
    }
    Graph G(N,node);
    Matrix M(G);
    for (int t=0; t < 4; ++t) {
      const bool submit1D        = (t && 1);
      const bool OptimizeStorage = (t && 2);
      if (submit1D) {
        G.set1DStructure(inds1D, begs, ends);
        M.set1DValues(vals1D);
        TEST_EQUALITY_CONST( G.is1DStructure(), true );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        TEST_EQUALITY_CONST( M.is1DStructure(), true );
        TEST_EQUALITY_CONST( M.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds;
        ArrayRCP<size_t> chkBegs, chkEnds;
        ArrayRCP<Scalar> chkVals;
        G.get1DStructure(chkInds, chkBegs, chkEnds);
        M.get1DValues(chkVals);
        TEST_EQUALITY( inds1D, chkInds );
        TEST_EQUALITY( begs, chkBegs );
        TEST_EQUALITY( ends, chkEnds );
        TEST_EQUALITY( vals1D, chkVals );
      }
      else {
        G.set2DStructure(inds2D, numper);
        M.set2DValues(vals2D);
        TEST_EQUALITY_CONST( G.is2DStructure(), true );
        TEST_EQUALITY_CONST( G.is1DStructure(), false );
        TEST_EQUALITY_CONST( M.is2DStructure(), true );
        TEST_EQUALITY_CONST( M.is1DStructure(), false );
        ArrayRCP<ArrayRCP<Ordinal> > chkInds;
        ArrayRCP<ArrayRCP<Scalar> > chkVals;
        ArrayRCP<size_t> chkNumPer;
        G.get2DStructure(chkInds, chkNumPer);
        M.get2DValues(chkVals);
        TEST_EQUALITY( inds2D, chkInds );
        TEST_EQUALITY( numper, chkNumPer );
        TEST_EQUALITY( vals2D, chkVals );
      }
      TEST_EQUALITY_CONST( G.isFinalized(), false );
      TEST_EQUALITY_CONST( G.isOptimized(), false );
      TEST_EQUALITY_CONST( G.getNumEntries(), testNumEntries );
      TEST_EQUALITY_CONST( M.isFinalized(), false );
      TEST_EQUALITY_CONST( M.isOptimized(), false );
      M.finalize(OptimizeStorage);
      TEST_EQUALITY_CONST( G.isFinalized(), true );
      TEST_EQUALITY_CONST( G.isEmpty(), false );
      TEST_EQUALITY_CONST( M.isFinalized(), true );
      TEST_EQUALITY_CONST( M.isEmpty(), false );
      if (OptimizeStorage) {
        TEST_EQUALITY_CONST( G.isOptimized(), true );
        TEST_EQUALITY_CONST( G.is1DStructure(), true );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        TEST_EQUALITY_CONST( M.isOptimized(), true );
        TEST_EQUALITY_CONST( M.is1DStructure(), true );
        TEST_EQUALITY_CONST( M.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<ArrayRCP<Scalar> >  chkVals2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        ArrayRCP<Scalar> chkVals1D;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        M.get1DValues(chkVals1D);
        M.get2DValues(chkVals2D);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkVals1D == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkVals2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
      else {
        TEST_EQUALITY_CONST( G.isOptimized(), false );
        TEST_EQUALITY( G.is1DStructure(),  submit1D );
        TEST_EQUALITY( G.is2DStructure(), !submit1D );
        TEST_EQUALITY_CONST( M.isOptimized(), false );
        TEST_EQUALITY( M.is1DStructure(),  submit1D );
        TEST_EQUALITY( M.is2DStructure(), !submit1D );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<Scalar>  chkVals1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<ArrayRCP<Scalar> >  chkVals2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        M.get1DValues(chkVals1D);
        M.get2DValues(chkVals2D);
        TEST_EQUALITY( chkInds1D == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkVals1D == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkBegs   == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkEnds   == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkInds2D == Teuchos::null,  submit1D );
        TEST_EQUALITY( chkVals2D == Teuchos::null,  submit1D );
        TEST_EQUALITY( chkNumPer == Teuchos::null,  submit1D );
      }
      {
        // test matrix and graph clear
        M.clear();
        TEST_EQUALITY( M.getNumRows(), N );
        TEST_EQUALITY( G.getNumRows(), N );
        TEST_EQUALITY_CONST( G.getNumEntries() != 0, true );
        TEST_EQUALITY_CONST( G.is1DStructure() || G.is2DStructure(), true );
        TEST_EQUALITY_CONST( M.isFinalized(), false );
        TEST_EQUALITY_CONST( G.isFinalized(), true );
        G.clear();
        TEST_EQUALITY( G.getNumRows(), N );
        TEST_EQUALITY_CONST( G.getNumEntries() == 0, true );
        TEST_EQUALITY_CONST( G.is1DStructure(), false );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        TEST_EQUALITY_CONST( G.isFinalized(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<Scalar>  chkVals1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<ArrayRCP<Scalar> >  chkVals2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        M.get1DValues(chkVals1D);
        M.get2DValues(chkVals2D);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkVals1D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkVals2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
    }
  }

  // bad ends
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraph, BadEnds, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<float,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                   Graph;
    const size_t N = 1;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a null-size graph
    Graph G(N,node);
    {
      ArrayRCP<size_t> begs, ends;
      ArrayRCP<Ordinal> inds1D;
      begs = arcp<size_t>(N+1);
      ends = arcp<size_t>(N);
      //
#ifdef HAVE_KOKKOS_DEBUG
      begs[0] = 0; begs[1] = 0;
      ends[0] = 1;  // ends causes rows to overlap; not consistent
      TEST_THROW( G.set1DStructure(inds1D, begs, ends), std::runtime_error );
#endif
      //
      begs[0] = 0; begs[1] = 2; // begs size is larger than allocation in inds1D
      ends[0] = 0;
      inds1D = arcp<Ordinal>(1);
      TEST_THROW( G.set1DStructure(inds1D, begs, ends), std::runtime_error );
      // 
      begs[0] = 0; begs[1] = 1;
      ends[0] = 0;
      // this allocation is allowed to be too large, w.r.t. begs
      inds1D = arcp<Ordinal>(2);
      G.set1DStructure(inds1D, begs, ends);
    }
  }

  // no rows
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraph, NoRows, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<float,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                   Graph;
    const size_t N = 0;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a null-size graph
    Graph G(N,node);
    {
      ArrayRCP<size_t> begs, ends;
      ArrayRCP<Ordinal> inds1D;
      begs = null;
      ends = null;
      inds1D = null;
      // begs is null; not allowed
      TEST_THROW( G.set1DStructure(inds1D, begs, ends), std::runtime_error );
      begs = arcp<size_t>(1);
      begs[0] = 0;
      G.set1DStructure(inds1D, begs, ends);
      begs = null;
      ends = null;
      inds1D = null;
      G.get1DStructure(inds1D, begs, ends);
      TEST_EQUALITY_CONST( begs == null, false );
      TEST_EQUALITY_CONST( ends == null, true );
      TEST_EQUALITY_CONST( inds1D == null, true );
    }
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, StaticGraph,  SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, DynamicGraph, SCALAR, ORDINAL )

 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraph, NoRows,  ORDINAL ) \
          TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraph, BadEnds,  ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
