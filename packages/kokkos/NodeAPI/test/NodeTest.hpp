#ifndef NODE_TEST_HPP_
#define NODE_TEST_HPP_ 

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>

#include <functional>
#include <algorithm>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "Kokkos_NodeHelpers.hpp"

#include "TestOps.hpp"

namespace NodeTest {
  extern int N;
  extern int numIters;
  extern int numThreads;
  extern int verbose;
  extern int cuda_dev;
}

namespace {

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::Time;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Kokkos::ReadyBufferHelper;
  using Teuchos::tuple;
  using Teuchos::ParameterList;

  template <class NODE>
  RCP<NODE> getNode() {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  template <class NODE>
  void initNode() {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  template <class SCALAR, class NODE>
  std::pair<double,double> nativeTimings(int N, int numIters, SCALAR &result) {
    // TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
    // FINISH
    return std::pair<double,double>(0.0,0.0);
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AAAAA_Is_First, InitNode, NODE )
  {
    out << "Initializing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    initNode<NODE>();
    TEST_EQUALITY(0,0);
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NodeAPI, MemoryInitTest, NODE )
  {
    out << "Testing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    ArrayRCP<int> x, y;
    RCP<NODE> node = getNode<NODE>();
    ReadyBufferHelper<NODE> rbh(node);
    const int N = 4096;
    // const int N = 1000000;
    const int Xoff  = 1024;
    const int Xsize = N + 2*Xoff;  // 6144
    const int Yoffs[5] = {0, N, 2*N, 3*N, 4*N};
    // allocate two compute buffers: 
    x = node->template allocBuffer<int>(Xsize);
    y = node->template allocBuffer<int>(Yoffs[4]);
    // set x = {?  fours  ?} using copyTo
    {
      Array<int> fours(N);
      std::fill(fours.begin(), fours.end(), 4);
      node->template copyToBuffer<int>(N,fours(),x+Xoff);
    }
    // set y = {? ? ? fours} using copyBuffers
    node->template copyBuffers<int>(N,x+Xoff,y+Yoffs[3]);
    // set y = {? ? threes fours} using local viewBuffer(WriteOnly)
    {
      ArrayRCP<int> view = node->viewBufferNonConst(Kokkos::WriteOnly, N, y+Yoffs[2]);
      TEST_EQUALITY_CONST( view.size(), N );
      std::fill(view.begin(), view.end(), 3);
    }
    // set y = {ones ones threes fours} using init
    {
      InitOp<int> iop;
      rbh.begin();
      iop.x = rbh.template addNonConstBuffer<int>(y);
      rbh.end();
      const int beg = Yoffs[0],
                end = Yoffs[2];
      node->parallel_for(beg,end,iop);
    }
    // set y = {ones twos threes fours} using viewBuffer(ReadWrite) and local arithmetic
    // over-size the view, to ensure that we don't affect any more changes than necessary
    {
      ArrayRCP<int> view = node->viewBufferNonConst(Kokkos::ReadWrite, 2*N, y+Yoffs[1]);
      // add 1 to all elements of view
      std::transform(view.begin(), view.begin()+N, view.begin(), std::bind2nd(std::plus<int>(), 1));
    }
    // check the result
    { 
      Array<int> expected(N);
      ArrayRCP<const int> yview = node->template viewBuffer<int>(4*N,y);
      // y should be {ones  twos  threes  fours}
      // check ones
      std::fill(expected.begin(), expected.end(), 1);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[0],N) ); 
      // check twos
      std::fill(expected.begin(), expected.end(), 2);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[1],N) ); 
      // check threes
      std::fill(expected.begin(), expected.end(), 3);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[2],N) ); 
      // check fours
      std::fill(expected.begin(), expected.end(), 4);
      TEST_COMPARE_ARRAYS( expected(), yview(Yoffs[3],N) ); 
    }
    // free the allocations
    x = Teuchos::null;
    y = Teuchos::null;
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( NodeAPI, TestAndTime, SCALAR, NODE )
  {
    using std::endl;
    out << "Testing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    Time tAlloc("Alloc Time"), tInit("Init Op"), tSum("Sum Op"), tFree("Free Time");
    Teuchos::ArrayRCP<SCALAR> x;
    RCP<NODE> node = getNode<NODE>();
    ReadyBufferHelper<NODE> rbh(node);
    SCALAR result;
    int N = NodeTest::N;
    int numIters = NodeTest::numIters;
    {
      TimeMonitor localTimer(tAlloc);
      x = node->template allocBuffer<SCALAR>(N);
    }
    ////////////////////////////////////////////////////////
    // timing loop around full vector init (parallel_for)
    {
      TimeMonitor localTimer(tInit);
      for (int t=0; t < numIters; ++t) {
        // set x[i] = 1, i=0:N-1
        InitOp<SCALAR> wdp;
        rbh.begin();
        wdp.x = rbh.template addNonConstBuffer<SCALAR>(x);
        rbh.end();
        node->parallel_for(0,N,wdp);
      }
    }
    // 
    // test partial sum of x[i], i=1:N-2
    SCALAR expectedResult;
    {
      SumOp<SCALAR> wdp;
      rbh.begin();
      wdp.x = rbh.template addConstBuffer<SCALAR>(x);
      rbh.end();
      result = node->parallel_reduce(1,N-1,wdp);
      expectedResult = (SCALAR)(N-2);
      TEST_EQUALITY(result, expectedResult);
    }
    ////////////////////////////////////////////////////////
    // timing loop around full vector sum (parallel_reduce)
    {
      TimeMonitor localTimer(tSum);
      for (int t=0; t < numIters; ++t) {
        // compute sum x[i], i=0:N-1
        SumOp<SCALAR> wdp;
        rbh.begin();
        wdp.x = rbh.template addConstBuffer<SCALAR>(x);
        rbh.end();
        result = node->parallel_reduce(0,N,wdp);
      }
    }
    expectedResult = (SCALAR)(N);
    TEST_EQUALITY(result, expectedResult);
    // free
    {
      TimeMonitor localTimer(tFree);
      x = Teuchos::null;
    }
    ////////////////////////////////////////////////////////
    // get native times
    SCALAR nativeResult;
    std::pair<double,double> ntimes = nativeTimings<SCALAR,NODE>(N, numIters, nativeResult);
    TEST_EQUALITY(nativeResult, expectedResult);
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    out << "Timing results for " << Teuchos::TypeNameTraits<NODE>::name() << endl
        << "  allocBuffer   time: " << std::scientific << std::setprecision(2) << tAlloc.totalElapsedTime()          << endl
        << "  freeBuffer    time: " << std::scientific << std::setprecision(2) << tFree.totalElapsedTime()           << endl
        << "  Kokkos InitOp time: " << std::scientific << std::setprecision(2) << tInit.totalElapsedTime()/numIters  << endl
        << "  Native init   time: " << std::scientific << std::setprecision(2) <<             ntimes.first/numIters  << endl
        << "  Kokkos SumOp  time: " << std::scientific << std::setprecision(2) <<  tSum.totalElapsedTime()/numIters  << endl
        << "  Native sum    time: " << std::scientific << std::setprecision(2) <<            ntimes.second/numIters  << endl;
  }

  // 
  // INSTANTIATIONS
  //

  #define TEST_NODE(NODE) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( AAAAA_Is_First, InitNode,    NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NodeAPI, MemoryInitTest,     NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, TestAndTime, int,   NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( NodeAPI, TestAndTime, float, NODE )

}

#endif
