/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

// To force use of OMP atomics instead of intrinsics
// #define KOKKOS_ATOMICS_USE_OMP31

#include <Kokkos_Core.hpp>

#include <Kokkos_Bitset.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_Vector.hpp>

//----------------------------------------------------------------------------
#include <TestBitset.hpp>
#include <TestUnorderedMap.hpp>
#include <TestVector.hpp>
#include <TestDualView.hpp>
#include <TestSegmentedView.hpp>


#include <iomanip>

namespace Test {

#ifdef KOKKOS_HAVE_OPENMP
class openmp : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;

    unsigned threads_count = 4 ;

    if ( Kokkos::hwloc::available() ) {
      threads_count = Kokkos::hwloc::get_available_numa_count() *
                      Kokkos::hwloc::get_available_cores_per_numa();
    }

    Kokkos::OpenMP::initialize( threads_count );
  }

  static void TearDownTestCase()
  {
    Kokkos::OpenMP::finalize();
  }
};

TEST_F( openmp, bitset )
{
  test_bitset<Kokkos::OpenMP>();
}

#define OPENMP_INSERT_TEST( name, num_nodes, num_inserts, num_duplicates, repeat, near )                                \
  TEST_F( openmp, UnorderedMap_insert_##name##_##num_nodes##_##num_inserts##_##num_duplicates##_##repeat##x) {   \
    for (int i=0; i<repeat; ++i)                                                                                \
      test_insert<Kokkos::OpenMP>(num_nodes,num_inserts,num_duplicates, near);                                   \
  }

#define OPENMP_FAILED_INSERT_TEST( num_nodes, repeat )                         \
  TEST_F( openmp, UnorderedMap_failed_insert_##num_nodes##_##repeat##x) {     \
    for (int i=0; i<repeat; ++i)                                               \
      test_failed_insert<Kokkos::OpenMP>(num_nodes);                             \
  }

#define OPENMP_ASSIGNEMENT_TEST( num_nodes, repeat )                             \
  TEST_F( openmp, UnorderedMap_assignment_operators_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_assignement_operators<Kokkos::OpenMP>(num_nodes);                     \
  }

#define OPENMP_DEEP_COPY( num_nodes, repeat )                             \
  TEST_F( openmp, UnorderedMap_deep_copy##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_deep_copy<Kokkos::OpenMP>(num_nodes);                     \
  }

#define OPENMP_VECTOR_COMBINE_TEST( size )                             \
  TEST_F( openmp, vector_combination##size##x) {       \
      test_vector_combinations<int,Kokkos::OpenMP>(size);                     \
  }

#define OPENMP_DUALVIEW_COMBINE_TEST( size )                             \
  TEST_F( openmp, dualview_combination##size##x) {       \
      test_dualview_combinations<int,Kokkos::OpenMP>(size);                     \
  }

#define OPENMP_SEGMENTEDVIEW_TEST( size )                             \
  TEST_F( openmp, segmentedview_##size##x) {       \
      test_segmented_view<double,Kokkos::OpenMP>(size);                     \
  }

OPENMP_INSERT_TEST(close, 100000, 90000, 100, 500, true)
OPENMP_INSERT_TEST(far, 100000, 90000, 100, 500, false)
OPENMP_FAILED_INSERT_TEST( 10000, 1000 )
OPENMP_DEEP_COPY( 10000, 1 )

OPENMP_VECTOR_COMBINE_TEST( 10 )
OPENMP_VECTOR_COMBINE_TEST( 3057 )
OPENMP_DUALVIEW_COMBINE_TEST( 10 )
OPENMP_SEGMENTEDVIEW_TEST( 10000 )

#undef OPENMP_INSERT_TEST
#undef OPENMP_FAILED_INSERT_TEST
#undef OPENMP_ASSIGNEMENT_TEST
#undef OPENMP_DEEP_COPY
#undef OPENMP_VECTOR_COMBINE_TEST
#undef OPENMP_DUALVIEW_COMBINE_TEST
#undef OPENMP_SEGMENTEDVIEW_TEST
#endif
} // namespace test

