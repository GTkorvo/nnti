/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <iostream>

#include <impl/Kokkos_MemoryView.hpp>
#include <Kokkos_Value.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_MDArray.hpp>
#include <Kokkos_CrsMap.hpp>

#include <Kokkos_Host.hpp>
#include <Kokkos_Cuda.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_Cuda_macros.hpp>

#include <TestMemoryManagement.hpp>
#include <TestValue.hpp>
#include <TestMultiVector.hpp>
#include <TestMDArray.hpp>
#include <TestCrsMap.hpp>

#include <impl/Kokkos_MDArrayIndexMapRight_macros.hpp>
#include <TestMDArrayIndexMap.hpp>

#include <TestReduce.hpp>
#include <TestMultiReduce.hpp>

#include <Kokkos_Clear_macros.hpp>

namespace Test {

void test_device_cuda_init() {
  Kokkos::Cuda::initialize();
}

void test_device_cuda_memory_management() {
  TestMemoryManagement< int, Kokkos::Cuda >();
  TestMemoryManagement< double, Kokkos::Cuda >();
}

void test_device_cuda_value() {
  TestValue< int, Kokkos::Cuda >();
  TestValue< double, Kokkos::Cuda >();
}

void test_device_cuda_multi_vector() {
  TestMultiVector< int, Kokkos::Cuda >();
  TestMultiVector< double, Kokkos::Cuda >();
}

void test_device_cuda_crsmap() {
  TestCrsMap< Kokkos::Cuda >();
}

void test_device_cuda_mdarray() {
  TestMDArray< double, Kokkos::Cuda >();
  TestMDArray< int, Kokkos::Cuda >();
}

void test_device_cuda_index_map() {
  TestMDArrayIndexMap< Kokkos::Cuda >();
}

void test_device_cuda_reduce() {
  TestReduce< long ,   Kokkos::Cuda >( 1000000 );
  TestReduce< double ,   Kokkos::Cuda >( 1000000 );
}

void test_device_cuda_multi_reduce() {
  TestReduceMulti< long , Kokkos::Cuda >( 1000000 , 7 );
}

}

