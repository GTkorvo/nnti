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

#include "Kokkos_ThrustGPUNode.hpp"

#ifdef HAVE_KOKKOS_THRUST

#include "NodeTest.hpp"
#include <thrust/device_vector.h>

void thrust_float_alloc(int N, thrust::device_vector<float> &buff);
void thrust_int_alloc(int N, thrust::device_vector<int> &buff);
void thrust_float_init(thrust::device_vector<float> &buff);
float thrust_float_sum(const thrust::device_vector<float> &buff);
void thrust_int_init(thrust::device_vector<int> &buff);
int thrust_int_sum(const thrust::device_vector<int> &buff);

namespace {

  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> thrustNode_;

  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    return thrustNode_;
  }

  template <>
  void initNode<ThrustGPUNode>() {
    Teuchos::ParameterList plist;
    plist.set<int>("Device Number",NodeTest::cuda_dev);
    plist.set<int>("Verbose",NodeTest::verbose);
    thrustNode_ = rcp(new ThrustGPUNode(plist));
  }

  template <>
  std::pair<double,double> nativeTimings<float,ThrustGPUNode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,ThrustGPUNode init"), sTime("float,ThrustGPUNode sum");
    thrust::device_vector<float> buff;
    thrust_float_alloc(N,buff);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        thrust_float_init(buff);
      }
      thrustNode_->sync();
    }
    float sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        sum = thrust_float_sum(buff);
      }
      thrustNode_->sync();
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,ThrustGPUNode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,ThrustGPUNode init"), sTime("int,ThrustGPUNode sum");
    thrust::device_vector<int> buff;
    thrust_int_alloc(N,buff);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        thrust_int_init(buff);
      }
      thrustNode_->sync();
    }
    int sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        sum = thrust_int_sum(buff);
      }
      thrustNode_->sync();
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(ThrustGPUNode)
}

#endif
