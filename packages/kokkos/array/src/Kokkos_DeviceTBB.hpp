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

#ifndef KOKKOS_DEVICETBB_HPP
#define KOKKOS_DEVICETBB_HPP

#include <Kokkos_DeviceHost.hpp>

#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>

#define KOKKOS_DEVICE_TBB  Kokkos::DeviceTBB

/*--------------------------------------------------------------------------*/

namespace Kokkos {
	
class DeviceTBB {
public:

  /** \brief  On the TBB device use size_t for indexing */
  typedef size_t               size_type ;

  /** \brief  The TBB device uses the Host memory space */
  typedef HostMemory           memory_space ;

  /** \brief  Default mdarray map is index from right */
  typedef Impl::MDArrayIndexMapRight  mdarray_map ;

  /*--------------------------------*/

  static void initialize( size_type nthreads );

  static void finalize();

  /*--------------------------------*/

  static void wait_functor_completion() {}

  /*--------------------------------*/
};

} // namespace Kokkos

#include <Kokkos_DeviceTBB_macros.hpp>
#include <impl/Kokkos_BasicFunctors_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

#endif /* #define KOKKOS_DEVICETBB_HPP */

