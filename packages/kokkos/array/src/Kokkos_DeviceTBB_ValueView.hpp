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

#ifndef KOKKOS_DEVICETBB_MDARRAYDEEPCOPY_HPP
#define KOKKOS_DEVICETBB_MDARRAYDEEPCOPY_HPP

#include <Kokkos_ValueView.hpp>

#include <Kokkos_DeviceTBB_macros.hpp>
#include <impl/Kokkos_ValueView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceTBB , DeviceHost > {
public:

  static void run( const ValueView< ValueType , DeviceTBB > & dst ,
                   const ValueType & src )
  { *dst = src ; }

  static void run( const ValueView< ValueType , DeviceTBB >  & dst ,
                   const ValueView< ValueType , DeviceHost > & src )
  { *dst = *src ; }
};

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceHost , DeviceTBB > {
public:

  static void run( ValueType & dst ,
                   const ValueView< ValueType , DeviceTBB >  & src )
  { dst = *src ; }

  static void run( const ValueView< ValueType , DeviceHost > & dst ,
                   const ValueView< ValueType , DeviceTBB >  & src )
  { *dst = *src ; }
};

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_DEVICETBB_MDARRAYDEEPCOPY_HPP */


