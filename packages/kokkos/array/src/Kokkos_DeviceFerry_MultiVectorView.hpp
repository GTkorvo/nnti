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

#ifndef KOKKOS_DEVICEFERRY_MULTIVECTORVIEW_HPP
#define KOKKOS_DEVICEFERRY_MULTIVECTORVIEW_HPP

#include <Kokkos_MultiVectorView.hpp>

#include <Kokkos_DeviceFerry_macros.hpp>
#include <impl/Kokkos_MultiVectorView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceFerry , DeviceHost ,
                           false /* different memory spaces */ ,
                           true  /* both are contiguous */ >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceFerry > & dst ,
                   const MultiVectorView< ValueType , DeviceHost > & src )
  {
    const ValueType * s = src.m_memory.ptr_on_device();
    long int d = (long int)dst.m_memory.ptr_on_device();
    int size = src.length() * src.count();
    #pragma offload target(mic) in(s:length(size)) in(size) in(d)
    {
    	ValueType * local_d = (ValueType*)d;
    	for(int i = 0; i < size ; i++) {
    		local_d[i] = s[i];	
    	}
    }
    
/*    Impl::copy_to_ferry_from_host( dst.m_ptr_on_device ,
                                  src.m_ptr_on_device,
                                  sizeof(ValueType),
                                  dst.length() * dst.count() ); */
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceHost , DeviceFerry ,
                           false /* different memory spaces */ ,
                           true  /* both are contiguous */ >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceHost > & dst ,
                   const MultiVectorView< ValueType , DeviceFerry > & src )
  {
    long int s = (long int)src.m_memory.ptr_on_device();
    ValueType * d = dst.m_memory.ptr_on_device();
    int size = src.length() * src.count();
    #pragma offload target(mic) in(s) in(size) in(d:length(size))
    {
    	ValueType * local_s = (ValueType*)s;
    	for(int i = 0; i < size ; i++) {
    		d[i] = local_s[i];	
    	}
    }
/*    Impl::copy_to_host_from_ferry( dst.m_ptr_on_device ,
                                  src.m_ptr_on_device,
                                  sizeof(ValueType),
                                  dst.length() * dst.count() ); */
  }
};

} // namespace Impl
} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICEFERRY_MULTIVECTORVIEW_HPP */


