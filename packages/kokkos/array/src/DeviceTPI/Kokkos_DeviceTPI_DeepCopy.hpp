/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_DEVICE_TPI_DEEP_COPY_HPP
#define KOKKOS_DEVICE_TPI_DEEP_COPY_HPP

#include <Kokkos_ArrayForwardDeclarations.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {

/*------------------------------------------------------------------------*/

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceTPI >
{
public:
  static void run( const ValueView< ValueType , DeviceTPI > & dst ,
                   const ValueType & src )
  { *dst = src ; }

  static void run( ValueType & dst ,
                   const ValueView< ValueType , DeviceTPI > & src )
  { dst = *src ; }
};

/*------------------------------------------------------------------------*/

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceTPI , DeviceHost >
{
public:
  typedef Impl::CopyFunctor< ValueType , DeviceTPI > functor_type ;

  static void run( const MultiVectorView< ValueType , DeviceTPI > & dst ,
                   const MultiVectorView< ValueType , DeviceHost > & src )
  {
    Impl::multivector_require_equal_dimension( dst , src );

    parallel_for( dst.size() ,
                  functor_type( dst.m_ptr_on_device ,
                                src.m_ptr_on_device ) );
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceHost , DeviceTPI >
{
  typedef Impl::CopyFunctor< ValueType , DeviceTPI > functor_type ;

  static void run( const MultiVectorView< ValueType , DeviceHost > & dst ,
                   const MultiVectorView< ValueType , DeviceTPI > & src )
  {
    Impl::multivector_require_equal_dimension( dst , src );

    parallel_for( dst.size() ,
                  functor_type( dst.m_ptr_on_device ,
                                src.m_ptr_on_device ) );
  }
};

/*------------------------------------------------------------------------*/
/** \brief  Copy Host to TPI specialization */
template< typename ValueType , class MapOpt >
class MDArrayDeepCopy< ValueType ,
                       DeviceTPI , MapOpt , true ,
                       DeviceHost , MapOpt , true >
{
public:
  typedef MDArrayView< ValueType , DeviceTPI , MapOpt > dst_type ;
  typedef MDArrayView< ValueType , DeviceHost , MapOpt > src_type ;

  typedef Impl::CopyFunctor< ValueType , DeviceTPI > functor_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};


/** \brief  Copy TPI to Host specialization */
template< typename ValueType , class MapOpt >
class MDArrayDeepCopy< ValueType ,
                       DeviceHost , MapOpt , true ,
                       DeviceTPI , MapOpt , true >
{
public:
  typedef MDArrayView< ValueType , DeviceHost , MapOpt > dst_type ;
  typedef MDArrayView< ValueType , DeviceTPI , MapOpt > src_type ;

  typedef Impl::CopyFunctor< ValueType , DeviceTPI > functor_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

/*------------------------------------------------------------------------*/

} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICE_TPI_DEEP_COPY_HPP */


