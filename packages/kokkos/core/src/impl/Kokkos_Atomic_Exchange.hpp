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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_ATOMIC_EXCHANGE_HPP )
#define KOKKOS_ATOMIC_EXCHANGE_HPP

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined( KOKKOS_ATOMICS_USE_CUDA )

KOKKOS_INLINE_FUNCTION
int atomic_exchange( volatile int * const dest , const int val )
{ return atomicExch( (int*) dest , val ); }

KOKKOS_INLINE_FUNCTION
unsigned int atomic_exchange( volatile unsigned int * const dest , const unsigned int val )
{ return atomicExch( (unsigned int*) dest , val ); }

KOKKOS_INLINE_FUNCTION
unsigned long long atomic_exchange( volatile unsigned long long * const dest , const unsigned long long val )
{ return atomicExch( (unsigned long long*) dest , val ); }

template < typename T >
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::UnionPair<T,int,unsigned long long int>::first_type
atomic_exchange( volatile T * const dest , const T val )
{
  typedef Kokkos::Impl::UnionPair<T,int,unsigned long long int> union_type ;
  typedef typename union_type::second_type type ;

  return union_type( atomicExch( (type *) union_type::cast( dest ) ,
                                 union_type::cast( val ) )
                   ).first ;
}

//----------------------------------------------------------------------------

#elif defined(KOKKOS_ATOMICS_USE_GCC) || defined(KOKKOS_ATOMICS_USE_INTEL)

template< typename T >
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::UnionPair<T,int,long>::first_type
atomic_exchange( volatile T * const dest , const T val )
{
  typedef Kokkos::Impl::UnionPair<T,int,long> union_type ;

  union_type assumed , old ;

  old.first = *dest ;
  do {
    assumed.second = old.second ;
    old.second = __sync_val_compare_and_swap( union_type::cast( dest ),
                                              assumed.second ,
                                              union_type::cast( val ) );
  } while ( assumed.second != old.second );

  return old.first ;
}

//----------------------------------------------------------------------------

#elif defined( KOKKOS_ATOMICS_USE_OMP31 )

template < typename T >
KOKKOS_INLINE_FUNCTION
T atomic_exchange( volatile T * const dest , const T val )
{
  T retval;
#pragma omp critical
  {
    retval = dest[0];
    dest[0] = val;
  }
  return retval;
}

#endif

//----------------------------------------------------------------------------

// Simpler version of atomic_exchange when the return value is not needed
template <typename T>
KOKKOS_INLINE_FUNCTION
void atomic_assign(volatile T * const dest, const T src)
{
  atomic_exchange(dest,src);
}

} // namespace Kokkos

#endif

//----------------------------------------------------------------------------

