/*
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
*/

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/* TR1 conformal compile-time array traits utilities.
 * Prefer to use TR1 when portably available.
 */

/** \brief  Rank of an array: rank<T>::value */
template< class T > struct rank ;

/** \brief  Extent (a.k.a. dimension) of an array: extent<T,I>::value.
 *          If rank I is not a dimension of the array then value == 0.
 *          If rank I is the leading unspecified dimension [] then value == 0.
 */
template< class T , unsigned I = 0 > struct extent ;

/** \brief  Remove the first (left most) dimension of the type. */
template< class T > struct remove_extent ;

/** \brief  Remove all dimensions of the array */
template< class T > struct remove_all_extents ;

//----------------------------------------------------------------------------

template< class T > struct rank { enum { value = 0 }; };
template< class T > struct rank<T[]> { enum { value = rank<T>::value + 1 }; };
template< class T , unsigned N > struct rank<T[N]>
  { enum { value = rank<T>::value + 1 }; };

//----------------------------------------------------------------------------

template< class T > struct remove_extent { typedef T type ; };
template< class T > struct remove_extent<T[]> { typedef T type ; };
template< class T , unsigned N > struct remove_extent<T[N]>
  { typedef T type ; };

//----------------------------------------------------------------------------

template< class T > struct remove_all_extents { typedef T type ; };
template< class T > struct remove_all_extents<T[]>
  { typedef typename remove_all_extents<T>::type type ; };
template< class T , unsigned N > struct remove_all_extents<T[N]>
  { typedef typename remove_all_extents<T>::type type ; };

//----------------------------------------------------------------------------

template< class T , unsigned I > struct extent         { enum { value = 0 }; };
template< class T >              struct extent<T[], 0> { enum { value = 0 }; };
template< class T , unsigned N > struct extent<T[N],0> { enum { value = N }; };
template< class T , unsigned I > struct extent<T[], I>
  { enum { value = extent<T,I-1>::value }; };
template< class T , unsigned N , unsigned I > struct extent<T[N],I>
  { enum { value = extent<T,I-1>::value }; };

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos


