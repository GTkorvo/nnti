// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_STATIC_ARRAY_TRAITS_HPP
#define STOKHOS_STATIC_ARRAY_TRAITS_HPP

#include <cstring>

#include "Sacado_Traits.hpp"

#include "KokkosArray_Macros.hpp"

namespace Stokhos {

  /*!
   * \brief Static array allocation class
   */
  template <typename T, typename node, 
	    bool isScalar = Sacado::IsScalarType<T>::value>
  struct StaticArrayTraits {};

  /*!
   * \brief Static array allocation class that works for any type
   */
  template <typename T, typename N>
  struct StaticArrayTraits<T, N, false> {

    typedef T value_type;
    typedef N node_type;
    
    //! Copy array from \c src to \c dest of length \c sz
    static
    KOKKOSARRAY_INLINE_FUNCTION
    void copy(const T* src, T*  dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    static 
    KOKKOSARRAY_INLINE_FUNCTION
    void zero(T* dest, std::size_t sz) {
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = T(0.);
    }

    //! Fill array \c dest of length \c sz with value \c v
    static 
    KOKKOSARRAY_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, const T& v) {
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = v;
    }

  };

  /*!
   * \brief Static array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T, typename N>
  struct StaticArrayTraits<T,N,true> {

    typedef T value_type;
    typedef N node_type;
    
    //! Copy array from \c src to \c dest of length \c sz
    static 
    KOKKOSARRAY_INLINE_FUNCTION
    void copy(const T* src, T* dest, std::size_t sz) {
      if (sz > 0) std::memcpy(dest,src,sz*sizeof(T));
    }
    
    //! Zero out array \c dest of length \c sz
    static 
    KOKKOSARRAY_INLINE_FUNCTION
    void zero(T* dest, std::size_t sz) {
      if (sz > 0) std::memset(dest,0,sz*sizeof(T));
    }

    //! Fill array \c dest of length \c sz with value \c v
    static 
    KOKKOSARRAY_INLINE_FUNCTION
    void fill(T* dest, std::size_t sz, T v) {
      //std::memset(dest,v,sz*sizeof(T)); // memset doesn't work if v != 0?
      for (std::size_t i=0; i<sz; ++i)
	*(dest++) = v;
    }

  };

} // namespace Stokhos

#endif // STOKHOS_STATIC_ARRAY_TRAITS_HPP
