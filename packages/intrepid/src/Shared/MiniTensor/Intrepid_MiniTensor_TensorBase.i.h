// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_TensorBase_i_h)
#define Intrepid_MiniTensor_TensorBase_i_h

namespace Intrepid
{

//
// Get total number of components
//
template<typename T>
inline
Index
TensorBase<T>::get_number_components() const
{
  return components_.size();
}

//
// Allocate space for components
//
template<typename T>
inline
void
TensorBase<T>::set_number_components(Index const number_components)
{
  components_.resize(number_components);

  return;
}

//
// Get dimension
//
template<typename T>
inline
Index
TensorBase<T>::get_dimension() const
{
  return dimension_;
}

//
// Set dimension
//
template<typename T>
inline
void
TensorBase<T>::set_dimension(Index const dimension)
{
  Index const
  order = get_order();

  dimension_ = dimension;

  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  return;
}

//
// Fill components from array defined by pointer.
//
template<typename T>
inline
void
TensorBase<T>::fill(T const * data_ptr)
{
  assert(data_ptr != NULL);

  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = data_ptr[i];
  }

  return;
}

//
// Default constructor
//
template<typename T>
inline
TensorBase<T>::TensorBase() :
dimension_(0)
{
  return;
}

//
// Construction that initializes to NaNs
//
template<typename T>
inline
TensorBase<T>::TensorBase(Index const dimension, Index const order) :
dimension_(dimension)
{
  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = not_a_number<T>();
  }

  return;
}

//
// Construction from a scalar
//
template<typename T>
inline
TensorBase<T>::TensorBase(
    Index const dimension,
    Index const order,
    T const & s) :
dimension_(dimension)
{
  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = s;
  }

  return;
}

//
// Construction from array
//
template<typename T>
inline
TensorBase<T>::TensorBase(
    Index const dimension,
    Index const order,
    T const * data_ptr) :
dimension_(dimension)
{
  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  fill(data_ptr);

  return;
}

//
// Copy constructor
//
template<typename T>
inline
TensorBase<T>::TensorBase(TensorBase<T> const & X) :
dimension_(0)
{
  Index const
  dimension = X.get_dimension();

  dimension_ = dimension;

  Index const
  order = X.get_order();

  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = X[i];
  }

  return;
}

//
// Simple destructor
//
template<typename T>
inline
TensorBase<T>::~TensorBase()
{
  return;
}

//
// Linear access to components
//
template<typename T>
inline
T const &
TensorBase<T>::operator[](Index const i) const
{
  assert(i < components_.size());
  return components_[i];
}

//
// Linear access to components
//
template<typename T>
inline
T &
TensorBase<T>::operator[](Index const i)
{
  assert(i < components_.size());
  return components_[i];
}

//
// Copy assignment
//
template<typename T>
inline
TensorBase<T> &
TensorBase<T>::operator=(TensorBase<T> const & X)
{
  if (this == &X) return *this;

  Index const
  order = X.get_order();

  assert(order == get_order());

  Index const
  dimension = X.get_dimension();

  dimension_ = dimension;

  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = X[i];
  }

  return *this;
}

//
// Component increment
//
template<typename T>
inline
TensorBase<T> &
TensorBase<T>::operator+=(TensorBase<T> const & X)
{
  assert(X.get_order() == get_order());
  assert(X.get_dimension() == get_dimension());

  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] += X[i];
  }

  return *this;
}

//
// Component decrement
//
template<typename T>
inline
TensorBase<T> &
TensorBase<T>::operator-=(TensorBase<T> const & X)
{
  assert(X.get_order() == get_order());
  assert(X.get_dimension() == get_dimension());

  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] -= X[i];
  }

  return *this;
}

//
// Fill with zeros
//
template<typename T>
inline
void
TensorBase<T>::clear()
{
  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = 0.0;
  }

  return;
}

//
// Frobenius norm
//
template<typename T>
T
norm_f(TensorBase<T> const & X)
{
  T
  s = 0.0;

  Index const
  number_components = X.get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    s += X[i] * X[i];
  }

  return std::sqrt(s);
}


} // namespace Intrepid

#endif // Intrepid_MiniTensor_TensorBase_i_h
