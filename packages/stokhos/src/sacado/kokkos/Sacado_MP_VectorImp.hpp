// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

namespace Sacado {
namespace MP {

template <typename T, typename Storage> 
Vector<T,Storage>::
Vector() :
  s(1)
{
}

template <typename T, typename Storage> 
Vector<T,Storage>::
Vector(const typename Vector<T,Storage>::value_type& x) :
  s(1)
{
  s.init(x);
}

template <typename T, typename Storage> 
Vector<T,Storage>::
Vector(ordinal_type sz, const value_type& x) :
  s(sz,x)
{
}

template <typename T, typename Storage> 
Vector<T,Storage>::
Vector(const Vector<T,Storage>& x) :
  s(x.s)
{
}

template <typename T, typename Storage> 
template <typename S>
Vector<T,Storage>::
Vector(const Expr<S>& xx) : s(xx.derived().size())
{
  const typename Expr<S>::derived_type& x = xx.derived();

  if (x.hasFastAccess(s.size())) {
    for (int i=0; i<s.size(); i++)
      s[i] = x.fastAccessCoeff(i);
  }
  else {
    for (int i=0; i<s.size(); i++)
      s[i] = x.coeff(i);
  }
}

template <typename T, typename Storage> 
void
Vector<T,Storage>::
reset(ordinal_type sz_new)
{
  int sz = this->size();
  s.resize(sz_new);
  if (sz == 1 && sz_new > sz)
    for (int i=1; i<sz_new; i++)
      s[i] = s[0];
}

template <typename T, typename Storage> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator=(const typename Vector<T,Storage>::value_type& v) 
{
  s.init(v);
  return *this;
}

template <typename T, typename Storage> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator=(const Vector<T,Storage>& x) 
{
  s = x.s;
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator=(const Expr<S>& xx) 
{
  const typename Expr<S>::derived_type& x = xx.derived();

  this->reset(x.size());
  if (x.hasFastAccess(s.size())) {
    for (int i=0; i<s.size(); i++)
      s[i] = x.fastAccessCoeff(i);
  }
  else {
    for (int i=0; i<s.size(); i++)
      s[i] = x.coeff(i);
  }
  return *this;
}

template <typename T, typename Storage> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator+=(const typename Vector<T,Storage>::value_type& v)
{
  for (int i=0; i<s.size(); i++)
    s[i] += v;
  return *this;
}

template <typename T, typename Storage> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator-=(const typename Vector<T,Storage>::value_type& v)
{
  for (int i=0; i<s.size(); i++)
    s[i] -= v;
  return *this;
}

template <typename T, typename Storage> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator*=(const typename Vector<T,Storage>::value_type& v)
{
  for (int i=0; i<s.size(); i++)
    s[i] *= v;
  return *this;
}

template <typename T, typename Storage> 
Vector<T,Storage>& 
Vector<T,Storage>::
operator/=(const typename Vector<T,Storage>::value_type& v)
{
  for (int i=0; i<s.size(); i++)
    s[i] /= v;
  return *this;
}

template <typename T, typename Storage>
std::ostream& 
operator << (std::ostream& os, const Vector<T,Storage>& a)
{
  typedef typename Vector<T,Storage>::ordinal_type ordinal_type;

  os << "[ ";
      
  for (ordinal_type i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

} // namespace MP
} // namespace Sacado
