// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// /////////////////////////////////////////////////////////////
// Teuchos_arrayArg.hpp

#ifndef TEUCHOS_ARRAY_ARG_HPP
#define TEUCHOS_ARRAY_ARG_HPP

/*! \file Teuchos_arrayArg.hpp
    \brief Utility that allows arrays to be passed into argument list
*/

#include "Teuchos_TestForException.hpp"

namespace Teuchos {

/** \defgroup Teuchos_Array_Arguments Utility functions for passing arrays into argument lists.

\brief The purpose of this utility is to make passing arrays into argument lists easier.

Declaring arrays outside of a function just to pass a (small) list of
values into a function can be tiresome. The templated function
<tt>arrayArg()</tt> simplifies this process.  With this function you
can construct (using stack memory not dynamically allocated memory) an
array of data to be passed into a function.

For example, consider the following function prototype:

\code
void f( const int x_size, const double x[] );
\endcode

which takes a <tt>const</tt> array of <tt>double</tt>s of length
<tt>x_size</tt>.  Generally, to call this function one would have to
first declare an array and then call the function as:

\code
void f()
{
  ...
  double x[] = { 1.0, 2.0, 3.0 };
  f( 3, x );
  ...
\endcode

Now, however, one can create the array in the call to <tt>f()</tt> as:

\code
void f()
{
  ...
  f( 3, arrayArg<double>(1.0,2.0,3.0)() );
  ...
}
\endcode

In the above situation, one may be able to write the call as:

\code
void f()
{
  ...
  f( 3, arrayArg(1.0,2.0,3.0) );
  ...
}
\endcode

but the former, slightly more verbose, version is to be preferred
since it makes explicit what type of array is being created and insures
that the compiler will not get confused about the final implicit conversion
to a raw <tt>const double*</tt> pointer.

The <tt>arrayArg()</tt> function is overloaded to accept 1, 2, 3, 4, 5 and 6 
arguments.  If more elements are needed, then more overrides are easy to add.

The utility functions <tt>arrayArg()</tt> can only be used
to pass in arrays of <tt>const</tt> objects.  Using this utility for
non-<tt>const</tt> objects is generally not possible due to the fact that
compilers are not allowed to pass in compiler-generated objects into
a function through non-<tt>const</tt> reference arguments.  For example,
to pass in a non-<tt>const</tt> array for the function:

\code
void f2( const int y_size, double y[] );
\endcode

one would have to call it as:

\code
void g2()
{
  double y[3]; // Output argument
  f2(3,y );
}
\endcode

*/

///
/** \brief Utility class that allows arrays to be passed into argument list.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<int N, class T>
class ArrayArg {
public:
  /// Basic constructor taking a copy of the \c array of length \c N
  ArrayArg( const T array[] ) { std::copy( array, array+N, array_ ); }

  /// Return a \c const pointer to the internal array
  const T* operator()() { return array_; }

  /// Return a \c const pointer to the internal array
  operator const T* () { return array_; }

private:
  T array_[N]; //  Can't be a const array!
}; // class Array1DArg

///
/** \brief Return an array with 1 member.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline ArrayArg<1,T> arrayArg( const T &t1 )
{
  const T array[] = { t1 };
  return ArrayArg<1,T>(array);
}

///
/** \brief Return an array with 2 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline ArrayArg<2,T> arrayArg( const T &t1, const T &t2 )
{
  const T array[] = { t1, t2 };
  return ArrayArg<2,T>(array);
}

///
/** \brief Return an array with 3 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline ArrayArg<3,T> arrayArg( const T &t1, const T &t2, const T &t3 )
{
  const T array[] = { t1, t2, t3 };
  return ArrayArg<3,T>(array);
}

///
/** \brief Return an array with 4 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline ArrayArg<4,T> arrayArg( const T &t1, const T &t2, const T &t3, const T &t4 )
{
  const T array[] = { t1, t2, t3, t4 };
  return ArrayArg<4,T>(array);
}

///
/** \brief Return an array with 5 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline ArrayArg<5,T> arrayArg( const T &t1, const T &t2, const T &t3, const T &t4, const T &t5 )
{
  const T array[] = { t1, t2, t3, t4, t5 };
  return ArrayArg<5,T>(array);
}

///
/** \brief Return an array with 6 members.
 *
 * \ingroup Teuchos_Array_Arguments
 */
template<class T>
inline ArrayArg<6,T> arrayArg( const T &t1, const T &t2, const T &t3, const T &t4, const T &t5, const T &t6 )
{
  const T array[] = { t1, t2, t3, t4, t5, t6 };
  return ArrayArg<6,T>(array);
}

} // namespace Teuchos

#endif // TEUCHOS_ARRAY_ARG_HPP
