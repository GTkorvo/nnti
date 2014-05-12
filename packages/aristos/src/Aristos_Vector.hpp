//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ARISTOS_VECTOR_H
#define ARISTOS_VECTOR_H

#include "Teuchos_RefCountPtr.hpp"

/** \class Aristos::Vector
    \brief Provides the interface to generic abstract vector libraries.

    The interfaced functionality is very basic and includes routines for:\n
    \li linear combinations,
    \li inner products,
    \li scaling and copying operations,
    \li the cloning of a vector.
*/


namespace Aristos {

class Vector {
public:

  virtual ~Vector() {}

  /** \brief Returns inner(*this,x).
  */
  virtual double innerProd( const Vector &x ) const = 0;

  /** \brief <tt>y = alpha*x + beta*y</tt> where <tt>y == *this</tt>.
  */
  virtual void linComb( const double &alpha, const Vector &x, const double &beta = 1.0 ) = 0;

  /** \brief <tt>y = alpha*y</tt> where <tt>y == *this</tt>.
  */
  virtual void Scale( const double &alpha ) = 0;

  /** \brief <tt>y = alpha</tt> where <tt>y == *this</tt>.
  */
  virtual void Set( const double &alpha ) = 0;

  /** \brief <tt>y = alpha*x</tt> where <tt>y == *this</tt>.
  */
  virtual void Set( const double &alpha, const Vector &x ) = 0;

  /** Clone to make a new (uninitialized) vector.
  */
  virtual Teuchos::RefCountPtr<Vector> createVector() const = 0;

}; // class Vector

} // namespace Aristos

#endif
