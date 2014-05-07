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

#ifndef ARISTOS_SLEDGEVECTOR_H
#define ARISTOS_SLEDGEVECTOR_H

#include "Aristos_Vector.hpp"
#include "Sledge++.h"

/** \class Aristos::SledgeVector
    \brief The Aristos::Vector / USEME::umDVec (Sledge) adapter class.

    Holds a pointer to an USEME::umDVec and implements the member functions of the Aristos::Vector class.
*/


namespace Aristos {

class SledgeVector : public Vector {

private:

  Teuchos::RefCountPtr<USEME::umDVec>  sledge_vec_;

public:

  SledgeVector( const Teuchos::RefCountPtr<USEME::umDVec> &sledge_vec );

  /** \name Overridden from Vector */
  //@{

  double innerProd( const Vector &x ) const;

  void linComb( const double &alpha, const Vector &x, const double &beta );

  void Scale( const double &alpha );

  void Set( const double &alpha );

  void Set( const double &alpha, const Vector &x );

  Teuchos::RefCountPtr<Vector> createVector() const;

  //@}

  /** Returns a reference counted pointer to the private data container.
  */
  Teuchos::RefCountPtr<const USEME::umDVec> getVector() const;




}; // class SledgeVector

} // namespace Aristos

#endif
