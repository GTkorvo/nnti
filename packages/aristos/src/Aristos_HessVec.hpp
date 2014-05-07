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

#ifndef ARISTOS_HESSVEC_H
#define ARISTOS_HESSVEC_H

#include "Aristos_Vector.hpp"

/** \class Aristos::HessVec
    \brief Provides the interface for the evaluation of Hessian-of-Lagrangian times vector products.
*/


namespace Aristos {

class HessVec {
public:

  virtual ~HessVec() {}
  
  /** \brief Returns Hessian-of-Lagrangian times vector product.

      \param x [in]     - Current SQP iterate vector.
      \param l [in]     - The vector of Lagrange multipliers.
      \param s [in]     - The vector to which the Hessian-of-Lagrangian operator is applied.
      \param Hs [out]   - Resulting vector.

      \return None.

      \par Detailed Description:

      Interface function that evaluates the vector obtained by applying the Hessian-of-Lagrangian operator
      to a vector. To be subclassed and implemented by the user.

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void getValue( const Vector &x, const Vector &l, const Vector &s, Vector &Hs) const = 0;

}; // class HessVec

} // namespace Aristos

#endif
