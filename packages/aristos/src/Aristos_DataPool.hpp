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

#ifndef ARISTOS_DATAPOOL_H
#define ARISTOS_DATAPOOL_H

#include "Aristos_Vector.hpp"

/** \class Aristos::DataPool
    \brief Provides the interface to a generic data pool.

    The only member function is computeAll, responsible for the computation of all data that can be
    stored and reused within one (or more) SQP iteration. The data is entirely problem specific.
*/


namespace Aristos {

class DataPool {
public:

  virtual ~DataPool() {}
  
  /** \brief Recompute all stored quantities.
      \param x [in]  - Current SQP iterate vector.
    
      \return None. 
    
      \par Detailed Description:

      Interface function that evaluates and stores problem-specific quantities
      that can be reused within one (or more) SQP iteration.
  
      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax:\n
      <tt>
      Teuchos::RefCountPtr<const umDVec> ex =
        (Teuchos::dyn_cast<Aristos::SledgeVector>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void computeAll( const Vector &x ) = 0;

}; // class DataPool

} // namespace Aristos

#endif
