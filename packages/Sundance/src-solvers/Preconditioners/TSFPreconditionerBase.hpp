/* @HEADER@ */
/* ***********************************************************************
// 
//           Playa: Trilinos Solver Framework Extended
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
// **********************************************************************/
/* @HEADER@ */

#ifndef PlayaPRECONDITIONERBASE_HPP
#define PlayaPRECONDITIONERBASE_HPP

#include "SundanceDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   * Base class for preconditioners. A general preconditioner object
   * is split into a left preconditioner M1^-1 and a right
   * preconditioner M2^-1. To solve A x = b, we define the auxiliary
   * system M2^-1 y = x, and solve M1^-1 A M2^-1 y = M1^-1 b to obtain y.
   * Having y, we can quickly recover x by applying M2^-1 to y.
   *
   * The base class implements neither a left nor a right preconditioner.
   */
  template <class Scalar>
  class PreconditionerBase : public Playa::Handleable<PreconditionerBase<Scalar> >
  {
  public:
    /** empty ctor */
    PreconditionerBase() {;}

    /** virtual dtor */
    virtual ~PreconditionerBase(){;}

    
    /** */
    virtual LinearOperator<Scalar> left() const = 0 ;

    /** */
    virtual LinearOperator<Scalar> right() const = 0 ;

    /** return true if this preconditioner has a nontrivial left component */
    virtual bool hasLeft() const = 0 ;

    /** return true if this preconditioner has
     * a nontrivial right component */
    virtual bool hasRight() const = 0 ;

  private:
  };

}

#endif
