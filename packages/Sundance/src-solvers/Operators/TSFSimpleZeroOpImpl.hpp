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

#ifndef Playa_SIMPLE_ZERO_OP_IMPL_HPP
#define Playa_SIMPLE_ZERO_OP_IMPL_HPP



#include "PlayaSimpleZeroOpDecl.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimplifiedLinearOpBaseImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;
using namespace Sundance;




/* ---- Zero op ------- */

template <class Scalar> inline
SimpleZeroOp<Scalar>::SimpleZeroOp(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : SimplifiedLinearOpWithSpaces<Scalar>(domain, range) {}



template <class Scalar> inline
void SimpleZeroOp<Scalar>::applyOp(const Thyra::EOpTransp M_trans,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(this->verb(), tab << "SimpleZeroOp::applyOp()");

  out.zero();

  SUNDANCE_MSG2(this->verb(), tab << "done SimpleZeroOp::applyOp()");
}

/* */
template <class Scalar> inline
std::string SimpleZeroOp<Scalar>::description() const 
{return "ZeroOp(domain=" 
    + this->domain()->description() 
    + ", range=" + this->range()->description() + ")";}



template <class Scalar> inline
LinearOperator<Scalar> zeroOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
{
  RCP<LinearOpBase<Scalar> > op 
    = rcp(new SimpleZeroOp<Scalar>(domain, range));

  return op;
}


}

#endif
