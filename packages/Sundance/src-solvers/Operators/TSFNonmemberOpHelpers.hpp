/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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

#ifndef TSFNONMEMBEROPHELPERS_HPP
#define TSFNONMEMBEROPHELPERS_HPP


#include "TSFLinearOperatorImpl.hpp"
#include "TSFMultiVectorOperator.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"


namespace TSFExtended
{

template <class Scalar>
LinearOperator<Scalar> zeroOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > op 
    = rcp(new DefaultZeroLinearOp<Scalar>(range.ptr(), domain.ptr()));

  return op;
}


template <class Scalar>
LinearOperator<Scalar> identityOperator(
  const VectorSpace<Scalar>& space)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > op 
    = rcp(new DefaultIdentityLinearOp<Scalar>(space.ptr()));

  return op;
}


template <class Scalar>
LinearOperator<Scalar> diagonalOperator(
  const Vector<Scalar>& vec)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > op 
    = rcp(new DefaultDiagonalLinearOp<Scalar>(vec.ptr()));

  return op;
}



template <class Scalar>
LinearOperator<Scalar> composedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  Array<RefCountPtr<LinearOpBase<Scalar, Scalar> > > ptrs;

  for (unsigned int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    const MultipliedLinearOpBase<Scalar>* constMultOp 
      = dynamic_cast<const MultipliedLinearOpBase<Scalar>* >(op_i.ptr().get());
    MultipliedLinearOpBase<Scalar>* multOp 
      = dynamic_cast<MultipliedLinearOpBase<Scalar>* >(op_i.ptr().get());
    if (multOp==0 && constMultOp==0) /* is not a MultipliedLinearOpBase */
    {
      ptrs.append(op_i.ptr());
    }
    else if (constMultOp!=0) /* is a const MultipliedLinearOpBase */
    {
      TEST_FOR_EXCEPT(true);
    }
    else if (multOp!=0) /* is a MultipliedLinearOpBase */
    {
      int n = multOp->numOps();
      for (int k=0; k<n; k++)
      {
        RefCountPtr<LinearOpBase<Scalar, Scalar> > A = multOp->getNonconstOp(k);
        ptrs.append(A);
      }
    }
  }
  
  RefCountPtr<LinearOpBase<Scalar, Scalar> > op 
    = rcp(new DefaultMultipliedLinearOp<Scalar>(ptrs.size(), &(ptrs[0])));

  return op;
}




template <class Scalar>
LinearOperator<Scalar> addedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  Array<RefCountPtr<LinearOpBase<Scalar, Scalar> > > ptrs;

  for (unsigned int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    const AddedLinearOpBase<Scalar>* constAddedOp 
      = dynamic_cast<const AddedLinearOpBase<Scalar>* >(op_i.ptr().get());
    AddedLinearOpBase<Scalar>* addedOp 
      = dynamic_cast<AddedLinearOpBase<Scalar>* >(op_i.ptr().get());
    if (addedOp==0 && constAddedOp==0) /* is not a AddedLinearOpBase */
    {
      ptrs.append(op_i.ptr());
    }
    else if (constAddedOp!=0) /* is a const AddedLinearOpBase */
    {
      TEST_FOR_EXCEPT(true);
    }
    else if (addedOp!=0) /* is a AddedLinearOpBase */
    {
      int n = addedOp->numOps();
      for (int k=0; k<n; k++)
      {
        RefCountPtr<LinearOpBase<Scalar, Scalar> > A = addedOp->getNonconstOp(k);
        ptrs.append(A);
      }
    }
  }
  
  RefCountPtr<LinearOpBase<Scalar, Scalar> > op 
    = rcp(new DefaultAddedLinearOp<Scalar>(ptrs.size(), &(ptrs[0])));

  return op;
}



template <class Scalar>
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > A 
    = rcp(new DefaultScaledAdjointLinearOp<Scalar>(scale, Thyra::NOTRANS, op.ptr()));

  return A;
}



template <class Scalar>
LinearOperator<Scalar> scaledTransposedOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > A 
    = rcp(new DefaultScaledAdjointLinearOp<Scalar>(scale, Thyra::TRANS, op.ptr()));

  return A;
}



template <class Scalar>
LinearOperator<Scalar> transposedOperator(
  const LinearOperator<Scalar>& op)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > A
    = rcp(new DefaultScaledAdjointLinearOp<Scalar>(Teuchos::ScalarTraits<Scalar>::one(), Thyra::TRANS, op.ptr()));

  return A;
}

  
template <class Scalar>
LinearOperator<Scalar> multiVectorOperator(
  const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain)
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > A
    = rcp(new MultiVectorOperator<Scalar>(cols, domain));

  return A;
}

 
  
}


#endif
