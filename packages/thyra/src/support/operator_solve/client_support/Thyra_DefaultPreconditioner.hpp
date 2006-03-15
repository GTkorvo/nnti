// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_DEFUALT_PRECONDITIONER_HPP
#define THYRA_DEFUALT_PRECONDITIONER_HPP

#include "Thyra_DefaultPreconditionerDecl.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template <class RangeScalar, class DomainScalar>
DefaultPreconditioner<RangeScalar,DomainScalar>::DefaultPreconditioner()
{}

template <class RangeScalar, class DomainScalar>
DefaultPreconditioner<RangeScalar,DomainScalar>::DefaultPreconditioner(
  const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
  ,const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( leftPrecOp.get()==NULL && rightPrecOp.get()==NULL ); 
#endif
  if(leftPrecOp.get()) leftPrecOp_.initialize(leftPrecOp);
  if(rightPrecOp.get()) rightPrecOp_.initialize(rightPrecOp);
}

template <class RangeScalar, class DomainScalar>
DefaultPreconditioner<RangeScalar,DomainScalar>::DefaultPreconditioner(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( leftPrecOp.get()==NULL && rightPrecOp.get()==NULL ); 
#endif
  leftPrecOp_.initialize(leftPrecOp);
  rightPrecOp_.initialize(rightPrecOp);
}

template <class RangeScalar, class DomainScalar>
DefaultPreconditioner<RangeScalar,DomainScalar>::DefaultPreconditioner(
  const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
  )
{
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}

template <class RangeScalar, class DomainScalar>
DefaultPreconditioner<RangeScalar,DomainScalar>::DefaultPreconditioner(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
  )
{
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeLeft(
  const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
  )
{
  uninitialize();
  leftPrecOp_.initialize(leftPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeLeft(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
  )
{
  uninitialize();
  leftPrecOp_.initialize(leftPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeRight(
  const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeRight(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeLeftRight(
  const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
  ,const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
  leftPrecOp_.initialize(leftPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeLeftRight(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &leftPrecOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >   &rightPrecOp
  )
{
  uninitialize();
  rightPrecOp_.initialize(rightPrecOp);
  leftPrecOp_.initialize(leftPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeUnspecified(
  const Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
  )
{
  uninitialize();
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::initializeUnspecified(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &unspecifiedPrecOp
  )
{
  uninitialize();
  unspecifiedPrecOp_.initialize(unspecifiedPrecOp);
}

template <class RangeScalar, class DomainScalar>
void DefaultPreconditioner<RangeScalar,DomainScalar>::uninitialize()
{
  leftPrecOp_.uninitialize();
  rightPrecOp_.uninitialize();
  unspecifiedPrecOp_.uninitialize();
}

// Overridden from PreconditionerBase

template <class RangeScalar, class DomainScalar>
bool DefaultPreconditioner<RangeScalar,DomainScalar>::isLeftPrecOpConst() const
{
  return leftPrecOp_.isConst();
}

template <class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >
DefaultPreconditioner<RangeScalar,DomainScalar>::getNonconstLeftPrecOp()
{
  return leftPrecOp_.getNonconstObj();
}

template <class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >
DefaultPreconditioner<RangeScalar,DomainScalar>::getLeftPrecOp() const
{
  return leftPrecOp_.getConstObj();
}

template <class RangeScalar, class DomainScalar>
bool DefaultPreconditioner<RangeScalar,DomainScalar>::isRightPrecOpConst() const
{
  return rightPrecOp_.isConst();
}

template <class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >
DefaultPreconditioner<RangeScalar,DomainScalar>::getNonconstRightPrecOp()
{
  return rightPrecOp_.getNonconstObj();
}

template <class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >
DefaultPreconditioner<RangeScalar,DomainScalar>::getRightPrecOp() const
{
  return rightPrecOp_.getConstObj();
}

template <class RangeScalar, class DomainScalar>
bool DefaultPreconditioner<RangeScalar,DomainScalar>::isUnspecifiedPrecOpConst() const
{
  return unspecifiedPrecOp_.isConst();
}

template <class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >
DefaultPreconditioner<RangeScalar,DomainScalar>::getNonconstUnspecifiedPrecOp()
{
  return unspecifiedPrecOp_.getNonconstObj();
}

template <class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >
DefaultPreconditioner<RangeScalar,DomainScalar>::getUnspecifiedPrecOp() const
{
  return unspecifiedPrecOp_.getConstObj();
}

} // namespace Thyra

#endif // THYRA_DEFUALT_PRECONDITIONER_HPP
