// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_BASE_HPP
#define THYRA_MULTI_VECTOR_BASE_HPP

#include "Thyra_MultiVectorBase_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorStdOps_decl.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

#include "RTOpPack_TOpAbs.hpp"

namespace Thyra {


// Provide access to the columns as VectorBase objects


template<class Scalar>
RCP<const VectorBase<Scalar> >
MultiVectorBase<Scalar>::colImpl(Ordinal j) const
{
  return const_cast<MultiVectorBase*>(this)->nonconstColImpl(j);
}


// Overridden methods from LinearOpBase


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
MultiVectorBase<Scalar>::clone() const
{
  return this->clone_mv();
}

// Overridden methods from RowStatLinearOpBase

template<class Scalar>
bool MultiVectorBase<Scalar>::
rowStatIsSupportedImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat) const
{
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      return true;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  return false; // will never be called
}

template<class Scalar>
void MultiVectorBase<Scalar>::
getRowStatImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat,
               const Ptr<VectorBase<Scalar> > &rowStatVec) const
{
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
      absRowSum(rowStatVec);
      ::Thyra::reciprocal<Scalar>(*rowStatVec,rowStatVec.ptr());
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      // compute absolute row sum
      absRowSum(rowStatVec);
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
      absColSum(rowStatVec);
      ::Thyra::reciprocal<Scalar>(*rowStatVec,rowStatVec.ptr());
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      // compute absolute row sum
      absColSum(rowStatVec);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}

// helper methods

template<class Scalar>
void MultiVectorBase<Scalar>::
absRowSum(const Teuchos::Ptr<Thyra::VectorBase<Scalar> > & output) const
{
  using Teuchos::RCP;
  using Teuchos::ptrFromRef;
  using Teuchos::tuple;

  // compute absolute value of multi-vector
  RTOpPack::TOpAbs<Scalar> abs_op;
  RCP<MultiVectorBase<Scalar> > abs_mv = createMembers(this->range(),this->domain());
  ::Thyra::applyOp<Scalar>( abs_op, tuple(ptrFromRef(*this)), tuple(abs_mv.ptr()), Teuchos::null );

  // compute sum over all rows
  RCP<VectorBase<Scalar> > ones = Thyra::createMember(this->domain());
  ::Thyra::put_scalar<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),ones.ptr());
  ::Thyra::apply<Scalar>(*abs_mv,Thyra::NOTRANS,*ones,output);
}

template<class Scalar>
void MultiVectorBase<Scalar>::
absColSum(const Teuchos::Ptr<Thyra::VectorBase<Scalar> > & output) const
{ 
  RTOpPack::SubVectorView<Scalar> view;
  output->acquireDetachedView(Thyra::Range1D(),&view);

  Thyra::norms_1<Scalar>(*this,view.values()());
  
  output->commitDetachedView(&view);
}


} // end namespace Thyra


#endif // THYRA_MULTI_VECTOR_BASE_HPP
