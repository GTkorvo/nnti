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

#ifndef TSFLINEAROPERATORIMPL_HPP
#define TSFLINEAROPERATORIMPL_HPP

#include "TSFConfigDefs.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFInverseOperator.hpp"
#include "TSFTransposeOperator.hpp"
#include "TSFComposedOperator.hpp"
#include "TSFSumOperator.hpp"
#include "TSFScaledOperator.hpp"
#include "TSFBlockOperator.hpp"
#include "TSFVectorType.hpp"



using namespace TSFExtended;
using namespace Thyra;
using namespace Teuchos;

template <class Scalar>
class InverseOperator;


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator() : Handle<SingleScalarTypeOpBase<Scalar> >() {;}



//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator(Handleable<SingleScalarTypeOpBase<Scalar> >* rawPtr) 
  : Handle<SingleScalarTypeOpBase<Scalar> >(rawPtr) {;}

//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator(const RefCountPtr<SingleScalarTypeOpBase<Scalar> >& smartPtr) 
  : Handle<SingleScalarTypeOpBase<Scalar> >(smartPtr) {;}




//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>::apply(const Vector<Scalar>& in,
                                   Vector<Scalar>& out,
                                   const Scalar& alpha,
                                   const Scalar& beta) const
{
  /* the result vector might not be initialized. If it's null,
   * create a new vector in the range space */
  if (out.ptr().get()==0)
    {
      out = this->range().createMember();
    }
  this->ptr()->generalApply(Thyra::NOTRANS, *(in.ptr().get()),
                            out.ptr().get(), alpha, beta);
}




//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>::applyTranspose(const Vector<Scalar>& in,
                                            Vector<Scalar>& out,
                                            const Scalar& alpha,
                                            const Scalar& beta) const
{
  /* the result vector might not be initialized. If it's null,
   * create a new vector in the domain space (i.e., the range space
   * of the transpose operator */
  if (out.ptr().get()==0)
    {
      out = this->domain().createMember();
    }
  this->ptr()->generalApply(Thyra::TRANS, *(in.ptr().get()),
                            out.ptr().get(), alpha, beta);
}


//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>
::generalApply(
               const Thyra::ETransp            M_trans
               ,const Thyra::VectorBase<Scalar>    &x
               ,Thyra::VectorBase<Scalar>          *y
               ,const Scalar            alpha
               ,const Scalar            beta
               ) const 
{
  this->ptr()->generalApply(M_trans, x, y, alpha, beta);
}


//=======================================================================
template <class Scalar>
RefCountPtr<Time>& LinearOperator<Scalar>::opTimer()
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("Low-level vector operations");
  return rtn;
}

//=======================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::transpose() const
{
  LinearOperator<Scalar> op = new TransposeOperator<Scalar>(*this);
  return op;
}

//=======================================================================
template <class Scalar>
LinearOperator<Scalar> 
LinearOperator<Scalar>::inverse(const LinearSolver<Scalar>& solver) const
{
  LinearOperator<Scalar> op = new InverseOperator<Scalar>(*this, solver);
  return op;
}



//=======================================================================
template <class Scalar>
LinearOperator<Scalar> 
LinearOperator<Scalar>::operator+(const LinearOperator<Scalar>& other) const
{
  LinearOperator<Scalar> op = new SumOperator<Scalar>(*this, other);
  return op;
}



//=======================================================================
template <class Scalar>
RefCountPtr<LoadableMatrix<Scalar> > LinearOperator<Scalar>::matrix()
{
  RefCountPtr<LoadableMatrix<Scalar> > rtn 
    = rcp_dynamic_cast<LoadableMatrix<Scalar> >(this->ptr());
  return rtn;
}

//=======================================================================
template <class Scalar>
void LinearOperator<Scalar>::getRow(const int& row, 
				    Teuchos::Array<int>& indices, 
				    Teuchos::Array<Scalar>& values) const
{
  const RowAccessibleOp<Scalar>* val = 
    dynamic_cast<const RowAccessibleOp<Scalar>* >(this->ptr().get());
  TEST_FOR_EXCEPTION(val == 0, runtime_error, 
		     "Operator not row accessible; getRow() not defined.");
  val->getRow(row, indices, values);
}

//=============================================================================
template <class Scalar>
int LinearOperator<Scalar>::numBlockRows() const
{
  const BlockOperator<Scalar>* b = dynamic_cast<const BlockOperator<Scalar>* >(this->ptr().get());
  if (b==0) return 1;
  return b->numBlockRows(); 
}

//=============================================================================
template <class Scalar>
int LinearOperator<Scalar>::numBlockCols() const
{
  const BlockOperator<Scalar>* b = dynamic_cast<const BlockOperator<Scalar>* >(this->ptr().get());
  if (b==0) return 1;
  return b->numBlockCols(); 
}


//=============================================================================
template <class Scalar>
const VectorSpace<Scalar> 
LinearOperator<Scalar>::range() const
{return this->ptr()->range();}
  

//=============================================================================
template <class Scalar>
void LinearOperator<Scalar>::setBlock(int i, int j, 
				      const LinearOperator<Scalar>& sub) 
{
  BlockOperator<Scalar>* b = 
    dynamic_cast<BlockOperator<Scalar>* >(this->ptr().get());
  
  TEST_FOR_EXCEPTION(b == 0, runtime_error, 
		     "Can't call setBlock since operator not BlockOperator");

  
  b->setBlock(i, j, sub.ptr());
} 



//=============================================================================
template <class Scalar>
const  VectorSpace<Scalar> 
LinearOperator<Scalar>::domain() const 
{return this->ptr()->domain();}



//=============================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::getBlock(const int &i, 
							const int &j) const 
{
  BlockOperator<Scalar>* b = 
    dynamic_cast<BlockOperator<Scalar>* >(this->ptr().get());
  
  if (b==0)
    {
      TEST_FOR_EXCEPTION(i != 0 || j != 0, runtime_error, 
                         "nonzero block index (" << i << "," << j << ") into "
                         "non-block operator");
      return *this;
    }
  RefCountPtr<LinearOpBase<Scalar, Scalar> > block = b->getBlock(i, j);
  return rcp_dynamic_cast<SingleScalarTypeOpBase<Scalar> >(block);
}






#endif
