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
#include "TSFBlockOperator.hpp"
#include "TSFVectorType.hpp"



using namespace TSFExtended;
using namespace Thyra;
using namespace Teuchos;

template <class Scalar>
class InverseOperator;


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator() : Handle<LinearOpBase<Scalar, Scalar> >() {;}


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator(const RefCountPtr<LinearOpBase<Scalar, Scalar> >& smartPtr) 
  : Handle<LinearOpBase<Scalar, Scalar> >(smartPtr) {;}




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
  this->ptr()->apply(Thyra::NONCONJ_ELE, *(in.ptr().get()),
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
  this->ptr()->applyTranspose(Thyra::NONCONJ_ELE, *(in.ptr().get()),
    out.ptr().get(), alpha, beta);
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
  LinearOperator<Scalar> op = transposedOperator(*this);
  return op;
}

//=======================================================================
template <class Scalar>
LinearOperator<Scalar> 
LinearOperator<Scalar>::inverse(const LinearSolver<Scalar>& solver) const
{
  RefCountPtr<LinearOpBase<Scalar, Scalar> > Ainv 
    = rcp(new InverseOperator<Scalar>(*this, solver));
  LinearOperator<Scalar> op = Ainv;
  return op;
}



//=======================================================================
template <class Scalar>
LinearOperator<Scalar> 
LinearOperator<Scalar>::operator+(const LinearOperator<Scalar>& other) const
{
  LinearOperator<Scalar> op = addedOperator(Array<LinearOperator<Scalar> >(tuple(*this, other)));
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
  TEST_FOR_EXCEPTION(val == 0, std::runtime_error, 
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
  
  TEST_FOR_EXCEPTION(b == 0, std::runtime_error, 
		     "Can't call setBlock since operator not BlockOperator");

  
  b->setNonconstBlock(i, j, sub.ptr());
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
      TEST_FOR_EXCEPTION(i != 0 || j != 0, std::runtime_error, 
                         "nonzero block index (" << i << "," << j << ") into "
                         "non-block operator");
      return *this;
    }
  RefCountPtr<const LinearOpBase<Scalar, Scalar> > block = b->getBlock(i, j);
  RefCountPtr<LinearOpBase<Scalar> > ncBlock = rcp_const_cast<LinearOpBase<Scalar> >(block);
  return ncBlock;
}

 

//=============================================================================
template <class Scalar>
void LinearOperator<Scalar>::endBlockFill() 
{
  Thyra::DefaultBlockedLinearOp<Scalar>* b = 
    dynamic_cast<Thyra::DefaultBlockedLinearOp<Scalar>* >(this->ptr().get());
  
  TEST_FOR_EXCEPTION(b == 0, std::runtime_error, 
		     "Can't call setBlock since operator not BlockOperator");

  
  b->endBlockFill();
} 






#endif
