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

#ifndef TSFMULTIVECTOROPERATOR_HPP
#define TSFMULTIVECTOROPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFSingleScalarTypeOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "TSFRowAccessibleOp.hpp"
#include "TSFHandleable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFVectorImpl.hpp"

namespace TSFExtended
{
  /** 
   * A MultiVectorOperator is a linear operator whose
   * rows or columns are represented as a multivector 
   */
  template <class Scalar> 
  class MultiVectorOperator 
    : public SingleScalarTypeOp<Scalar>,
      public Handleable<SingleScalarTypeOpBase<Scalar> >,
      public RowAccessibleOp<Scalar>
  {
  public:
    GET_RCP(SingleScalarTypeOpBase<Scalar>);

    /**
     * Construct from an array of vectors and a specifier for the 
     * domain space. 
     */
    MultiVectorOperator(const Teuchos::Array<Vector<Scalar> >& cols,
      const VectorSpace<Scalar>& domain)
      : cols_(cols),
        domain_(domain.ptr()),
        range_()
    {
      TEST_FOR_EXCEPTION(cols.size() == 0U, runtime_error,
        "empty multivector given to MultiVectorOperator ctor");
      range_ = cols[0].space();
      for (unsigned int i=1; i<cols.size(); i++)
      {
        TEST_FOR_EXCEPTION(cols[i].space() != range_, runtime_error,
          "inconsistent vector spaces in  MultiVectorOperator ctor");
      }
    }

    /** Virtual dtor */
    virtual ~MultiVectorOperator(){;}


    /** 
     * Apply does an element-by-element multiply between the input 
     * vector, x, and the diagonal values.
     */
    virtual void generalApply(
			      const Thyra::ETransp            M_trans
			      ,const Thyra::VectorBase<Scalar>    &x
			      ,Thyra::VectorBase<Scalar>          *y
			      ,const Scalar            alpha = 1.0
			      ,const Scalar            beta  = 0.0
			      ) const 
    {
      if (M_trans == Thyra::NOTRANS)
      {
        Vector<Scalar> vx 
          = rcp(const_cast<Thyra::VectorBase<Scalar>*>(&x), false);
        Vector<Scalar> vy = rcp(y, false);

        if (beta != 0.0) vy.scale(beta);
        else vy.zero();

        for (unsigned int i=0; i<cols_.size(); i++)
        {
          vy.update(alpha * vx.getElement(i), cols_[i]);
        }
      }
      else
      {
        Vector<Scalar> vx 
          = rcp(const_cast<Thyra::VectorBase<Scalar>*>(&x), false);
        Vector<Scalar> vy = rcp(y, false);

        if (beta != 0.0) vy.scale(beta);
        else vy.zero();

        for (unsigned int i=0; i<cols_.size(); i++)
        {
          vy.addToElement(i, alpha * vx.dot(cols_[i]));
        }
      }
    }



    /** Return the domain of the operator */
    RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {return domain_.ptr();}
 



    /** Return the range of the operator */
    RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > range() const 
    {return range_.ptr();}



    /** Return the kth row  */
    void getRow(const int& k, 
                Teuchos::Array<int>& indices, 
                Teuchos::Array<Scalar>& values) const
    {
      indices.resize(cols_.size());
      values.resize(cols_.size());
      for (unsigned int j=0; j<cols_.size(); j++)
      {
        indices[j] = j;
        values[j] = cols_[j].getElement(k);
      }
    }



  private:

    Teuchos::Array<Vector<Scalar> > cols_;
    VectorSpace<Scalar> domain_;
    VectorSpace<Scalar> range_;
  };
}

#endif
