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

#ifndef TSFIDENTITYOPERATOR_HPP
#define TSFIDENTITYOPERATOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFSingleScalarTypeOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "TSFHandleable.hpp"
#include "TSFRowAccessibleOp.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace TSFExtended
{
  using namespace Teuchos;
  /** 
   * TSFIdentityOperator is the identity ("I") operator on a vector space.
   */
  template <class Scalar> 
  class IdentityOperator : public SingleScalarTypeOp<Scalar>,
                           public Handleable<SingleScalarTypeOp<Scalar> >,
                           public RowAccessibleOp<Scalar>
  {
  public:
    GET_RCP(SingleScalarTypeOp<Scalar>);

    /** The domain and range spaces for an identity operator
     * are equivalent, so the ctor needs only a single space
     */
    IdentityOperator(const VectorSpace<Scalar>& space)
    : space_(space) {;}
    

    /** Virtual dtor */
    virtual ~IdentityOperator(){;}

    /** 
     * apply returns the input vector
     */
    virtual void generalApply(
                       const Thyra::ETransp            M_trans
                       ,const Thyra::VectorBase<Scalar>    &x
                       ,Thyra::VectorBase<Scalar>          *y
                       ,const Scalar            alpha = 1.0
                       ,const Scalar            beta  = 0.0
                       ) const 
    {
      if (beta == 0.0)
	{
	  if (alpha == 1.0)
	    {
	      assign(y, x);
	      return;
	    }
	  else
	    {
	      assign(y, x);
	      Vt_S(y, alpha);
	      return;
	    }
	}
      else
	{
	  //Thyra::linear_combination(1, &alpha, &&x, beta, y);
	  const Thyra::VectorBase<Scalar>* px = &x;
	  Thyra::linear_combination(1, &alpha, &px, beta, y);
	  return;
	}
    }

    /** Return the domain of the operator */
    RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {return space_.ptr();}
    

    /** Return the range of the operator */
    RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > range() const 
    {return space_.ptr();}

    /** Return the kth row  */
    void getRow(const int& k, 
		Teuchos::Array<int>& indices, 
		Teuchos::Array<Scalar>& values) const
    {
      indices.resize(1);
      indices[0] = k;
      values.resize(1);
      values[0] = 1.0;
    }

  private:
    /** The vector space on which this operator works. Note that the range and
     * domain spaces are identical for the identity operator */
    VectorSpace<Scalar>  space_;
    
  };
}

#endif
