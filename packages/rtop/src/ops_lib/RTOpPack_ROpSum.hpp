// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
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

#ifndef RTOPPACK_ROP_SUM_HPP
#define RTOPPACK_ROP_SUM_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Summation reduction operator: <tt>result = sum( v0[i], i=1...n )</tt>.
 */
template<class Scalar>
class ROpSum : public ROpScalarReductionBase<Scalar> {
public:
  /** \brief . */
  ROpSum() : RTOpT<Scalar>("ROpSum") {}
  /** \brief . */
  Scalar operator()(const ReductTarget& reduct_obj ) const { return this->getRawVal(reduct_obj); }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
    ,ReductTarget *_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      ReductTargetScalar<Scalar> &reduct_obj = dyn_cast<ReductTargetScalar<Scalar> >(*_reduct_obj); 
      RTOP_APPLY_OP_1_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar sum = reduct_obj.get();
      if( v0_s == 1 ) {
        for( Teuchos_Index i = 0; i < subDim; ++i )
          sum += *v0_val++;
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s )
          sum += *v0_val;
      }
      reduct_obj.set(sum);
    }
  //@}
}; // class ROpSum

} // namespace RTOpPack

#endif // RTOPPACK_ROP_SUM_HPP
