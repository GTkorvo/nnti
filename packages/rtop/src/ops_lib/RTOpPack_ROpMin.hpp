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

#ifndef RTOPPACK_ROP_MIN_HPP
#define RTOPPACK_ROP_MIN_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Returns the minimum element: <tt>result = min{ v0[i], i=0...n-1 }</tt>.
 */
template<class Scalar>
class ROpMin : public ROpScalarReductionBase<Scalar> {
public:
  /** \brief . */
  ROpMin() : RTOpT<Scalar>("ROpMin"), ROpScalarReductionBase<Scalar>(+Teuchos::ScalarTraits<Scalar>::rmax()) {}
  /** \brief . */
  Scalar operator()(const ReductTarget& reduct_obj ) const { return this->getRawVal(reduct_obj); }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
    ) const
    {
      const Scalar in_min_ele    = this->getRawVal(in_reduct_obj);
      const Scalar inout_min_ele = this->getRawVal(*inout_reduct_obj);
      this->setRawVal( in_min_ele < inout_min_ele ? in_min_ele : inout_min_ele, inout_reduct_obj );
    }
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const
    {
      RTOP_APPLY_OP_1_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar min_ele = this->getRawVal(*reduct_obj);
      if( v0_s == 1 ) {
        for( Teuchos_Index i = 0; i < subDim; ++i ) {
          const Scalar &v0_i = *v0_val++;
          min_ele = ( v0_i < min_ele ? v0_i : min_ele );
        }
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s ) {
          const Scalar &v0_i = *v0_val;
          min_ele = ( v0_i < min_ele ? v0_i : min_ele );
        }
      }
      this->setRawVal(min_ele,reduct_obj);
    }
  //@}
}; // class ROpMin

} // namespace RTOpPack

#endif // RTOPPACK_ROP_MIN_HPP
