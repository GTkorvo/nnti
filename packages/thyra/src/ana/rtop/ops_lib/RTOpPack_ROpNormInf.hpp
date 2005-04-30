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

#ifndef RTOPPACK_ROP_NORMINF_HPP
#define RTOPPACK_ROP_NORMINF_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Infinity norm reduction operator: <tt>result = sum( |v0[i]|, i=1...n )</tt>.
 */
template<class Scalar>
class ROpNormInf : public ROpScalarReductionBase<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> {
public:
  /** \brief . */
  ROpNormInf() : RTOpT<Scalar>("ROpNormInf") {}
  /** \brief . */
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  operator()(const ReductTarget& reduct_obj ) const { return this->getRawVal(reduct_obj); }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void reduce_reduct_objs(
    const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      const ReductTargetScalar<Scalar>
        &in_reduct_obj = dyn_cast<const ReductTargetScalar<Scalar> >(_in_reduct_obj);
      ReductTargetScalar<Scalar>
        &inout_reduct_obj = dyn_cast<ReductTargetScalar<Scalar> >(*_inout_reduct_obj);
      if( in_reduct_obj.get() > inout_reduct_obj.get() ) inout_reduct_obj.set(in_reduct_obj.get());
    }
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
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm_inf = reduct_obj.get();
      if( v0_s == 1 ) {
        for( RTOp_index_type i = 0; i < subDim; ++i ) {
          const typename Teuchos::ScalarTraits<Scalar>::magnitudeType
            mag = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val++);
          norm_inf = mag > norm_inf ? mag : norm_inf;
        }
      }
      else {
        for( RTOp_index_type i = 0; i < subDim; ++i, v0_val += v0_s ) {
          const typename Teuchos::ScalarTraits<Scalar>::magnitudeType
            mag = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val);
          norm_inf = mag > norm_inf ? mag : norm_inf;
        }
      }
      reduct_obj.set(norm_inf);
    }
  //@}
}; // class ROpNormInf

} // namespace RTOpPack

#endif // RTOPPACK_ROP_NORMINF_HPP
