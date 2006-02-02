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

#ifndef RTOPPACK_ROP_MIN_INDEX_GREATER_THAN_BOUND_HPP
#define RTOPPACK_ROP_MIN_INDEX_GREATER_THAN_BOUND_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace RTOpPack {

/** \brief Returns the minimum element greater than some bound along
 * with its index: <tt>result.scalar = x(k)</tt> and <tt>result.index
 * = k</tt> such that <tt>x(k) <= x(i)</tt> for all <tt>i</tt> where
 * <tt>x(i) > bound</tt> and <tt>k</tt> is the minimum index to break
 * ties.
 *
 * If no element is greater than <tt>bound</tt> then <tt>results.index
 * < 0</tt>.
 *
 * Warning, this class can only be used in serial and SPMD mode as it
 * does not yet support the externalization and internalization of
 * operator object data.
 */
template<class Scalar>
class ROpMinIndexGreaterThanBound : public ROpScalarIndexReductionBase<Scalar> {
public:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, bound )
  /** \brief . */
  ROpMinIndexGreaterThanBound( const Scalar &bound = Teuchos::ScalarTraits<Scalar>::zero() )
    :RTOpT<Scalar>("ROpMinIndexGreaterThanBound")
    ,ROpScalarIndexReductionBase<Scalar>(+Teuchos::ScalarTraits<Scalar>::rmax(),-1)
    ,bound_(bound)
    {}
  /** \brief . */
  ScalarIndex<Scalar> operator()(const ReductTarget& reduct_obj ) const { return this->getRawVal(reduct_obj); }
  /** @name Overridden from RTOpT */
  //@{
  /// This RTOp is NOT coordinate invariant!
  bool coord_invariant() const { return false; }
  /** \brief . */
  void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
    ) const
    {
      const ScalarIndex<Scalar> in    = this->getRawVal(in_reduct_obj);
      const ScalarIndex<Scalar> inout = this->getRawVal(*inout_reduct_obj);
      if( in.scalar < inout.scalar || (in.scalar == inout.scalar && in.index < inout.index) )
        this->setRawVal(in,inout_reduct_obj);
    }
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const
    {
      RTOP_APPLY_OP_1_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      const Index globalOffset = sub_vecs[0].globalOffset();
      ScalarIndex<Scalar> minEle = this->getRawVal(*reduct_obj);
      if( v0_s == 1 ) {
        for( Teuchos_Index i = 1; i <= subDim; ++i ) {
          const Scalar &v0_i = *v0_val++;
          if( v0_i > bound() && (v0_i < minEle.scalar || ( v0_i == minEle.scalar && globalOffset + i < minEle.index ) ) )
            minEle = ScalarIndex<Scalar>(v0_i,globalOffset+i);
        }
      }
      else {
        for( Teuchos_Index i = 1; i <= subDim; ++i, v0_val += v0_s ) {
          const Scalar &v0_i = *v0_val;
          if( v0_i > bound() && (v0_i < minEle.scalar || ( v0_i == minEle.scalar && globalOffset + i < minEle.index ) ) )
            minEle = ScalarIndex<Scalar>(v0_i,globalOffset+i);
        }
      }
      this->setRawVal(minEle,reduct_obj);
    }
  //@}
}; // class ROpMinIndexGreaterThanBound

} // namespace RTOpPack

#endif // RTOPPACK_ROP_MIN_INDEX_GREATER_THAN_BOUND_HPP
