// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
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

// ///////////////////////////////
// RTOpPack_TOpRandomize.hpp

#ifndef RTOPPACK_TOP_RANDOMIZE_HPP
#define RTOPPACK_TOP_RANDOMIZE_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

///
/** Generate a random vector in the range [l,u]: <tt>z0[i] = 0.5*((u-l)*Teuchos::ScalarTraits<Scalar>::random()+(u+l)), i=1...n</tt>.
 *
 * The seed for the random number generator can be set by:
 \code
  Teuchos::ScalarTraits<Scalar>::seedrandom(s)
 \endcode
 * where <tt>s</tt> is some unsigned integer
 */
template<class Scalar>
class TOpRandomize : public ROpScalarScalarTransformationBase<Scalar> {
public:
  ///
  void set_bounds( const Scalar& l, const Scalar& u ) { this->scalarData1(l); this->scalarData2(u); }
  ///
  TOpRandomize(
    const Scalar& l   = -Teuchos::ScalarTraits<Scalar>::one()
    ,const Scalar& u  = +Teuchos::ScalarTraits<Scalar>::one()
    )
    : RTOpT<Scalar>("TOpRandomize"), ROpScalarScalarTransformationBase<Scalar>(l,u)
    {}
  /** @name Overridden from RTOpT */
  //@{
  ///
	void apply_op(
		const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
		,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
		,ReductTarget *reduct_obj
		) const
    {
      const Scalar l = this->scalarData1(), u = this->scalarData2();
      const Scalar a = Scalar(0.5)*(u-l), b = Scalar(0.5)*(u+l) ; // Linear coefficients for translating from [-1,+1] to [l,b]
      RTOP_APPLY_OP_0_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      for( RTOp_index_type i = 0; i < subDim; ++i, z0_val += z0_s ) {
        *z0_val = a * Teuchos::ScalarTraits<Scalar>::random() + b; // Should be in the range [l,b]
      }
    }
  //@}
}; // class TOpRandomize

} // namespace RTOpPack

#endif // RTOPPACK_TOP_RANDOMIZE_HPP
