// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreLinearOp.hpp

#ifndef TSFCORE_LINEAR_OP_HPP
#define TSFCORE_LINEAR_OP_HPP

#include "TSFCoreLinearOpDecl.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreAssertOp.hpp"

namespace TSFCore {

// Virtual functions with default implemenations

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> > 
LinearOp<Scalar>::clone() const
{
	return Teuchos::null;
}

template<class Scalar>
void LinearOp<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	const VectorSpace<Scalar> &space_mv_rows = *Y->domain();
	const Index               num_mv_cols    = space_mv_rows.dim();
	for( Index j = 1; j <= num_mv_cols; ++j )
		this->apply(M_trans,*X.col(j),Y->col(j).get(),alpha,beta);
}

}	// end namespace TSFCore

#endif // TSFCORE_LINEAR_OP_HPP
