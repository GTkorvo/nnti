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
#include "TSFCoreOpBase.hpp"
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

// Overridden from Teuchos::Describable

template<class Scalar>
std::ostream& LinearOp<Scalar>::describe(
		std::ostream                         &out
		,const Teuchos::EVerbosityLevel      verbLevel
		,const std::string                   leadingIndent
		,const std::string                   indentSpacer
		) const
{
	typedef Teuchos::ScalarTraits<Scalar> ST;
	const Index dimDomain = this->domain()->dim(), dimRange = this->range()->dim();
	out << leadingIndent << indentSpacer << "type = \'" << this->description()
			<< "\', rangeDim = " << dimRange
			<< ", domainDim = " << dimDomain << "\n";
	if(verbLevel >= Teuchos::VERB_EXTREME) {
		// We will extract by column if op==NOTRANS is supported and by row otherwise
		const ETransp opM = ( this->opSupported(NOTRANS) ? NOTRANS : TRANS );
		// Copy into dense matrix (by column or row)
		Teuchos::RefCountPtr<Vector<Scalar> >
			e_j = ( opM==NOTRANS ? this->domain() : this->range()  )->createMember(),
			t   = ( opM==NOTRANS ? this->range()  : this->domain() )->createMember(); // temp column or row
		const Index
			dimOpMDomain = ( opM==NOTRANS ? dimDomain : dimRange  ),
			dimOpMRange  = ( opM==NOTRANS ? dimRange  : dimDomain );
		RTOpPack::SubVectorT<Scalar> sv;
		std::vector<Scalar>  Md( dimOpMRange * dimOpMDomain ); // Column major
		const Index
			cs = ( opM==NOTRANS ? 1         : dimRange ),  // stride for columns or rows 
			rs = ( opM==NOTRANS ? dimRange  : 1        );  // stride for rows or columns
		Index i, j;
		for( j = 1; j <= dimOpMDomain; ++j ) {
			TSFCore::assign( e_j.get(), ST::zero() );
			TSFCore::set_ele( j, ST::one(), e_j.get() );
			this->apply(opM,*e_j,t.get());  // extract the ith column or row
			t->getSubVector(Range1D(),&sv);
			for( i = 1; i <= dimOpMRange; ++i ) Md[ (i-1)*cs + (j-1)*rs ] = sv(i);
			t->freeSubVector(&sv);
		}
		// Print the matrix
		for( i = 1; i <= dimRange; ++i ) {
			out << leadingIndent << indentSpacer << indentSpacer;
			for( j = 1; j <= dimDomain; ++j )
				out << " " << i << ":" << j << ":" << Md[ (i-1) + (j-1)*dimRange ];
			out << std::endl;
		}
	}
	return out;
}

}	// end namespace TSFCore

#endif // TSFCORE_LINEAR_OP_HPP
