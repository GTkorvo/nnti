// ////////////////////////////////////////////////////////////////////////////////////////////////////
// MatrixIdentConcat.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include "ConstrainedOptimizationPack/include/MatrixIdentConcat.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/AbstractLinAlgPackAssertOp.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StV;
}

namespace {

// Get a view of a vector (two versions)

inline
ReferenceCountingPack::ref_count_ptr<const AbstractLinAlgPack::VectorWithOp>
get_view(
	const AbstractLinAlgPack::VectorWithOp  &v
	,const RangePack::Range1D               &rng
	)
{
 	return v.sub_view(rng);
}

inline
ReferenceCountingPack::ref_count_ptr<const AbstractLinAlgPack::SpVectorSlice>
get_view(
	const AbstractLinAlgPack::SpVectorSlice &v
	,const RangePack::Range1D               &rng
	)
{
	return ReferenceCountingPack::ref_count_ptr<const AbstractLinAlgPack::SpVectorSlice>(
		new AbstractLinAlgPack::SpVectorSlice( v(rng) ) );
}

// Matrix-vector multiplication

template<class V>
void mat_vec(
	AbstractLinAlgPack::VectorWithOpMutable        *y
	,AbstractLinAlgPack::value_type                a
	,AbstractLinAlgPack::value_type                alpha
	,const AbstractLinAlgPack::MatrixWithOp        &D
	,BLAS_Cpp::Transp                              D_trans
	,const LinAlgPack::Range1D                     &D_rng
	,const LinAlgPack::Range1D                     &I_rng
	,BLAS_Cpp::Transp                              M_trans
	,const V                                       &x
	,AbstractLinAlgPack::value_type                b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;

	//
	// M = [ alpha*op(D) ] or [      I       ]
	//     [     I       ]    [ alpha*op(D)  ]
	//
	if( M_trans == no_trans ) {
		// 
		// y = b*y + a*M*x
		// =>
		// y(D_rng) = b*y(D_rng) + a*alpha*op(D)*x
		// y(I_rng) = b*y(I_rng) + a*x
		//
		AbstractLinAlgPack::VectorSpace::vec_mut_ptr_t
			y_D = y->sub_view(D_rng),
			y_I = y->sub_view(I_rng);
		// y(D_rng) = b*y(D_rng) + a*alpha*op(D)*x
		AbstractLinAlgPack::Vp_StMtV( y_D.get(), a*alpha, D, D_trans, x, b );
		// y(I_rng) = b*y(I_rng) + a*x
		AbstractLinAlgPack::Vt_S(y_I.get(),b);
		AbstractLinAlgPack::Vp_StV( y_I.get(), a, x );
	}
	else {
		//
		// y = b*y + a*M'*x
		// =>
		// y = b*y + a*alpha*op(D')*x(D_rng) + a*x(I_rng)
		//
		AbstractLinAlgPack::Vp_StMtV( y, a*alpha, D, trans_not(D_trans), *get_view(x,D_rng), b );
		AbstractLinAlgPack::Vp_StV( y, a, *get_view(x,I_rng) );
	}
}

} // end namespace

namespace ConstrainedOptimizationPack {

// Overridden from MatrixBase

size_type MatrixIdentConcat::rows() const
{
	return this->D_rng().size() + this->I_rng().size();
}

size_type MatrixIdentConcat::cols() const
{
	return this->I_rng().size();
}

size_type MatrixIdentConcat::nz() const
{
	const MatrixWithOp& D = this->D();
	return D.nz() + BLAS_Cpp::cols(D.rows(),D.cols(),D_trans()); // D and I
}

// Overridden from MatrixWithOp

std::ostream& MatrixIdentConcat::output(std::ostream& out) const
{
	const Range1D           D_rng   = this->D_rng();
	const BLAS_Cpp::Transp  D_trans = this->D_trans();
	if( D_rng.lbound() == 1 ) {
		if( D_trans == BLAS_Cpp::no_trans )
			out << "[ alpha*D; I ]";
		else
			out << "[ alpha*D'; I ]";
	}
	else {
		if( D_trans == BLAS_Cpp::no_trans )
			out << "[ I; alpha*D ]";
		else
			out << "[ I; alpha*D' ]";
	}
	out << "\nalpha = " << this->alpha();
	out << "\nD =\n" << D();
	return out;
}

void MatrixIdentConcat::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorWithOp& x, value_type b
	) const
{
	AbstractLinAlgPack::Vp_MtV_assert_compatibility(y,*this,M_trans,x);
	mat_vec( y, a, alpha(), D(), D_trans(), D_rng(), I_rng(), M_trans, x, b );
}

void MatrixIdentConcat::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b
	) const
{
	AbstractLinAlgPack::Vp_MtV_assert_compatibility(y,*this,M_trans,x);
//	mat_vec( y, a, alpha(), D(), D_trans(), D_rng(), I_rng(), M_trans, x, b );
	assert(0); // ToDo: Must implement Vp_StV(...) for rhs SpVectorSlice
}

} // end namespace ConstrainedOptimizationPack
