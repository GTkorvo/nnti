// //////////////////////////////////////////////////////////////////////////////////
// MatrixVarReductImplicit.cpp
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

#include "ConstrainedOptimizationPack/include/MatrixVarReductImplicit.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
//#include "SparseLinAlgPack/include/dense_Vp_StPtMtV.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/GenPermMatrixSlice.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "WorkspacePack.h"

namespace {

/*

//
// Implicit matrix-vector multiplication:
//
// y = b*y + a*op(inv(C)*N)*x
//
template<class V>
void imp_Vp_StMtV_implicit(
	LinAlgPack::VectorSlice                                             *y
	,LinAlgPack::value_type                                             a
	,const ConstrainedOptimizationPack::DecompositionSystemVarReduct    &decomp_sys
	,BLAS_Cpp::Transp                                                   D_trans
	,const V                                                            &x
	,LinAlgPack::value_type                                             b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	const LinAlgPack::size_type
		r   = decomp_sys.C().rows(),
		dof = decomp_sys.N().cols();

	if( D_trans == no_trans ) {
		//
		// y += a * inv(C) * ( N * x )
		//
		wsp::Workspace<LinAlgPack::value_type> t1_ws(wss,r), t2_ws(wss,r);
		LinAlgPack::VectorSlice t1(&t1_ws[0],t1_ws.size()), t2(&t2_ws[0],t2_ws.size());
		// t1 = N*x
		LinAlgOpPack::V_MtV( &t1, decomp_sys.N(), no_trans, x );
		// t2 = inv(C) * t1
		decomp_sys.solve_C( t1, no_trans, &t2 );
		// y = b*y
		if(b==0.0)       *y = 0.0;
		else if(b!=1.0)  LinAlgPack::Vt_S(y,b);
		// y += a*t2
		LinAlgPack::Vp_StV( y, -a, t2 );
	}
	else {
		//
		// y = b*y + a * N' * ( inv(C') * x )
		//
		wsp::Workspace<LinAlgPack::value_type> t1_ws(wss,r);
		LinAlgPack::VectorSlice t1(&t1_ws[0],t1_ws.size());
		// t1 = inv(C')*x
		decomp_sys.solve_C( x, trans, &t1 );
		// y = b*y + a*N'*t1
	    SparseLinAlgPack::Vp_StMtV( y, a,  decomp_sys.N(), trans, t1, b );
	}
}

//
// Generate a row of inv(C)*N if not already computed.
//
void imp_compute_InvCtN_row(
	LinAlgPack::size_type                                                 r
	,const ConstrainedOptimizationPack::DecompositionSystemVarReduct      &decomp_sys
	,LinAlgPack::size_type                                                j
	,LinAlgPack::VectorSlice                                              *e_j  // Set to all zeros on input and output!
	,ConstrainedOptimizationPack::MatrixVarReductImplicit::InvCtN_rows_t  *InvCtN_rows
	)
{
	typedef  LinAlgPack::value_type  value_type;
	using LinAlgPack::VectorSlice;
	if( (*InvCtN_rows)[j-1] == NULL ) {
		// Generate row j of inv(C)*N
		value_type *vec = (*InvCtN_rows)[j-1] = new value_type[r]; // ToDo: We may want to allocate more vectors at once!
		VectorSlice row_j(vec,r);
		// row_j = N'*inv(C')*e_j
		(*e_j)(j) = 1.0;
		imp_Vp_StMtV_implicit( &row_j, 1.0, decomp_sys, BLAS_Cpp::trans, *e_j, 0.0 );
		(*e_j)(j) = 0.0;
	}
}

//
// Perform the matrix-vector multiplication:
// 
// y = b*y -a * op(P) * [inv(C) * N] * x
//
// by generating rows [inv(C)*N](j,:) for each nonzero entry op(P)(i,j).
//
template<class V>
void imp_Vp_StPtMtV_by_row(
	LinAlgPack::VectorSlice                                               *y
	,LinAlgPack::value_type                                               a
	,const ConstrainedOptimizationPack::GenPermMatrixSlice                &P
	,BLAS_Cpp::Transp                                                     P_trans
	,const ConstrainedOptimizationPack::DecompositionSystemVarReduct      &decomp_sys
	,const V                                                              &x
	,LinAlgPack::value_type                                               b
	,ConstrainedOptimizationPack::MatrixVarReductImplicit::InvCtN_rows_t *InvCtN_rows
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using LinAlgPack::dot;
	using LinAlgPack::VectorSlice;
	using SparseLinAlgPack::dot;
	using SparseLinAlgPack::Vp_StMtV;
	using SparseLinAlgPack::GenPermMatrixSlice;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	
	const LinAlgPack::size_type
		D_rows = decomp_sys.C().rows(),
		D_cols = decomp_sys.N().cols();
	// y = b*y
	if(b==0.0)       *y = 0.0;
	else if(b!=1.0)  LinAlgPack::Vt_S(y,b);
	// Compute t = N'*inv(C')*e(j) then y(i) += -a*t'*x where op(P)(i,j) = 1.0
	wsp::Workspace<LinAlgPack::value_type>   e_j_ws(wss,D_rows);
	VectorSlice                              e_j(&e_j_ws[0],e_j_ws.size());
	e_j = 0.0;
	for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
		const LinAlgPack::size_type
			i = P_trans == no_trans ? itr->row_i() : itr->col_j(),
			j = P_trans == no_trans ? itr->col_j() : itr->row_i();
		// t = op(M') * e(j)
		imp_compute_InvCtN_row(D_rows,decomp_sys,j,&e_j,InvCtN_rows);
		VectorSlice t((*InvCtN_rows)[j-1],D_cols);
		// y(i) += -a*t'*x
		(*y)(i) += (-a) * dot(t,x);
	}
}

*/

} // end namespace

namespace ConstrainedOptimizationPack {

void MatrixVarReductImplicit::initialize(
	const mat_nonsing_ptr_t          &C
	,const mat_ptr_t                 &N
	,const mat_ptr_t                 &D_direct
	)
{
/*
	// Validate the inputs
	if( !decomp_sys )
		throw std::invalid_argument(
			"MatrixVarReductImplicit::initialize(...): Error, "
			"decomp_sys must not be NULL" );
	if( D_dense && (D_dense->rows() != decomp_sys->C().rows() || D_dense->cols() != decomp_sys->N().cols() ) )
		throw std::invalid_argument(
			"MatrixVarReductImplicit::initialize(...): Error, "
			"*D_dense does not match the size of the decomposition in *decomp_sys" );
	// Set the members
	decomp_sys_ = decomp_sys;
	if( D_dense )
		D_dense_.bind( const_cast<GenMatrixSlice&>(*D_dense) );
	use_dense_mat_vec_ = true; // ToDo: We should use a timer to determine if this is faster than sparse or not!
	release_resource_ptr_ = release_resource_ptr;
	if(InvCtN_rows_.size()) { // Free previously allocated vectors
		for( InvCtN_rows_t::iterator itr = InvCtN_rows_.begin(); itr != InvCtN_rows_.end(); ++itr ) {
			if( *itr )
				delete [] *itr; // ToDo: We may want to allocate vectors in larger chuncks
			*itr = (value_type*)NULL;
		}
	}
*/
	assert(0); // ToDo: Initialize the above!
}

void MatrixVarReductImplicit::set_uninitialized()
{
	C_        = NULL;
	N_        = NULL;
	D_direct_ = NULL;
}

// Overridden from MatrixBase

size_type MatrixVarReductImplicit::rows() const
{
	return C_.get() ? C_->rows() : 0;
}

size_type MatrixVarReductImplicit::cols() const
{
	return N_.get() ? N_->cols() : 0;
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixVarReductImplicit::space_cols() const
{
	assert_initialized();
	return C_->space_cols();
}

const VectorSpace& MatrixVarReductImplicit::space_rows() const
{
	assert_initialized();
	return N_->space_rows();
}

MatrixWithOp& MatrixVarReductImplicit::operator=(const MatrixWithOp& M)
{
	assert_initialized();
	assert(0); // ToDo: Finish!
	return *this;
}

std::ostream& MatrixVarReductImplicit::output(std::ostream& o)
{
	return MatrixWithOp::output(o); // ToDo: Specialize!
}

void MatrixVarReductImplicit::Vp_StMtV(
	VectorWithOpMutable* y, value_type a
	,BLAS_Cpp::Transp D_trans
	,const VectorWithOp& x, value_type b
	) const
{
	assert_initialized();
/*
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),rows(),cols(),D_trans,x.size());
	if( use_dense_mat_vec_ && D_dense_.rows() > 0 ) {
		LinAlgOpPack::Vp_StMtV( y, a, D_dense_, D_trans, x, b );
	}
	else {
		imp_Vp_StMtV_implicit( y, a, *decomp_sys_, D_trans, x, b );
	}
*/
	assert(0); // ToDo: Update above!
}

void MatrixVarReductImplicit::Vp_StMtV(
	VectorWithOpMutable* y, value_type a
	,BLAS_Cpp::Transp D_trans
	,const SpVectorSlice& x, value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using LinAlgPack::Vt_S;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	assert_initialized();
/*
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),D_rows,D_cols,D_trans,x.size());
	if( use_dense_mat_vec_ && D_dense_.rows() > 0 ) {
		LinAlgOpPack::Vp_StMtV( y, a, D_dense_, D_trans, x, b );
	}
	else {
		if( x.nz() == x.size() ) {  // This is B.S.  Should use MatrixWithOpFactorized object for C!
			VectorSlice dx = SparseLinAlgPack::dense_view(x);
			imp_Vp_StMtV_implicit( y, -a, *decomp_sys_, D_trans, dx, b );
		}
		else if( D_trans == BLAS_Cpp::trans && x.nz() < D_cols ) {
			//
			// y = b*y + (-a)*[N'*inv(C')]*x
			//
			// We can do something crafty here.  We can generate columns of N'*inv(C')
			// and then perform y += -a*[N'*inv(C')](:,j)*x(j) for nonzero x(j)
			//
			wsp::Workspace<LinAlgPack::value_type>   e_j_ws(wss,D_rows);
			VectorSlice                              e_j(&e_j_ws[0],e_j_ws.size());
			e_j = 0.0;
			// y = b*y
			if(b==0.0)       *y = 0.0;
			else if(b!=1.0)  Vt_S(y,b);
			// y += -a*[N'*inv(C')](:,j)*x(j), for j <: { j | x(j) != 0.0 }
			const SpVectorSlice::difference_type o = x.offset();
			for( SpVectorSlice::const_iterator itr = x.begin(); itr != x.end(); ++itr ) {
				const size_type j = itr->indice() + o;
				imp_compute_InvCtN_row(D_rows,*decomp_sys_,j,&e_j,&InvCtN_rows_);
				LinAlgPack::Vp_StV( y, -a * itr->value(), VectorSlice(InvCtN_rows_[j-1],D_cols) );
			}
		}
		else {   // This is B.S.  Should use MatrixWithOpFactorized object for C!
			Vector dx;
			LinAlgOpPack::assign( &dx, x );
			imp_Vp_StMtV_implicit( y, -a, *decomp_sys_, D_trans, dx(), b );
		}
	}
*/
	assert(0); // ToDo: Update above!
}

void MatrixVarReductImplicit::Vp_StPtMtV(
	VectorWithOpMutable* y, value_type a
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,BLAS_Cpp::Transp D_trans
	,const VectorWithOp& x, value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	assert_initialized();
/*
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans,opD_rows);
	LinAlgPack::Vp_MtV_assert_sizes(cols(P.rows(),P.cols(),P_trans),D_rows,D_cols,D_trans,x.size());
	if( D_dense_.rows() > 0 ) {
		SparseLinAlgPack::dense_Vp_StPtMtV(y,a,P,P_trans,D_dense_,D_trans,x,b);
	}
	else if( P.nz() > D_cols || D_trans == BLAS_Cpp::trans ) {
		// Just use the default implementation
		MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b);
	}
	else {
		imp_Vp_StPtMtV_by_row(y,a,P,P_trans,*decomp_sys_,x,b,&InvCtN_rows_);
	}
*/
	assert(0); // ToDo: Update above!
}

void MatrixVarReductImplicit::Vp_StPtMtV(
	VectorWithOpMutable* y, value_type a
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,BLAS_Cpp::Transp D_trans
	,const SpVectorSlice& x, value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	assert_initialized();
/*
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	LinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans,opD_rows);
	LinAlgPack::Vp_MtV_assert_sizes(cols(P.rows(),P.cols(),P_trans),D_rows,D_cols,D_trans,x.size());
	if( D_dense_.rows() > 0 ) {
		SparseLinAlgPack::dense_Vp_StPtMtV(y,a,P,P_trans,D_dense_,D_trans,x,b);
	}
	else if( P.nz() > D_cols || D_trans == BLAS_Cpp::trans ) {
		// Just use the default implementation
		MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b);
	}
	else {
		imp_Vp_StPtMtV_by_row(y,a,P,P_trans,*decomp_sys_,x,b,&InvCtN_rows_);
	}
*/
	assert(0); // ToDo: Update above!
}

// Private member functions

void MatrixVarReductImplicit::assert_initialized() const
{
	THROW_EXCEPTION(
		C_.get() == NULL, std::logic_error
		,"MatrixVarReductImplicit::assert_initialized(): Error, "
		"initialize(...) has not been called yet!" );
}

}	// end namespace ConstrainedOptimizationPack 
