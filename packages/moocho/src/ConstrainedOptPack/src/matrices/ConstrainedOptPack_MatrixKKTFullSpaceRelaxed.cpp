// //////////////////////////////////////
// MatrixKKTFullSpaceRelaxed.cpp
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
//
// ToDo: 6/2/00: Implement the relaxed version in the future.

#include <assert.h>

#include <sstream>
#include <typeinfo>

#include "ConstrainedOptimizationPack/src/MatrixKKTFullSpaceRelaxed.hpp"
#include "AbstractLinAlgPack/src/DirectSparseFortranCompatibleSolver.h"
#include "AbstractLinAlgPack/src/MatrixConvertToSparseFortranCompatible.hpp"
#include "AbstractLinAlgPack/test/TestMatrixConvertToSparseFortranCompatible.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DenseLinAlgPackAssertOp.hpp"

namespace ConstrainedOptimizationPack {

MatrixKKTFullSpaceRelaxed::MatrixKKTFullSpaceRelaxed(
		const direct_solver_ptr_t& direct_solver		)
	:
		  direct_solver_(direct_solver)
		, initialized_(false)
		, n_(0), m_(0)
		, G_(NULL), convG_(0), G_nz_(0)
		, A_(NULL), convA_(0), A_nz_(0)
{}

void MatrixKKTFullSpaceRelaxed::initialize(
	  const MatrixOp& G, const MatrixOp& A
	, std::ostream* out, EPrintMoreOrLess print_what, ERunTests test_what )
{
	// Set some members first
	out_			= out;
	print_what_		= print_what;
	test_what_		= test_what;

	// Validate that the matrices check out and get the conversion
	// interfaces.
	validate_and_set_matrices(G,A);

	use_relaxation_ = false;

	// Factor this matrix.
	//
	// If the structure of G and A looks to be the same we will reuse the
	// previous factorization structure if we have one.

	// ToDo: Finish this

	// Now we are initialized and ready to go.
	initialized_ = true;
}

void MatrixKKTFullSpaceRelaxed::initialize_relaxed(
	  const MatrixOp& G, const MatrixOp& A
	, const DVectorSlice& c, value_type M
	, std::ostream* out, EPrintMoreOrLess print_what, ERunTests test_what )
{
	// ToDo: implement the relaxation in the future!
	assert(0);
}

void MatrixKKTFullSpaceRelaxed::set_uninitialized()
{
	G_ = NULL;
	A_ = NULL;
	initialized_ = false;
}

void MatrixKKTFullSpaceRelaxed::release_memory()
{
	direct_solver().release_memory();
	n_ = m_ = 0;
	G_nz_ = A_nz_ = 0;
	set_uninitialized();
}

// Overridden from Matrix

size_type MatrixKKTFullSpaceRelaxed::rows() const
{
	assert_initialized();
	return n_ + m_ + (use_relaxation_ ? 1 : 0 );
}

size_type MatrixKKTFullSpaceRelaxed::cols() const
{
	assert_initialized();
	return n_ + m_ + (use_relaxation_ ? 1 : 0 );
}

// Overridden from MatrixOp

std::ostream& MatrixKKTFullSpaceRelaxed::output(std::ostream& out) const
{
	assert_initialized();
	// ToDo: Finish me!
	assert(0);
	return out;
}

MatrixOp& MatrixKKTFullSpaceRelaxed::operator=(const MatrixOp& m)
{
	assert_initialized();
	// ToDo: Finish me!
	assert(0);
	return *this;
}

void MatrixKKTFullSpaceRelaxed::Vp_StMtV(
	  DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2, value_type beta) const
{
	using AbstractLinAlgPack::Vp_StMtV;

	assert_initialized();

	DenseLinAlgPack::Vp_MtV_assert_sizes(vs_lhs->size(),rows(),cols(),BLAS_Cpp::no_trans
		,vs_rhs2.size());

	if( use_relaxation_ ) {
		// ToDo: Implement the relaxation in the future
		assert(0);
	}
	else {
		// y = b*y + a*K*x
		//
		// [y1] = b * [y1] + a * [ G  A ] * [x1]    
		// [y2]       [y2]       [ A'   ]   [x2]
		//
		// y1 = b*y1 + a*G*x1 + a*A*x2
		//
		// y2 = b*y2 + a*A'*x1

		DVectorSlice
			y1 = (*vs_lhs)(1,n_),
			y2 = (*vs_lhs)(n_+1,n_+m_);
		const DVectorSlice
			x1 = vs_rhs2(1,n_),
			x2 = vs_rhs2(n_+1,n_+m_);

		// y1 = b*y1 + a*G*x1 + a*A*x2

		// y1 = a*G*x1 + b*y1
		Vp_StMtV( &y1, alpha, *G_, BLAS_Cpp::no_trans, x1, beta );
		// y1 += a*A*x2
		Vp_StMtV( &y1, alpha, *A_, BLAS_Cpp::no_trans, x2 );

		// y2 = a*A'*x1 + b*y2
		Vp_StMtV( &y2, alpha, *A_, BLAS_Cpp::trans, x1, beta );
	}
}

// Overridden from MatrixFactorized

void MatrixKKTFullSpaceRelaxed::V_InvMtV( DVectorSlice* v_lhs, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2) const
{
	assert_initialized();
	// ToDo: Finish me!
	assert(0);
}

// Overridden from MatrixConvertToSparseFortranCompatible

FortranTypes::f_int
MatrixKKTFullSpaceRelaxed::num_nonzeros( EExtractRegion extract_region ) const
{
	assert_matrices_set();

	FortranTypes::f_int
		nz;
	if( use_relaxation_ ) {
		// ToDo: Add support for the relaxation.
		assert(0);
	}
	else {
		// Return the number of nonzeros in a region of :
		//
		// K = [ G  A ]
		//     [ A'   ]
		typedef MatrixConvertToSparseFortranCompatible MCTSFC_t;
		const FortranTypes::f_int
			A_nz = convA_->num_nonzeros(MCTSFC_t::EXTRACT_FULL_MATRIX);
		nz = convG_->num_nonzeros( extract_region )
			+ ( extract_region == MCTSFC_t::EXTRACT_FULL_MATRIX ? 2 * A_nz : A_nz );
	}
	return nz;
}

void MatrixKKTFullSpaceRelaxed::coor_extract_nonzeros(
	  EExtractRegion extract_region
	, const FortranTypes::f_int len_Aval
		, FortranTypes::f_dbl_prec Aval[]
	, const FortranTypes::f_int len_Aij
		, FortranTypes::f_int Arow[]
		, FortranTypes::f_int Acol[]
		, const FortranTypes::f_int row_offset
		, const FortranTypes::f_int col_offset
	 ) const
{
	assert_matrices_set();

	if( len_Aval == 0 && len_Aij == 0 ) {
		if(*out_)
			*out_
				<< "\n*** MatrixKKTFullSpaceRelaxed::coor_extract_nonzeros(...) : "
				<< "Warning, nothing to compute: len_Aval==0 && len_Aij==0\n";
		return;
	}

	// Validate the conversion interfaces if asked to.
	if( test_what_ == RUN_TESTS ) {

		typedef MatrixConvertToSparseFortranCompatible
			MCTSFC_t;
		namespace slaptp = AbstractLinAlgPack::TestingPack;
		using slaptp::TestMatrixConvertToSparseFortranCompatible;

		const value_type
			warning_tol		= 1e-10,	// There should be very little roundoff error
			error_tol		= 1e-6;

		// Test G
		{
			slaptp::ETestSparseFortranFormat
				sparse_format = slaptp::TEST_COOR_FORMAT;

			if(out_)
				*out_
					<< "\n*** Testing conversion interface for submatrix G ...\n";

			bool result = TestMatrixConvertToSparseFortranCompatible(
				 extract_region,sparse_format,*convG_,*G_
				,warning_tol,error_tol,true,out_
				,print_what_==PRINT_MORE );

			if( !result) {
				char omsg[] = "MatrixKKTFullSpaceRelaxed : Error, the conversion "
					"interface for G did not check out\n";
				if(out_)
					*out_
						<< std::endl << omsg;
				throw std::logic_error(omsg);
			}
			else {
				if(out_)
					*out_
						<< "\n*** Conversion interface for G checked out!\n";
			}
		}

		// Test A
		{
			MCTSFC_t::EExtractRegion
				extract_region = MCTSFC_t::EXTRACT_FULL_MATRIX;
			slaptp::ETestSparseFortranFormat
				sparse_format = slaptp::TEST_COOR_FORMAT;

			if(out_)
				*out_
					<< "\n*** Testing conversion interface for submatrix A ...\n";

			bool result = TestMatrixConvertToSparseFortranCompatible(
				 extract_region,sparse_format,*convA_,*A_
				,warning_tol,error_tol,true,out_
				,print_what_==PRINT_MORE );

			if( !result) {
				char omsg[] = "MatrixKKTFullSpaceRelaxed : Error, the conversion "
					"interface for A did not check out\n";
				if(out_)
					*out_
						<< std::endl << omsg;
				throw std::logic_error(omsg);
			}
			else {
				if(out_)
					*out_
						<< "\n*** Conversion interface for A checked out!\n";
			}
		}

	}

	// Extract the nonzero entries
	if( use_relaxation_ ) {
		// ToDo: Add support for the relaxation.
		assert(0);
	}
	else {
		// Extract the nonzero entries in a region of :
		//
		// K = [ G  A ]
		//     [ A'   ]

		typedef MatrixConvertToSparseFortranCompatible MCTSFC_t;

		switch(extract_region) {
			case MCTSFC_t::EXTRACT_FULL_MATRIX:
				assert(0);	// We don't support this yet!
				break;
			case MCTSFC_t::EXTRACT_UPPER_TRIANGULAR:
				assert(0);	// We don't support this yet!
				break;
			case MCTSFC_t::EXTRACT_LOWER_TRIANGULAR:
			{
				// Set elements for
				//
				// K_lower = [ G_lower     ] n
				//           [ A'          ] m
				//             n        m

				// Get the numbers of nonzeros
				const FortranTypes::f_int
					G_lo_nz = convG_->num_nonzeros( MCTSFC_t::EXTRACT_LOWER_TRIANGULAR ),
					A_nz = convA_->num_nonzeros( MCTSFC_t::EXTRACT_FULL_MATRIX);
				assert( (len_Aval == 0 || len_Aval == G_lo_nz + A_nz)
						&& (len_Aij == 0 || len_Aij == G_lo_nz + A_nz) );
				// Set the elements for G
				convG_->coor_extract_nonzeros(
					 MCTSFC_t::EXTRACT_LOWER_TRIANGULAR
					, (len_Aval > 0 ? G_lo_nz : 0), Aval
					, (len_Aij > 0 ? G_lo_nz : 0), Arow, Acol, 0, 0 );
				// Set the elements for A'
				// To do this we will reverse the row and column indices
				// and the row offset will be n (col_offset in argument list)
				convA_->coor_extract_nonzeros(
					 MCTSFC_t::EXTRACT_FULL_MATRIX
					, (len_Aval > 0 ? A_nz : 0), Aval+G_lo_nz
					, (len_Aij > 0 ? A_nz: 0), Acol+G_lo_nz, Arow+G_lo_nz, 0, n_ );
				break;
			}
			default:
				assert(0);
				break;
		}
	}
}

// Private member functions

void MatrixKKTFullSpaceRelaxed::assert_initialized() const
{
	if(!initialized_) {
		throw NotInitializedException("MatrixKKTFullSpaceRelaxed::assert_initialized() : "
			"Error, called a member function without initialize(...) or "
			"initialize_relaxed(...) being called first probably" );
	}
}

void MatrixKKTFullSpaceRelaxed::assert_matrices_set() const
{
	if( !G_ || !convG_ || !A_ || !convA_ ) {
		throw NotInitializedException("MatrixKKTFullSpaceRelaxed::assert_matrices_set() : "
			"Error, the matrices G and A have not been set yet" );
	}
}

void MatrixKKTFullSpaceRelaxed::validate_and_set_matrices(
	  const MatrixOp& G, const MatrixOp& A )
{
	const size_type
		n = G.rows(),
		m = A.cols();

	if( G.cols() != n ) {
		std::ostringstream omsg;
		omsg
			<< "MatrixKKTFullSpaceRelaxed::validate_and_set_matrices(...) : Error, "
			<< "The matrix G with rows = " << n
			<< " is not square with cols = " << G.cols();
		throw std::length_error(omsg.str());
	}
	if( A.rows() != n ) {
		std::ostringstream omsg;
		omsg
			<< "MatrixKKTFullSpaceRelaxed::validate_and_set_matrices(...) : Error, "
			<< "The number of rows in the matrix A with rows = " << A.rows()
			<< ", cols = " << m
			<< " does not match the size of G with rows = cols = " << n;
		throw std::length_error(omsg.str());
	}
	if( !(m < n) ||  n == 0 || m == 0 ) {
		std::ostringstream omsg;
		omsg
			<< "MatrixKKTFullSpaceRelaxed::validate_and_set_matrices(...) : Error, "
			<< "the size of the matrices G nxn and A nxm is not valid with "
			<< "n = " << n << " and m = " << m;
		throw std::length_error(omsg.str());
	}

	const MatrixConvertToSparseFortranCompatible
		*convG = dynamic_cast<const MatrixConvertToSparseFortranCompatible*>(&G);
	if( ! convG ) {
		std::ostringstream omsg;
		omsg
			<< "MatrixKKTFullSpaceRelaxed::validate_and_set_matrices(...) : Error, "
			<< "The matrix G with concrete type " << typeid(G).name()
			<< " does not support the MatrixConvertToSparseFortranCompatible "
			<< "interface";
		throw InvalidMatrixType(omsg.str());
	}
	const MatrixConvertToSparseFortranCompatible
		*convA = dynamic_cast<const MatrixConvertToSparseFortranCompatible*>(&A);
	if( ! convA ) {
		std::ostringstream omsg;
		omsg
			<< "MatrixKKTFullSpaceRelaxed::validate_and_set_matrices(...) : Error, "
			<< "The matrix A with concrete type " << typeid(A).name()
			<< " does not support the MatrixConvertToSparseFortranCompatible "
			<< "interface";
		throw InvalidMatrixType(omsg.str());
	}

	// The input matrices checkout so set this stuff now
	G_ = &G;
	convG_ = convG;
	G_nz_ = G.nz();
	A_ = &A;
	convA_ = convA;
	A_nz_ = A.nz();
	n_ = n;
	m_ = m;
}


}	// end namespace ConstrainedOptimizationPack
