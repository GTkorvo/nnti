// /////////////////////////////////////////////////////////
// GenPermMatrixSliceOp.cpp

#include <assert.h>

#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

void SparseLinAlgPack::V_StMtV(
	  SpVector* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const VectorSlice& x )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	using LinAlgPack::MtV_assert_sizes;

	MtV_assert_sizes( P.rows(), P.cols(), P_trans, x.size() );

	y->resize( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ), P.nz() );

	typedef SpVector::element_type ele_t;

	if( P.is_identity() ) {
		for( size_type i = 1; i <= P.nz(); ++i ) {
			const value_type x_j = x(i);
			if( x_j != 0.0 )
				y->add_element( ele_t( i, a * x_j ) );
		}
	}		
	else if( P_trans == no_trans ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				i = itr->row_i(),
				j = itr->col_j();
			const value_type x_j = x(j);
			if( x_j != 0.0 )
				y->add_element( ele_t( i, a * x_j ) );
		}
	}
	else {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				j = itr->row_i(),
				i = itr->col_j();
			const value_type x_j = x(j);
			if( x_j != 0.0 )
				y->add_element( ele_t( i, a * x_j ) );
		}
	}
	if( P.ordered_by() == GPMSIP::BY_ROW_AND_COL
		|| ( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
		||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL )	)
	{
		y->assume_sorted(true);
	}
}

void SparseLinAlgPack::V_StMtV(
	  SpVector* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const SpVectorSlice& x )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	using LinAlgPack::MtV_assert_sizes;
	MtV_assert_sizes( P.rows(), P.cols(), P_trans, x.size() );

	y->resize( BLAS_Cpp::rows( P.rows(), P.cols(), P_trans ), P.nz() );

	typedef SpVector::element_type ele_t;
	const SpVectorSlice::element_type *ele_ptr;

	if( P.is_identity() ) {
		SparseLinAlgPack::add_elements( y, x(1,P.nz()) );
		SparseLinAlgPack::Vt_S( &(*y)(), a );
	}		
	else if( x.is_sorted() ) {
		if( P_trans == no_trans ) {
			for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
				const size_type
					i = itr->row_i(),
					j = itr->col_j();
				if( ele_ptr = x.lookup_element(j) )
					y->add_element( ele_t( i, a * ele_ptr->value() ) );
			}
		}
		else {
			for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
				const size_type
					j = itr->row_i(),
					i = itr->col_j();
				if( ele_ptr = x.lookup_element(j) )
					y->add_element( ele_t( i, a * ele_ptr->value() ) );
			}
		}
	}
	else {
		assert(0);	// ToDo: Implement the other cases!
	}
	if(	 P.ordered_by() == GPMSIP::BY_ROW_AND_COL
		|| 	( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
		||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL )	)
	{
		y->assume_sorted(true);
	}
}

void SparseLinAlgPack::Vp_StMtV(
	  SpVector* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const VectorSlice& x )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	using LinAlgPack::Vp_MtV_assert_sizes;

	Vp_MtV_assert_sizes( y->size(), P.rows(), P.cols(), P_trans, x.size() );

	typedef SpVector::element_type ele_t;

	const bool was_sorted = y->is_sorted();

	if( P.is_identity() ) {
		for( size_type i = 1; i <= P.nz(); ++i ) {
			const value_type x_j = x(i);
			if( x_j != 0.0 )
				y->add_element( ele_t( i, a * x_j ) );
		}
	}		
	else if( P_trans == no_trans ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				i = itr->row_i(),
				j = itr->col_j();
			const value_type x_j = x(j);
			if( x_j != 0.0 )
				y->add_element( ele_t( i, a * x_j ) );
		}
	}
	else {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				j = itr->row_i(),
				i = itr->col_j();
			const value_type x_j = x(j);
			if( x_j != 0.0 )
				y->add_element( ele_t( i, a * x_j ) );
		}
	}
	if( was_sorted &&
		( P.ordered_by() == GPMSIP::BY_ROW_AND_COL
		  || ( P_trans == no_trans	&& P.ordered_by() == GPMSIP::BY_ROW )
		  ||	( P_trans == trans		&& P.ordered_by() == GPMSIP::BY_COL ) )	)
	{
		y->assume_sorted(true);
	}
}

void SparseLinAlgPack::Vp_StMtV(
	  VectorSlice* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const VectorSlice& x, value_type b )
{
	using LinAlgPack::Vt_S;
	using LinAlgPack::Vp_MtV_assert_sizes;
	Vp_MtV_assert_sizes( y->size(), P.rows(), P.cols(), P_trans, x.size() );
	// y = b*y
	if( b == 0.0 )
		*y = 0.0;
	else
		Vt_S(y,b);	
	// y += a*op(P)*x
	if( P.is_identity() ) {
		if( b == 0.0 )
			*y = 0.0;
		else
			LinAlgPack::Vt_S( y, b );
		LinAlgPack::Vp_StV( &(*y)(1,P.nz()), a, x(1,P.nz()) );
	}		
	else if( P_trans == BLAS_Cpp::no_trans ) {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				i = itr->row_i(),
				j = itr->col_j();
			(*y)(i) += a * x(j);
		}
	}
	else {
		for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
			const size_type
				j = itr->row_i(),
				i = itr->col_j();
			(*y)(i) += a * x(j);
		}
	}
}

void SparseLinAlgPack::Vp_StMtV(
	  VectorSlice* y, value_type a, const GenPermMatrixSlice& P
	, BLAS_Cpp::Transp P_trans, const SpVectorSlice& x, value_type b )
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	using LinAlgPack::Vt_S;
	using LinAlgPack::Vp_MtV_assert_sizes;
	
	Vp_MtV_assert_sizes( y->size(), P.rows(), P.cols(), P_trans, x.size() );
	// y = b*y
	if( b == 0.0 )
		*y = 0.0;
	else
		Vt_S(y,b);
	// y += a*op(P)*x
	if( P.is_identity() ) {
		if( b == 0.0 )
			*y = 0.0;
		else
			LinAlgPack::Vt_S( y, b );
		SparseLinAlgPack::Vp_StV( &(*y)(1,P.nz()), a, x(1,P.nz()) );
	}		
	else if( x.is_sorted() ) {
		const SpVectorSlice::difference_type x_off = x.offset();
		if( P_trans == no_trans && P.ordered_by() == GPMSIP::BY_COL ) {
			assert(0);	// ToDo: implement this!
		}
		else if( P_trans == trans && P.ordered_by() == GPMSIP::BY_ROW ) {
			GenPermMatrixSlice::const_iterator
				P_itr = P.begin(),
				P_end = P.end();
			SpVectorSlice::const_iterator
				x_itr = x.begin(),
				x_end = x.end();
			while( P_itr != P_end && x_itr != x_end ) {
				const size_type
					j = P_itr->row_i(),
					i = P_itr->col_j();
				if( j < x_itr->indice() + x_off ) {
					++P_itr;
					continue;
				}
				else if( j > x_itr->indice() + x_off ) {
					++x_itr;
					continue;
				}
				else {	// they are equal
					(*y)(i) += a * x_itr->value();
					++P_itr;
					++x_itr;
				}
			}
		}
	}
	else {
		// Since things do not match up we will have to create a
		// temporary dense copy of x to operate on.
		assert(0);	// ToDo: Implement this!
	}
}

namespace {

SparseLinAlgPack::GenPermMatrixSliceIteratorPack::EOrderedBy
ordered_by(
	SparseLinAlgPack::GenPermMatrixSliceIteratorPack::EOrderedBy P_ordered_by
	, BLAS_Cpp::Transp P_trans
	)
{
	using BLAS_Cpp::no_trans;
	namespace GPMSIP = SparseLinAlgPack::GenPermMatrixSliceIteratorPack;
	GPMSIP::EOrderedBy
		opP_ordered_by;
	switch( P_ordered_by ) {
	    case GPMSIP::BY_ROW_AND_COL:
			opP_ordered_by = GPMSIP::BY_ROW_AND_COL;
			break;
	    case GPMSIP::BY_ROW:
			opP_ordered_by = P_trans == no_trans ? GPMSIP::BY_ROW : GPMSIP::BY_COL;
			break;
	    case GPMSIP::BY_COL:
			opP_ordered_by = P_trans == no_trans ? GPMSIP::BY_COL : GPMSIP::BY_COL;
			break;
	    case GPMSIP::UNORDERED:
			opP_ordered_by = GPMSIP::UNORDERED;
			break;
  	    default:
			assert(0); // Should never happen
	}
	return opP_ordered_by;
}

} // end namespace

void SparseLinAlgPack::intersection(
	const GenPermMatrixSlice     &P1
	,BLAS_Cpp::Transp            P1_trans
	,const GenPermMatrixSlice    &P2
	,BLAS_Cpp::Transp            P2_trans
	,size_type                   *Q_nz
	,const size_type             Q_max_nz
	,size_type                   Q_row_i[]
	,size_type                   Q_col_j[]
	,GenPermMatrixSlice          *Q
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	namespace GPMSIP = GenPermMatrixSliceIteratorPack;
	//
	// Q = op(P1)*op(P2)
	//
	// There are several different possibilities for how to compute this
	// intersection.
	//
	LinAlgPack::MtM_assert_sizes(
		P1.rows(), P1.cols() , P1_trans, P2.rows(), P2.cols() , P2_trans );
	//
	const size_type
		opP1_rows = rows(P1.rows(),P1.cols(),P1_trans),
		opP1_cols = cols(P1.rows(),P1.cols(),P1_trans),
		opP2_rows = rows(P2.rows(),P2.cols(),P2_trans),
		opP2_cols = cols(P2.rows(),P2.cols(),P2_trans);
	GPMSIP::EOrderedBy
		opP1_ordered_by = ordered_by(P1.ordered_by(),P1_trans),
		opP2_ordered_by = ordered_by(P2.ordered_by(),P2_trans);
	//
	*Q_nz = 0;
	// Either is zero?
	if( !P1.nz() || !P2.nz() ) {
		*Q_nz = 0;
		if(Q)
			Q->initialize(opP1_rows,opP2_cols,GenPermMatrixSlice::ZERO_MATRIX);
		return;
	}
	// Both are identity?
	if( P1.is_identity() && P2.is_identity() ) {
		assert(0); // ToDo: think about this but it is tricky!
		return;
	}
	// One is identity?
	if( P1.is_identity() || P2.is_identity() ) {
		assert(0); // ToDo: Implement this but it is a little tricking?
		return;
	}
	//
	// Both are not identity or zero
	//
	if( ( opP1_ordered_by == GPMSIP::BY_COL || opP1_ordered_by == GPMSIP::BY_ROW_AND_COL )
		&& ( opP2_ordered_by == GPMSIP::BY_ROW || opP2_ordered_by == GPMSIP::BY_ROW_AND_COL ) )
	{
		// This is great!  Both of the matrices are sorted and we don't need any temparary storage
		if( Q_max_nz ) {
			assert(0); // ToDo: Implement initializing Q_row_i, Q_col_j
		}
		else {
			GenPermMatrixSlice::const_iterator
				P1_itr = P1.begin(), // Should not throw exception!
				P1_end = P1.end(),
				P2_itr = P2.begin(), // Should not throw exception!
				P2_end = P2.end();
			while( P1_itr != P1_end && P2_itr != P2_end ) {
				const size_type
					opP1_col_j = cols(P1_itr->row_i(),P1_itr->col_j(),P1_trans),
					opP2_row_i = rows(P2_itr->row_i(),P2_itr->col_j(),P2_trans);
				if( opP1_col_j < opP2_row_i ) {
					++P1_itr;
					continue;
				}
				if( opP1_col_j > opP2_row_i ) {
					++P2_itr;
					continue;
				}
				++(*Q_nz);
				++P1_itr;
				++P2_itr;
			}
		}
	}
	else {
		//
		// We will load the row indices of op(P1) or the column indices op(P1)
		// indexed by column or row indices (whichever is smaller)
		// into a temporary sorted buffer and then loop through the nonzeros of the other.
		//
		// First let's get reorder P1 and P2 so that we put the rows of P2 into a buffer
		//
		const GenPermMatrixSlice
			&oP1      = opP1_cols > opP2_rows ? P1 : P2,
			&oP2      = opP1_cols > opP2_rows ? P2 : P1;
		const BLAS_Cpp::Transp
			oP1_trans = opP1_cols > opP2_rows ? P1_trans : trans_not(P1_trans),
			oP2_trans = opP1_cols > opP2_rows ? P2_trans : trans_not(P2_trans);
		// Load the column indices of op(op(P2)) into a lookup array
		typedef std::vector<size_type> oP2_col_j_lookup_t;      // Todo: use tempoarary workspace
		oP2_col_j_lookup_t oP2_col_j_lookup(rows(oP2.rows(),oP2.rows(),oP2_trans));
		std::fill( oP2_col_j_lookup.begin(), oP2_col_j_lookup.end(), 0 );
		{
			GenPermMatrixSlice::const_iterator
				itr     = oP2.begin(), // Should not throw exception!
				itr_end = oP2.end();
			while( itr != itr_end ) {
				oP2_col_j_lookup[rows(itr->row_i(),itr->col_j(),oP2_trans)]
					= cols(itr->row_i(),itr->col_j(),oP2_trans);
			}
		}
		// Loop through the columns of op(op(P1)) and look for matches
		if( Q_max_nz ) {
			assert(0); // ToDo: Implement initializing Q_row_i, Q_col_j
		}
		else {
			GenPermMatrixSlice::const_iterator
				itr     = oP1.begin(), // Should not throw exception!
				itr_end = oP1.end();
			while( itr != itr_end ) {
				const size_type
					oP2_col_j = oP2_col_j_lookup[cols(itr->row_i(),itr->col_j(),oP1_trans)];
				if(oP2_col_j)
					++(*Q_nz);
			}
		}

	}
	// Setup Q
	assert(Q == NULL); // ToDo: Implement initializing Q
}
