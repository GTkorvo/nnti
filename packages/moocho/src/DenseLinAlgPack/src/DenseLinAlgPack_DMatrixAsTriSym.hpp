// //////////////////////////////////////////////////////////////////////////////////
// GenMatrixAsTriSym.h

#ifndef GEN_MATRIX_AS_TRI_SYM_H
#define GEN_MATRIX_AS_TRI_SYM_H

#include "GenMatrixClass.h"

namespace LinAlgPack {

/** @name	Packaging arguments for GenMatrixSlices treated as triangular
  *			and symmetric matrices in BLAS-like linear algebra operations.
  *
  */

//@{

// /////////////////////////////////////////////////////////////////////////////////////
///
/** Aggregate information for a triangular matrix (element-wise) stored in a GenMatrix.
  *
  * This is the type to be used as lhs and rhs arguments in element-wise
  * linear algebra operations like assignment and binary arithmetic..
  */
class tri_ele_gms {
public:
	///
	tri_ele_gms(const GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
		: gms_(const_cast<GenMatrixSlice&>(gms)), uplo_(uplo)
	{
		#ifdef LINALGPACK_CHECK_RHS_SIZES
			assert_gms_square(gms);
		#endif
	}

	///
	size_type rows() const {
		return gms_.rows();
	}

	///
	size_type cols() const {
		return gms_.cols();
	}

	///
	GenMatrixSlice& gms() {
		return gms_;
	}

	///
	const GenMatrixSlice& gms() const {
		return gms_;
	}

	///
	BLAS_Cpp::Uplo	uplo() const {
		return uplo_;
	}

private:	
	GenMatrixSlice	gms_;
	BLAS_Cpp::Uplo	uplo_;
};	// end class tri_ele_gms

inline
/// Return a triangular element-wise matrix
tri_ele_gms tri_ele(GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
	return tri_ele_gms(gms, uplo);
}

inline
///
const tri_ele_gms tri_ele(const GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
	return tri_ele_gms(gms, uplo);
}

// ////////////////////////////////////////////////////////////////////////////////
///
/** Aggregate information for a triangular matrix (structure dependent) stored in a GenMatrix.
  *
  * This is the type to be used as a rhs argument in linear algebra operations
  * that are structure specific like the BLAS operations.
  */
class tri_gms {
public:
	///
	tri_gms(const GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
		: gms_(const_cast<GenMatrixSlice&>(gms)), uplo_(uplo), diag_(diag)
	{
		#ifdef LINALGPACK_CHECK_RHS_SIZES
			assert_gms_square(gms);
		#endif
	}

	///
	size_type rows() const {
		return gms_.rows();
	}

	///
	size_type cols() const {
		return gms_.cols();
	}

	///
	GenMatrixSlice& gms() {
		return gms_;
	}

	///
	const GenMatrixSlice& gms() const {
		return gms_;
	}

	///
	BLAS_Cpp::Uplo	uplo() const {
		return uplo_;
	}

	///
	BLAS_Cpp::Diag	diag() const {
		return diag_;
	}

private:	
	GenMatrixSlice	gms_;
	BLAS_Cpp::Uplo	uplo_;
	BLAS_Cpp::Diag	diag_;
};	// end class tri_gms

inline
/// Return a triangular matrix
tri_gms tri(GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
{
	return tri_gms(gms, uplo, diag);
}

inline
///
const tri_gms tri(const GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
{
	return tri_gms(gms, uplo, diag);
}

// /////////////////////////////////////////////////////////////////////////////////////////
///
/** Aggregate information for a symmetric matrix stored in a GenMatrix.
  *
  * This is the type to be used as both lhs and rhs arguments in linear algebra operations.
  */
class sym_gms {
public:
	///
	sym_gms(const GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
		: gms_(const_cast<GenMatrixSlice&>(gms)), uplo_(uplo)
	{
		#ifdef LINALGPACK_CHECK_RHS_SIZES
			assert_gms_square(gms);
		#endif
	}

	///
	size_type rows() const {
		return gms_.rows();
	}

	///
	size_type cols() const {
		return gms_.cols();
	}

	///
	GenMatrixSlice& gms() {
		return gms_;
	}

	///
	const GenMatrixSlice& gms() const {
		return gms_;
	}

	///
	BLAS_Cpp::Uplo	uplo() const {
		return uplo_;
	}

private:	
	GenMatrixSlice	gms_;
	BLAS_Cpp::Uplo	uplo_;
};	// end class sym_gms

inline
/// Return a symmetric matrix
sym_gms sym(GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
	return sym_gms(gms, uplo);
}

inline
///
const sym_gms sym(const GenMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
	return sym_gms(gms, uplo);
}

//@}

} // end namespace LinAlgPack

#endif	// GEN_MATRIX_AS_TRI_SYM_H
