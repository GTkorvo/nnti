// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_BLAS_HPP
#define SACADO_FAD_BLAS_HPP

#include "Teuchos_BLAS.hpp"
#include "Sacado.hpp"

namespace Sacado {

  namespace Fad {

    template <typename OrdinalType, typename FadType>
    class ArrayTraits {

      typedef typename Sacado::ValueType<FadType>::type ValueType;
      
    public:
      
      ArrayTraits(bool use_dynamic = true,
		  OrdinalType workspace_size = 0);

      ArrayTraits(const ArrayTraits& a);

      ~ArrayTraits();
      
      void unpack(const FadType& a, OrdinalType& n_dot, ValueType& val, 
		  const ValueType*& dot) const;
      
      void unpack(const FadType* a, OrdinalType n, OrdinalType inc,
		  OrdinalType& n_dot, OrdinalType& inc_val, 
		  OrdinalType& inc_dot,
		  const ValueType*& val, const ValueType*& dot) const;
      
      void unpack(const FadType* A, OrdinalType m, OrdinalType n, 
		  OrdinalType lda, OrdinalType& n_dot, 
		  OrdinalType& lda_val, OrdinalType& lda_dot,
		  const ValueType*& val, const ValueType*& dot) const;

      void unpack(FadType& a, OrdinalType& n_dot, OrdinalType& final_n_dot, 
		  ValueType& val, ValueType*& dot) const;
      
      void unpack(FadType* a, OrdinalType n, OrdinalType inc,
		  OrdinalType& n_dot, OrdinalType& final_n_dot, 
		  OrdinalType& inc_val, OrdinalType& inc_dot,
		  ValueType*& val, ValueType*& dot) const;
      
      void unpack(FadType* A, OrdinalType m, OrdinalType n, OrdinalType lda, 
		  OrdinalType& n_dot, OrdinalType& final_n_dot, 
		  OrdinalType& lda_val, OrdinalType& lda_dot,
		  ValueType*& val, ValueType*& dot) const;

      void pack(FadType& a, OrdinalType n_dot, const ValueType& val, 
		const ValueType* dot) const;
      
      void pack(FadType* a, OrdinalType n, OrdinalType inc,
		OrdinalType n_dot, OrdinalType inc_val, OrdinalType inc_dot,
		const ValueType* val, const ValueType* dot) const;
      
      void pack(FadType* A, OrdinalType m, OrdinalType n, 
		OrdinalType lda, OrdinalType n_dot, 
		OrdinalType lda_val, OrdinalType lda_dot,
		const ValueType* val, const ValueType* dot) const;

      void free(const FadType& a, OrdinalType n_dot, 
		const ValueType* dot) const;
      
      void free(const FadType* a, OrdinalType n, OrdinalType n_dot,
		OrdinalType inc_val, OrdinalType inc_dot,
		const ValueType* val, const ValueType* dot) const;

      void free(const FadType* A, OrdinalType m, OrdinalType n, 
		OrdinalType n_dot, OrdinalType lda_val, OrdinalType lda_dot,
		const ValueType* val, const ValueType* dot) const;

      ValueType* allocate_array(OrdinalType size) const;

      void free_array(const ValueType* ptr, OrdinalType size) const;

      bool is_array_contiguous(const FadType* a, OrdinalType n, 
			       OrdinalType n_dot) const;

    protected:

      //! Use dynamic memory allocation
      bool use_dynamic;

      //! Size of static workspace
      OrdinalType workspace_size;

      //! Workspace for holding contiguous values/derivatives
      mutable ValueType *workspace;

      //! Pointer to current free entry in workspace
      mutable ValueType *workspace_pointer;
		
    };

    //! Fad specializations for Teuchos::BLAS wrappers
    template <typename OrdinalType, typename FadType>
    class BLAS : public Teuchos::DefaultBLASImpl<OrdinalType,FadType> {    
      
      typedef typename Teuchos::ScalarTraits<FadType>::magnitudeType MagnitudeType;
      typedef FadType ScalarType;
      typedef typename Sacado::ValueType<FadType>::type ValueType;
      typedef Teuchos::DefaultBLASImpl<OrdinalType,ScalarType> BLASType;
    
    public:
      //! @name Constructor/Destructor.
      //@{ 
    
      //! Default constructor.
      BLAS(bool use_default_impl = true,
	   bool use_dynamic = true, OrdinalType static_workspace_size = 0);

      //! Copy constructor.

      BLAS(const BLAS& x);

      //! Destructor.
      virtual ~BLAS();

      //@}
      
      //! @name Level 1 BLAS Routines.
      //@{ 
      
      //! Computes a Givens plane rotation.
      void ROTG(ScalarType* da, ScalarType* db, MagnitudeType* c, 
		ScalarType* s) const { 
	BLASType::ROTG(da,db,c,s); 
      }
      
      //! Applies a Givens plane rotation.
      void ROT(const OrdinalType n, ScalarType* dx, const OrdinalType incx, 
	       ScalarType* dy, const OrdinalType incy, MagnitudeType* c, 
	       ScalarType* s) const { 
	BLASType::ROT(n,dx,incx,dy,incy,c,s); 
      }
      
      //! Scale the std::vector \c x by the constant \c alpha.
      void SCAL(const OrdinalType n, const ScalarType& alpha, ScalarType* x, 
		const OrdinalType incx) const;

      //! Copy the std::vector \c x to the std::vector \c y.
      void COPY(const OrdinalType n, const ScalarType* x, 
		const OrdinalType incx, ScalarType* y, 
		const OrdinalType incy) const;

      //! Perform the operation: \c y \c <- \c y+alpha*x.
      void AXPY(const OrdinalType n, const ScalarType& alpha, 
		const ScalarType* x, const OrdinalType incx, ScalarType* y, 
		const OrdinalType incy) const;

      //! Perform the operation: \c y \c <- \c y+alpha*x.
      /*!
       * Overload for the case when \c x is constant.
       */
      void AXPY(const OrdinalType n, const ScalarType& alpha, 
		const ValueType* x, const OrdinalType incx, ScalarType* y, 
		const OrdinalType incy) const;

      //! Sum the absolute values of the entries of \c x.
      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
      ASUM(const OrdinalType n, const ScalarType* x, 
	   const OrdinalType incx) const {
	return BLASType::ASUM(n,x,incx);
      }

      //! Form the dot product of the vectors \c x and \c y.
      ScalarType DOT(const OrdinalType n, const ScalarType* x, 
		     const OrdinalType incx, const ScalarType* y, 
		     const OrdinalType incy) const;
      
      //! Form the dot product of the vectors \c x and \c y.
      /*!
       * Overload for the case when \c x is constant.
       */
      ScalarType DOT(const OrdinalType n, const ValueType* x, 
		     const OrdinalType incx, const ScalarType* y, 
		     const OrdinalType incy) const;

      //! Form the dot product of the vectors \c x and \c y.
      /*!
       * Overload for the case when \c y is constant.
       */
      ScalarType DOT(const OrdinalType n, const ScalarType* x, 
		     const OrdinalType incx, const ValueType* y, 
		     const OrdinalType incy) const;

      //! Compute the 2-norm of the std::vector \c x.
      MagnitudeType NRM2(const OrdinalType n, const ScalarType* x, 
			 const OrdinalType incx) const;

      //! Return the index of the element of \c x with the maximum magnitude.
      OrdinalType IAMAX(const OrdinalType n, const ScalarType* x, 
			const OrdinalType incx) const {
	return BLASType::IAMAX(n,x,incx); 
      }

      //@}
      
      //! @name Level 2 BLAS Routines.
      //@{ 
      
      /*! 
       * \brief Performs the matrix-std::vector operation:  
       * \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y 
       * where \c A is a general \c m by \c n matrix.
       */
      /*!
       * If \c alpha or \c beta are constant, they will be automatically
       * promoted to Fad types by the compiler with little computational cost.
       */
      void GEMV(Teuchos::ETransp trans, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType alpha, const ScalarType* A, 
		const OrdinalType lda, const ScalarType* x, 
		const OrdinalType incx, const ScalarType beta, 
		ScalarType* y, const OrdinalType incy) const;

      /*! 
       * \brief Performs the matrix-std::vector operation:  
       * \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y 
       * where \c A is a general \c m by \c n matrix.
       */
      /*!
       * Overload for case when matrix \c A is constant
       */
      void GEMV(Teuchos::ETransp trans, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType alpha, const ValueType* A, 
		const OrdinalType lda, const ScalarType* x, 
		const OrdinalType incx, const ScalarType beta, 
		ScalarType* y, const OrdinalType incy) const;

      /*! 
       * \brief Performs the matrix-std::vector operation:  
       * \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y 
       * where \c A is a general \c m by \c n matrix.
       */
      /*!
       * Overload for case when vector \c x is constant
       */
      void GEMV(Teuchos::ETransp trans, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType alpha, const ScalarType* A, 
		const OrdinalType lda, const ValueType* x, 
		const OrdinalType incx, const ScalarType beta, 
		ScalarType* y, const OrdinalType incy) const;

      /*! 
       * \brief Performs the matrix-std::vector operation:  
       * \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y 
       * where \c A is a general \c m by \c n matrix.
       */
      /*!
       * Overload for case when \C A and \c x are constant
       */
      void GEMV(Teuchos::ETransp trans, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType alpha, const ValueType* A, 
		const OrdinalType lda, const ValueType* x, 
		const OrdinalType incx, const ScalarType beta, 
		ScalarType* y, const OrdinalType incy) const;

      /*!
       * \brief Performs the matrix-std::vector operation:  
       * \c x \c <- \c A*x or \c x \c <- \c A'*x where \c A is a unit/non-unit 
       * \c n by \c n upper/lower triangular matrix.
       */
      void TRMV(Teuchos::EUplo uplo, Teuchos::ETransp trans, 
		Teuchos::EDiag diag, const OrdinalType n, 
		const ScalarType* A, const OrdinalType lda, ScalarType* x, 
		const OrdinalType incx) const;

      /*!
       * \brief Performs the matrix-std::vector operation:  
       * \c x \c <- \c A*x or \c x \c <- \c A'*x where \c A is a unit/non-unit 
       * \c n by \c n upper/lower triangular matrix.
       */
      /*!
       * Overload for case when matrix \c A is constant
       */
      void TRMV(Teuchos::EUplo uplo, Teuchos::ETransp trans, 
		Teuchos::EDiag diag, const OrdinalType n, 
		const ValueType* A, const OrdinalType lda, ScalarType* x, 
		const OrdinalType incx) const;

      //! Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A. 
      void GER(const OrdinalType m, const OrdinalType n, 
	       const ScalarType& alpha, 
	       const ScalarType* x, const OrdinalType incx, 
	       const ScalarType* y, const OrdinalType incy, 
	       ScalarType* A, const OrdinalType lda) const;

      //! Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A. 
      /*!
       * Overload for case when vector \c x is constant
       */
      void GER(const OrdinalType m, const OrdinalType n, 
	       const ScalarType& alpha, 
	       const ValueType* x, const OrdinalType incx, 
	       const ScalarType* y, const OrdinalType incy, 
	       ScalarType* A, const OrdinalType lda) const;

      //! Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A. 
      /*!
       * Overload for case when vector \c y is constant
       */
      void GER(const OrdinalType m, const OrdinalType n, 
	       const ScalarType& alpha, 
	       const ScalarType* x, const OrdinalType incx, 
	       const ValueType* y, const OrdinalType incy, 
	       ScalarType* A, const OrdinalType lda) const;

      //! Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A. 
      /*!
       * Overload for case when vectors \c x and \c y are constant
       */
      void GER(const OrdinalType m, const OrdinalType n, 
	       const ScalarType& alpha, 
	       const ValueType* x, const OrdinalType incx, 
	       const ValueType* y, const OrdinalType incy, 
	       ScalarType* A, const OrdinalType lda) const;
      
      //@}
      
      //! @name Level 3 BLAS Routines. 
      //@{ 
      
      /*! 
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*op(B)+beta*C where \c op(A) is either \c A 
       * or \c A', \c op(B) is either \c B or \c B', and C is an \c m by \c k 
       * matrix.
       */
      /*!
       * If \c alpha or \c beta are constant, they will be automatically
       * promoted to Fad types by the compiler with little computational cost.
       */
      void GEMM(Teuchos::ETransp transa, Teuchos::ETransp transb, 
		const OrdinalType m, 
		const OrdinalType n, const OrdinalType k, 
		const ScalarType& alpha, 
		const ScalarType* A, const OrdinalType lda, 
		const ScalarType* B, 
		const OrdinalType ldb, const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;

      /*! 
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*op(B)+beta*C where \c op(A) is either \c A 
       * or \c A', \c op(B) is either \c B or \c B', and C is an \c m by \c k 
       * matrix.
       */
      /*!
       * Overload for case when matrix \c A is constant
       */
      void GEMM(Teuchos::ETransp transa, Teuchos::ETransp transb, 
		const OrdinalType m, 
		const OrdinalType n, const OrdinalType k, 
		const ScalarType& alpha, 
		const ValueType* A, const OrdinalType lda, 
		const ScalarType* B, 
		const OrdinalType ldb, const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;

      /*! 
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*op(B)+beta*C where \c op(A) is either \c A 
       * or \c A', \c op(B) is either \c B or \c B', and C is an \c m by \c k 
       * matrix.
       */
      /*!
       * Overload for case when matrix \c B is constant
       */
      void GEMM(Teuchos::ETransp transa, Teuchos::ETransp transb, 
		const OrdinalType m, 
		const OrdinalType n, const OrdinalType k, 
		const ScalarType& alpha, 
		const ScalarType* A, const OrdinalType lda, 
		const ValueType* B, 
		const OrdinalType ldb, const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;

      /*! 
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*op(B)+beta*C where \c op(A) is either \c A 
       * or \c A', \c op(B) is either \c B or \c B', and C is an \c m by \c k 
       * matrix.
       */
      /*!
       * Overload for case when matrices \c A and \c B are constant
       */
      void GEMM(Teuchos::ETransp transa, Teuchos::ETransp transb, 
		const OrdinalType m, 
		const OrdinalType n, const OrdinalType k, 
		const ScalarType& alpha, 
		const ValueType* A, const OrdinalType lda, 
		const ValueType* B, 
		const OrdinalType ldb, const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;
      
      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where 
       * \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a 
       * general matrix.
       */
      /*!
       * If \c alpha or \c beta are constant, they will be automatically
       * promoted to Fad types by the compiler with little computational cost.
       */
      void SYMM(Teuchos::ESide side, Teuchos::EUplo uplo, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType& alpha, const ScalarType* A, 
		const OrdinalType lda, const ScalarType* B, 
		const OrdinalType ldb,
		const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;

      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where 
       * \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a 
       * general matrix.
       */
      /*!
       * Overload for case when matrix \c A is constant
       */
      void SYMM(Teuchos::ESide side, Teuchos::EUplo uplo, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType& alpha, const ValueType* A, 
		const OrdinalType lda, const ScalarType* B, 
		const OrdinalType ldb,
		const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;

      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where 
       * \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a 
       * general matrix.
       */
      /*!
       * Overload for case when matrix \c B is constant
       */
      void SYMM(Teuchos::ESide side, Teuchos::EUplo uplo, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType& alpha, const ScalarType* A, 
		const OrdinalType lda, const ValueType* B, 
		const OrdinalType ldb,
		const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;

      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where 
       * \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a 
       * general matrix.
       */
      /*!
       * Overload for case when matrices \c A and \c B are constant
       */
      void SYMM(Teuchos::ESide side, Teuchos::EUplo uplo, const OrdinalType m, 
		const OrdinalType n, 
		const ScalarType& alpha, const ValueType* A, 
		const OrdinalType lda, const ValueType* B, 
		const OrdinalType ldb,
		const ScalarType& beta, ScalarType* C, 
		const OrdinalType ldc) const;
      
      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*B+beta*C or 
       * \c C \c <- \c alpha*B*op(A)+beta*C where \c op(A) is an unit/non-unit, 
       * upper/lower triangular matrix and \c B is a general matrix.
       */
      /*!
       * If \c alpha or \c beta are constant, they will be automatically
       * promoted to Fad types by the compiler with little computational cost.
       */
      void TRMM(Teuchos::ESide side, Teuchos::EUplo uplo, 
		Teuchos::ETransp transa, Teuchos::EDiag diag, 
		const OrdinalType m, const OrdinalType n, 
		const ScalarType& alpha, 
		const ScalarType* A, const OrdinalType lda, 
		ScalarType* B, const OrdinalType ldb) const;

      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*B+beta*C or 
       * \c C \c <- \c alpha*B*op(A)+beta*C where \c op(A) is an unit/non-unit, 
       * upper/lower triangular matrix and \c B is a general matrix.
       */
      /*!
       * Overload for case when matrix \c A is constant
       */
      void TRMM(Teuchos::ESide side, Teuchos::EUplo uplo, 
		Teuchos::ETransp transa, Teuchos::EDiag diag, 
		const OrdinalType m, const OrdinalType n, 
		const ScalarType& alpha, 
		const ValueType* A, const OrdinalType lda, 
		ScalarType* B, const OrdinalType ldb) const;

      /*! 
       * \brief Solves the matrix equations:  
       * \c op(A)*X=alpha*B or \c X*op(A)=alpha*B where \c X and \c B are \c m 
       * by \c n matrices, \c A is a unit/non-unit, upper/lower triangular 
       * matrix and \c op(A) is \c A or \c A'.  The matrix \c X is overwritten 
       * on \c B.
       */
      /*!
       * If \c alpha or \c beta are constant, they will be automatically
       * promoted to Fad types by the compiler with little computational cost.
       */
      void TRSM(Teuchos::ESide side, Teuchos::EUplo uplo, 
		Teuchos::ETransp transa, Teuchos::EDiag diag, 
		const OrdinalType m, const OrdinalType n, 
		const ScalarType& alpha, 
		const ScalarType* A, const OrdinalType lda, 
		ScalarType* B, const OrdinalType ldb) const;

      /*! 
       * \brief Solves the matrix equations:  
       * \c op(A)*X=alpha*B or \c X*op(A)=alpha*B where \c X and \c B are \c m 
       * by \c n matrices, \c A is a unit/non-unit, upper/lower triangular 
       * matrix and \c op(A) is \c A or \c A'.  The matrix \c X is overwritten 
       * on \c B.
       */
      /*!
       * Overload for case when matrix \c A is constant
       */
      void TRSM(Teuchos::ESide side, Teuchos::EUplo uplo, 
		Teuchos::ETransp transa, Teuchos::EDiag diag, 
		const OrdinalType m, const OrdinalType n, 
		const ScalarType& alpha, 
		const ValueType* A, const OrdinalType lda, 
		ScalarType* B, const OrdinalType ldb) const;

      //@}

    protected:

      //! ArrayTraits for packing/unpacking value/derivative arrays
      ArrayTraits<OrdinalType,FadType> arrayTraits;

      //! BLAS for values
      Teuchos::BLAS<OrdinalType, ValueType> blas;

      //! Use custom or default implementation
      bool use_default_impl;
      
      //! Temporary array for GEMV
      mutable std::vector<ValueType> gemv_Ax;

      //! Temporary array for GEMM
      mutable std::vector<ValueType> gemm_AB;

    protected:

      //! Implementation of DOT
      void Fad_DOT(const OrdinalType n,
		   const ValueType* x, 
		   const OrdinalType incx, 
		   const OrdinalType n_x_dot,  
		   const ValueType* x_dot, 
		   const OrdinalType incx_dot,
		   const ValueType* y, 
		   const OrdinalType incy, 
		   const OrdinalType n_y_dot, 
		   const ValueType* y_dot,
		   const OrdinalType incy_dot,
		   ValueType& z,
		   const OrdinalType n_z_dot,
		   ValueType* zdot) const;

      //! Implementation of GEMV
      void Fad_GEMV(Teuchos::ETransp trans, 
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const ValueType& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const ValueType* alpha_dot,
		    const ValueType* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const ValueType* A_dot,
		    const OrdinalType lda_dot,
		    const ValueType* x, 
		    const OrdinalType incx, 
		    const OrdinalType n_x_dot, 
		    const ValueType* x_dot, 
		    const OrdinalType incx_dot, 
		    const ValueType& beta, 
		    const OrdinalType n_beta_dot, 
		    const ValueType* beta_dot,
		    ValueType* y, 
		    const OrdinalType incy, 
		    const OrdinalType n_y_dot, 
		    ValueType* y_dot,
		    const OrdinalType incy_dot,
		    const OrdinalType n_dot) const;

      //! Implementation of GEMV
      void Fad_GER(const OrdinalType m, 
		   const OrdinalType n, 
		   const ValueType& alpha, 
		   const OrdinalType n_alpha_dot, 
		   const ValueType* alpha_dot,
		   const ValueType* x, 
		   const OrdinalType incx, 
		   const OrdinalType n_x_dot, 
		   const ValueType* x_dot, 
		   const OrdinalType incx_dot, 
		   const ValueType* y, 
		   const OrdinalType incy, 
		   const OrdinalType n_y_dot, 
		   const ValueType* y_dot,
		   const OrdinalType incy_dot,
		   ValueType* A, 
		   const OrdinalType lda, 
		   const OrdinalType n_A_dot, 
		   ValueType* A_dot,
		   const OrdinalType lda_dot,
		   const OrdinalType n_dot) const;

      //! Implementation of GEMM
      void Fad_GEMM(Teuchos::ETransp transa,
		    Teuchos::ETransp transb,
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const OrdinalType k,
		    const ValueType& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const ValueType* alpha_dot,
		    const ValueType* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const ValueType* A_dot,
		    const OrdinalType lda_dot,
		    const ValueType* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    const ValueType* B_dot, 
		    const OrdinalType ldb_dot, 
		    const ValueType& beta, 
		    const OrdinalType n_beta_dot, 
		    const ValueType* beta_dot,
		    ValueType* C, 
		    const OrdinalType ldc, 
		    const OrdinalType n_C_dot, 
		    ValueType* C_dot,
		    const OrdinalType ldc_dot,
		    const OrdinalType n_dot) const;

      //! Implementation of SYMM
      void Fad_SYMM(Teuchos::ESide side, 
		    Teuchos::EUplo uplo,
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const ValueType& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const ValueType* alpha_dot,
		    const ValueType* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const ValueType* A_dot,
		    const OrdinalType lda_dot,
		    const ValueType* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    const ValueType* B_dot, 
		    const OrdinalType ldb_dot, 
		    const ValueType& beta, 
		    const OrdinalType n_beta_dot, 
		    const ValueType* beta_dot,
		    ValueType* C, 
		    const OrdinalType ldc, 
		    const OrdinalType n_C_dot, 
		    ValueType* C_dot,
		    const OrdinalType ldc_dot,
		    const OrdinalType n_dot) const;

      //! Implementation of TRMM
      void Fad_TRMM(Teuchos::ESide side, 
		    Teuchos::EUplo uplo,
		    Teuchos::ETransp transa, 
		    Teuchos::EDiag diag, 
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const ValueType& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const ValueType* alpha_dot,
		    const ValueType* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const ValueType* A_dot,
		    const OrdinalType lda_dot,
		    ValueType* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    ValueType* B_dot, 
		    const OrdinalType ldb_dot, 
		    const OrdinalType n_dot) const;

      //! Implementation of TRMM
      void Fad_TRSM(Teuchos::ESide side, 
		    Teuchos::EUplo uplo,
		    Teuchos::ETransp transa, 
		    Teuchos::EDiag diag, 
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const ValueType& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const ValueType* alpha_dot,
		    const ValueType* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const ValueType* A_dot,
		    const OrdinalType lda_dot,
		    ValueType* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    ValueType* B_dot, 
		    const OrdinalType ldb_dot, 
		    const OrdinalType n_dot) const;

    }; // class FadBLAS

  }  // namespace Fad

} // namespace Sacado

// Here we provide partial specializations for Teuchos::BLAS for each Fad type
#define TEUCHOS_BLAS_FAD_SPEC(FADTYPE)					\
namespace Teuchos {							\
  template <typename OrdinalType, typename ValueT, typename ScalarT>	\
  class BLAS< OrdinalType, FADTYPE<ValueT,ScalarT> > :			\
    public Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,ScalarT> > {	\
  public:								\
    BLAS(bool use_default_impl = true,	bool use_dynamic = true,	\
	 OrdinalType static_workspace_size = 0) :			\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,ScalarT> >(	\
	use_default_impl, use_dynamic,static_workspace_size) {}		\
    BLAS(const BLAS& x) :						\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,ScalarT> >(x) {}	\
    virtual ~BLAS() {}							\
  };									\
}
#define TEUCHOS_BLAS_SFAD_SPEC(FADTYPE)					\
namespace Teuchos {							\
  template <typename OrdinalType, typename ValueT, int Num,		\
	    typename ScalarT>						\
  class BLAS< OrdinalType, FADTYPE<ValueT,Num,ScalarT> > :		\
    public Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,Num,ScalarT> > { \
  public:								\
    BLAS(bool use_default_impl = true,	bool use_dynamic = true,	\
	 OrdinalType static_workspace_size = 0) :			\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,Num,ScalarT> >(	\
	use_default_impl, use_dynamic, static_workspace_size) {}	\
    BLAS(const BLAS& x) :						\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,Num,ScalarT> >(x) {} \
    virtual ~BLAS() {}							\
  };									\
}
TEUCHOS_BLAS_FAD_SPEC(Sacado::Fad::DFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::Fad::SFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::Fad::SLFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::Fad::DMFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::Fad::DVFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::ELRFad::DFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::ELRFad::SFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::ELRFad::SLFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::CacheFad::DFad)

#undef TEUCHOS_BLAS_FAD_SPEC
#undef TEUCHOS_BLAS_SFAD_SPEC

#include "Sacado_Fad_BLASImp.hpp"

#endif // SACADO_FAD_BLAS_HPP 
