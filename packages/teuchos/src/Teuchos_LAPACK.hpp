// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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


// Kris
// 06.11.03 -- Format cleanup
// 06.13.03 -- Initial templatization (Tpetra_LAPACK.cpp is no longer needed)
// 06.17.03 -- Changed LAPACK wrapper calls from XXXXX_F77() to xxxxx_()
//          -- Added warning messages to default function calls
//          -- Added LAPY2 and GEES by request
// 06.18.03 -- Renamed GEES parameters for default implementation
//          -- Changed LAPACK wrapper calls back to XXXXX_F77() from xxxxx_() (Oops!)
//          -- Streamlined character macro stuff
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_LAPACK_HPP_
#define _TEUCHOS_LAPACK_HPP_

/*! \file Teuchos_LAPACK.hpp
    \brief Templated interface class to LAPACK routines.
*/
/** \example LAPACK/cxx_main.cpp
    This is an example of how to use the Teuchos::LAPACK class.
*/

/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
#if defined (INTEL_CXML)
        unsigned int one=1;
#endif
*/

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, 1
#else
#define CHAR_MACRO(char_var) &char_var
#endif

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_LAPACK_wrappers.hpp"

/*! \class Teuchos::LAPACK
    \brief The Templated LAPACK Wrapper Class.

    The Teuchos::LAPACK class is a wrapper that encapsulates LAPACK
    (Linear Algebra Package).  LAPACK provides portable, high-
    performance implementations of linear, eigen, SVD, etc solvers.

    The standard LAPACK interface is Fortran-specific.  Unfortunately, the
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Teuchos::LAPACK class provides C++ wrappers for the LAPACK
    kernels in order to insulate the rest of Teuchos from the details of C++ to Fortran
    translation.  A Teuchos::LAPACK object is essentially nothing, but allows access to 
    the LAPACK wrapper functions.

    Teuchos::LAPACK is a serial interface only.  This is appropriate since the standard
    LAPACK are only specified for serial execution (or shared memory parallel).

    \note 
	<ol>
	    	<li>These templates are specialized to use the Fortran LAPACK routines for
		scalar types \c float and \c double.
    
		<li>If Teuchos is configured with \c --enable-teuchos-complex then these templates
		are specialized for scalar types \c complex<float> and \c complex<double> also.

		<li>A short description is given for each method.  For more detailed documentation, see the
		LAPACK website (\c http://www.netlib.org/lapack/ ).
	</ol>
*/

namespace Teuchos
{
  template<typename OrdinalType, typename ScalarType>
  class LAPACK
  {    
  public:
    //@{ \name Constructors/Destructors.

    //! Default Constructor.
    inline LAPACK(void) {};

    //! Copy Constructor.
    inline LAPACK(const LAPACK<OrdinalType, ScalarType>& LAPACK) {};

    //! Destructor.
    inline virtual ~LAPACK(void) {};
    //@}

    //@{ \name Symmetric Positive Definite Linear System Routines.

    //! Computes Cholesky factorization of a real symmetric positive definite matrix \c A.    
    void POTRF(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix factored by POTRF and the \c nrhs solutions are returned in \c B.
    void POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes the inverse of a real symmetric positive definite matrix \c A using the Cholesky factorization \c A from POTRF.
    void POTRI(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const;

    //! Estimates the reciprocal of the condition number (1-norm) of a real symmetric positive definite matrix \c A using the Cholesky factorization from POTRF.

    void POCON(const char UPLO, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Computes the solution to a real system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix and the \c nrhs solutions are returned in \c B. 
    void POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes row and column scalings intended to equilibrate a symmetric positive definite matrix \c A and reduce its condition number (w.r.t. 2-norm).
    void POEQU(const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* S, ScalarType* scond, ScalarType* amax, OrdinalType* info) const;

    //! Improves the computed solution to a system of linear equations when the coefficient matrix is symmetric positive definite, and provides error bounds and backward error estimates for the solution.
    void PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Uses the Cholesky factorization to compute the solution to a real system of linear equations \c A*X=B, where \c A is symmetric positive definite.  System can be equilibrated by POEQU and iteratively refined by PORFS, if requested.
    void POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, char EQUED, ScalarType* S, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    //@}

    //@{ \name General Linear System Routines.

    //! Solves an over/underdetermined real \c m by \c n linear system \c A using QR or LQ factorization of A.
    void GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Computes an LU factorization of a general \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GETRF(const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general \c n by \c n matrix \c A using the LU factorization computed by GETRF.
    void GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes the inverse of a matrix \c A using the LU factorization computed by GETRF.
    void GETRI(const OrdinalType n, ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Estimates the reciprocal of the condition number of a general real matrix \c A, in either the 1-norm or the infinity-norm, using the LU factorization computed by GETRF.
    void GECON(const char NORM, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Computes the solution to a real system of linear equations \c A*X=B, where \c A is factored through GETRF and the \c nrhs solutions are computed through GETRS.
    void GESV(const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes row and column scalings intended to equilibrate an \c m by \c n matrix \c A and reduce its condition number.
    void GEEQU(const OrdinalType m, const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* R, ScalarType* C, ScalarType* rowcond, ScalarType* colcond, ScalarType* amax, OrdinalType* info) const;

    //! Improves the computed solution to a system of linear equations and provides error bounds and backward error estimates for the solution.  Use after GETRF/GETRS.
    void GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Uses the LU factorization to compute the solution to a real system of linear equations \c A*X=B, returning error bounds on the solution and a condition estimate.
    void GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, ScalarType* R, ScalarType* C, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Reduces a real general matrix \c A to upper Hessenberg form by orthogonal similarity transformations.
    void GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;
    //@}

    //@{ \name Hessenberg Eigenproblem Routines
    //! Computes the eigenvalues of a real upper Hessenberg matrix \c H and, optionally, the matrices \c T and \c Z from the Schur decomposition, where T is an upper quasi-triangular matrix and Z contains the Schur vectors. 
    void HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* H, const OrdinalType ldh, ScalarType* WR, ScalarType* WI, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;
    
    //! Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues, the real Schur form \c T, and, optionally, the matrix of Schur vectors \c Z.
    void GEES(const char JOBVS, const char SORT, OrdinalType* SELECT, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* WR, ScalarType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const;    

    //! Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
    void GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;
    //@}

    //@{ \name Orthogonal matrix routines
    //! Generates a real orthogonal matrix \c Q which is the product of \c ihi-ilo elementary reflectors of order \c n, as returned by GEHRD.  On return \c Q is stored in \c A.
    void ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Overwrites the general real \c m by \c n matrix \c C with the product of \c C and \c Q, which is a product of \c ihi-ilo elementary reflectors, as returned by GEHRD.
    void ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;
    //@}

    //@{ \name Triangular Matrix Routines
    //! Computes some or all of the right and/or left eigenvectors of a real upper quasi-triangular matrix \c T.
    void TREVC(const char SIDE, const char HOWMNY, OrdinalType* SELECT, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, OrdinalType* info) const;

    //! Reorders the real Schur factorization of a real matrix via orthogonal similarity transformations so that the diagonal block of \c T with row index \c ifst is moved to row \c ilst.
    void TREXC(const char COMPQ, const OrdinalType n, ScalarType* T, const OrdinalType ldt, ScalarType* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, ScalarType* WORK, OrdinalType* info) const;
    //@}

    //@{ \name Random number generators
    //! Returns a random number from a uniform or normal distribution.
    ScalarType LARND( const OrdinalType idist, OrdinalType* seed ) const;

    //! Returns a vector of random numbers from a chosen distribution.
    void LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, ScalarType* v ) const;    
    //@}

    //@{ \name Machine Characteristics Routines.
    //! Determines machine parameters for floating point characteristics.
    ScalarType LAMCH(const char CMACH) const;
    //@}

    //@{ \name Miscellaneous Utilities.
    //! Computes x^2 + y^2 safely, to avoid overflow.
    ScalarType LAPY2(const ScalarType x, const ScalarType y) const;
    //@}
  };

  // END GENERAL TEMPLATE DECLARATION //

  // BEGIN GENERAL TEMPLATE IMPLEMENTATION //


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRF(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POTRF() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POTRS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRI(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POTRI() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POCON(const char UPLO, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POCON() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POSV() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POEQU(const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* S, ScalarType* scond, ScalarType* amax, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POEQU() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::PORFS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, char EQUED, ScalarType* S, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::POSVX() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GELS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GETRF(const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GETRF() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GETRS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GETRI(const OrdinalType n, ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GETRI() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GECON(const char NORM, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GECON() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GESV(const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GESV() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GEEQU(const OrdinalType m, const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* R, ScalarType* C, ScalarType* rowcond, ScalarType* colcond, ScalarType* amax, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GEEQU() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GERFS() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, ScalarType* R, ScalarType* C, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GESVX() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GEHRD() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* H, const OrdinalType ldh, ScalarType* WR, ScalarType* WI, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::HSEQR() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEES(const char JOBVS, const char SORT, OrdinalType* SELECT, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* WR, ScalarType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const    
  {
    std::cout << "Warning: default LAPACK::GEES() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::GEEV() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::ORGHR() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::ORMHR() not yet implemented" << std::endl;
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREVC(const char SIDE, const char HOWMNY, OrdinalType* SELECT, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::TREVC() not yet implemented" << std::endl;
//     if(HOWMNY=='S')
//       {
// 	*info = -3; // We do not support 'S' since it requires a logical array
//       }
//     else
//       {

//       }
  }
  
  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREXC(const char COMPQ, const OrdinalType n, ScalarType* T, const OrdinalType ldt, ScalarType* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, ScalarType* WORK, OrdinalType* info) const
  {
    std::cout << "Warning: default LAPACK::TREXC() not yet implemented" << std::endl;
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAMCH(const char CMACH) const
  {
    std::cout << "Warning: default LAPACK::LAMCH() not yet implemented" << std::endl;
    return 0.0;
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAPY2(const ScalarType x, const ScalarType y) const
  {
    std::cout << "Warning: default LAPACK::LAPY2() not yet implemented" << std::endl;
    return 0.0;
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LARND( const OrdinalType idist, OrdinalType* seed ) const
  {
    std::cout << "Warning: default LAPACK::LARND() not yet implemented" << std::endl;
    return 0.0;
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, ScalarType* v ) const    
  {
    std::cout << "Warning: default LAPACK::LARNV() not yet implemented" << std::endl;
    return 0.0;
  }

  // END GENERAL TEMPLATE IMPLEMENTATION //

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  // BEGIN FLOAT PARTIAL SPECIALIZATION DECLARATION //

  template<typename OrdinalType>
  class LAPACK<OrdinalType, float>
  {    
  public:
    inline LAPACK(void) {};
    inline LAPACK(const LAPACK<OrdinalType, float>& LAPACK) {};
    inline virtual ~LAPACK(void) {};

    // Symmetric positive definite linear system routines
    void POTRF(const char UPLO, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* info) const;
    void POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const float* A, const OrdinalType lda, float* B, const OrdinalType ldb, OrdinalType* info) const;
    void POTRI(const char UPLO, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* info) const;
    void POCON(const char UPLO, const OrdinalType n, const float* A, const OrdinalType lda, const float anorm, float* rcond, float* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* B, const OrdinalType ldb, OrdinalType* info) const;
    void POEQU(const OrdinalType n, const float* A, const OrdinalType lda, float* S, float* scond, float* amax, OrdinalType* info) const;
    void PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, const float* AF, const OrdinalType ldaf, const float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* AF, const OrdinalType ldaf, char EQUED, float* S, float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* rcond, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const; 

    // General Linear System Routines
    void GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* B, const OrdinalType ldb, float* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void GETRF(const OrdinalType m, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const;
    void GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const float* A, const OrdinalType lda, const OrdinalType* IPIV, float* B, const OrdinalType ldb, OrdinalType* info) const;
    void GETRI(const OrdinalType n, float* A, const OrdinalType lda, const OrdinalType* IPIV, float* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void GECON(const char NORM, const OrdinalType n, const float* A, const OrdinalType lda, const float anorm, float* rcond, float* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void GESV(const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, OrdinalType* IPIV, float* B, const OrdinalType ldb, OrdinalType* info) const;
    void GEEQU(const OrdinalType m, const OrdinalType n, const float* A, const OrdinalType lda, float* R, float* C, float* rowcond, float* colcond, float* amax, OrdinalType* info) const;
    void GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const float* A, const OrdinalType lda, const float* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, float* R, float* C, float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* rcond, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, float* A, const OrdinalType lda, float* TAU, float* WORK, const OrdinalType lwork, OrdinalType* info) const;

    // Hessenberg eigenvalue routines.
    void HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, float* H, const OrdinalType ldh, float* WR, float* WI, float* Z, const OrdinalType ldz, float* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void GEES(const char JOBVS, const char SORT, OrdinalType* SELECT, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* sdim, float* WR, float* WI, float* VS, const OrdinalType ldvs, float* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const;    
    void GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, float* A, const OrdinalType lda, float* WR, float* WI, float* VL, const OrdinalType ldvl, float* VR, const OrdinalType ldvr, float* WORK, const OrdinalType lwork, OrdinalType* info) const;

    // Orthogonal matrix routines.
    void ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, float* A, const OrdinalType lda, const float* TAU, float* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const float* A, const OrdinalType lda, const float* TAU, float* C, const OrdinalType ldc, float* WORK, const OrdinalType lwork, OrdinalType* info) const;

    // Triangular matrix routines.
    void TREVC(const char SIDE, const char HOWMNY, OrdinalType* SELECT, const OrdinalType n, const float* T, const OrdinalType ldt, float* VL, const OrdinalType ldvl, float* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, float* WORK, OrdinalType* info) const;
    void TREXC(const char COMPQ, const OrdinalType n, float* T, const OrdinalType ldt, float* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, float* WORK, OrdinalType* info) const;

    // Random number generators
    float LARND( const OrdinalType idist, OrdinalType* seed ) const;
    void LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, float* v ) const;    

    // Machine characteristics.
    float LAMCH(const char CMACH) const;

    // Miscellaneous routines.
    float LAPY2(const float x, const float y) const;
  };

  // END FLOAT PARTIAL SPECIALIZATION DECLARATION //

  // BEGIN FLOAT PARTIAL SPECIALIZATION IMPLEMENTATION //

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRF(const char UPLO, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* info) const
  {
    SPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const float* A, const OrdinalType lda, float* B, const OrdinalType ldb, OrdinalType* info) const
  {
    SPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POTRI(const char UPLO, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* info) const
  {
    SPOTRI_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POCON(const char UPLO, const OrdinalType n, const float* A, const OrdinalType lda, const float anorm, float* rcond, float* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    SPOCON_F77(CHAR_MACRO(UPLO), &n, A, &lda, &anorm, rcond, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* B, const OrdinalType ldb, OrdinalType* info) const
  {
    SPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POEQU(const OrdinalType n, const float* A, const OrdinalType lda, float* S, float* scond, float* amax, OrdinalType* info) const
  {
    SPOEQU_F77(&n, A, &lda, S, scond, amax, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, const float* AF, const OrdinalType ldaf, const float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    SPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* AF, const OrdinalType ldaf, char EQUED, float* S, float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* rcond, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    SPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, CHAR_MACRO(EQUED), S, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* B, const OrdinalType ldb, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GETRF(const OrdinalType m, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
  {
    SGETRF_F77(&m, &n, A, &lda, IPIV, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const float* A, const OrdinalType lda, const OrdinalType* IPIV, float* B, const OrdinalType ldb, OrdinalType* info) const
  {
    SGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, X, &ldx, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GETRI(const OrdinalType n, float* A, const OrdinalType lda, const OrdinalType* IPIV, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SGETRI_F77(&n, A, &lda, IPIV, WORK, lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GECON(const char NORM, const OrdinalType n, const float* A, const OrdinalType lda, const float anorm, float* rcond, float* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    SGECON_F77(CHAR_MACRO(NORM), &n, A, &lda, &anorm, rcond, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GESV(const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, OrdinalType* IPIV, float* B, const OrdinalType ldb, OrdinalType* info) const
  {
    SGESV_F77(&n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GEEQU(const OrdinalType m, const OrdinalType n, const float* A, const OrdinalType lda, float* R, float* C, float* rowcond, float* colcond, float* amax, OrdinalType* info) const
  {
    SGEEQU_F77(&m, &n, A, &lda, R, C, rowcond, colcond, amax, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const float* A, const OrdinalType lda, const float* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    SGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, float* A, const OrdinalType lda, float* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, float* R, float* C, float* B, const OrdinalType ldb, float* X, const OrdinalType ldx, float* rcond, float* FERR, float* BERR, float* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    SGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, CHAR_MACRO(EQUED), R, C, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,float>::GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, float* A, const OrdinalType lda, float* TAU, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SGEHRD_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, float* H, const OrdinalType ldh, float* WR, float* WI, float* Z, const OrdinalType ldz, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, H, &ldh, WR, WI, Z, &ldz, WORK, &lwork, info);
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEES(const char JOBVS, const char SORT, OrdinalType* SELECT, const OrdinalType n, float* A, const OrdinalType lda, OrdinalType* sdim, float* WR, float* WI, float* VS, const OrdinalType ldvs, float* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const    
  {
    SGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), SELECT, &n, A, &lda, sdim, WR, WI, VS, &ldvs, WORK, &lwork, BWORK, info);
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, float* A, const OrdinalType lda, float* WR, float* WI, float* VL, const OrdinalType ldvl, float* VR, const OrdinalType ldvr, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(SORT), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, float* A, const OrdinalType lda, const float* TAU, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SORGHR_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const float* A, const OrdinalType lda, const float* TAU, float* C, const OrdinalType ldc, float* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    SORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &ilo, &ihi, A, &lda, TAU, C, &ldc, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::TREVC(const char SIDE, const char HOWMNY, OrdinalType* SELECT, const OrdinalType n, const float* T, const OrdinalType ldt, float* VL, const OrdinalType ldvl, float* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, float* WORK, OrdinalType* info) const
  {
    if(HOWMNY=='S')
      {
	*info = -3; // We do not support 'S' since it requires a logical array (yuck!)
      }
    else
      {
	STREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, info);
      }
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::TREXC(const char COMPQ, const OrdinalType n, float* T, const OrdinalType ldt, float* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, float* WORK, OrdinalType* info) const
  {
    STREXC_F77(CHAR_MACRO(COMPQ), &n, T, &ldt, Q, &ldq, &ifst, &ilst, WORK, info);
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LARND( const OrdinalType idist, OrdinalType* seed ) const
  {
    return(SLARND_F77(&idist, seed));
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, float>::LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, float* v ) const
  {
    SLARNV_F77(&idist, seed, &n, v);
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LAMCH(const char CMACH) const
  {
    return(SLAMCH_F77(CHAR_MACRO(CMACH)));
  }

  template<typename OrdinalType>
  float LAPACK<OrdinalType, float>::LAPY2(const float x, const float y) const
  {
    return SLAPY2_F77(&x, &y);
  }

  // END FLOAT PARTIAL SPECIALIZATION IMPLEMENTATION //

  // BEGIN DOUBLE PARTIAL SPECIALIZATION DECLARATION //

  template<typename OrdinalType>
  class LAPACK<OrdinalType, double>
  {    
  public:
    inline LAPACK(void) {};
    inline LAPACK(const LAPACK<OrdinalType, double>& LAPACK) {};
    inline virtual ~LAPACK(void) {};

    // Symmetric positive definite linear system routines
    void POTRF(const char UPLO, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* info) const;
    void POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const double* A, const OrdinalType lda, double* B, const OrdinalType ldb, OrdinalType* info) const;
    void POTRI(const char UPLO, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* info) const;
    void POCON(const char UPLO, const OrdinalType n, const double* A, const OrdinalType lda, const double anorm, double* rcond, double* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* B, const OrdinalType ldb, OrdinalType* info) const;
    void POEQU(const OrdinalType n, const double* A, const OrdinalType lda, double* S, double* scond, double* amax, OrdinalType* info) const;
    void PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, const double* AF, const OrdinalType ldaf, const double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* AF, const OrdinalType ldaf, char EQUED, double* S, double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* rcond, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const; 

    // General linear system routines
    void GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* B, const OrdinalType ldb, double* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void GETRF(const OrdinalType m, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const;
    void GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const double* A, const OrdinalType lda, const OrdinalType* IPIV, double* B, const OrdinalType ldb, OrdinalType* info) const;
    void GETRI(const OrdinalType n, double* A, const OrdinalType lda, const OrdinalType* IPIV, double* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void GECON(const char NORM, const OrdinalType n, const double* A, const OrdinalType lda, const double anorm, double* rcond, double* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void GESV(const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, OrdinalType* IPIV, double* B, const OrdinalType ldb, OrdinalType* info) const;
    void GEEQU(const OrdinalType m, const OrdinalType n, const double* A, const OrdinalType lda, double* R, double* C, double* rowcond, double* colcond, double* amax, OrdinalType* info) const;
    void GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const double* A, const OrdinalType lda, const double* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, double* R, double* C, double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* rcond, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    void GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, double* A, const OrdinalType lda, double* TAU, double* WORK, const OrdinalType lwork, OrdinalType* info) const;

    // Hessenberg eigenproblem routines.
    void HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, double* H, const OrdinalType ldh, double* WR, double* WI, double* Z, const OrdinalType ldz, double* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void GEES(const char JOBVS, const char SORT, OrdinalType* SELECT, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* sdim, double* WR, double* WI, double* VS, const OrdinalType ldvs, double* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const;    
    void GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, double* A, const OrdinalType lda, double* WR, double* WI, double* VL, const OrdinalType ldvl, double* VR, const OrdinalType ldvr, double* WORK, const OrdinalType lwork, OrdinalType* info) const;

    // Orthogonal matrix routines.
    void ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, double* A, const OrdinalType lda, const double* TAU, double* WORK, const OrdinalType lwork, OrdinalType* info) const;
    void ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const double* A, const OrdinalType lda, const double* TAU, double* C, const OrdinalType ldc, double* WORK, const OrdinalType lwork, OrdinalType* info) const;

    // Triangular matrix routines.
    void TREVC(const char SIDE, const char HOWMNY, OrdinalType* SELECT, const OrdinalType n, const double* T, const OrdinalType ldt, double* VL, const OrdinalType ldvl, double* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, double* WORK, OrdinalType* info) const;
    void TREXC(const char COMPQ, const OrdinalType n, double* T, const OrdinalType ldt, double* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, double* WORK, OrdinalType* info) const;

    // Random number generators
    double LARND( const OrdinalType idist, OrdinalType* seed ) const;
    void LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, double* v ) const;    

    // Machine characteristic routines.
    double LAMCH(const char CMACH) const;

    // Miscellaneous routines.
    double LAPY2(const double x, const double y) const;
  };

  // END DOUBLE PARTIAL SPECIALIZATION DECLARATION //

  // BEGIN DOUBLE PARTIAL SPECIALIZATION IMPLEMENTATION //


  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRF(const char UPLO, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* info) const
  {
    DPOTRF_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const double* A, const OrdinalType lda, double* B, const OrdinalType ldb, OrdinalType* info) const
  {
    DPOTRS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POTRI(const char UPLO, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* info) const
  {
    DPOTRI_F77(CHAR_MACRO(UPLO), &n, A, &lda, info);
  }
  
  template<typename OrdinalType>
    void LAPACK<OrdinalType, double>::POCON(const char UPLO, const OrdinalType n, const double* A, const OrdinalType lda, const double anorm, double* rcond, double* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    DPOCON_F77(CHAR_MACRO(UPLO), &n, A, &lda, &anorm, rcond, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* B, const OrdinalType ldb, OrdinalType* info) const
  {
    DPOSV_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::POEQU(const OrdinalType n, const double* A, const OrdinalType lda, double* S, double* scond, double* amax, OrdinalType* info) const
  {
    DPOEQU_F77(&n, A, &lda, S, scond, amax, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, const double* AF, const OrdinalType ldaf, const double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    DPORFS_F77(CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
    void LAPACK<OrdinalType, double>::POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* AF, const OrdinalType ldaf, char EQUED, double* S, double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* rcond, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const 
  {
    DPOSVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(UPLO), &n, &nrhs, A, &lda, AF, &ldaf, CHAR_MACRO(EQUED), S, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* B, const OrdinalType ldb, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DGELS_F77(CHAR_MACRO(TRANS), &m, &n, &nrhs, A, &lda, B, &ldb, WORK, &lwork, &info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GETRF(const OrdinalType m, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
  {
    DGETRF_F77(&m, &n, A, &lda, IPIV, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const double* A, const OrdinalType lda, const OrdinalType* IPIV, double* B, const OrdinalType ldb, OrdinalType* info) const
  {
    DGETRS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GETRI(const OrdinalType n, double* A, const OrdinalType lda, const OrdinalType* IPIV, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DGETRI_F77(&n, A, &lda, IPIV, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GECON(const char NORM, const OrdinalType n, const double* A, const OrdinalType lda, const double anorm, double* rcond, double* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    DGECON_F77(CHAR_MACRO(NORM), &n, A, &lda, &anorm, rcond, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GESV(const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, OrdinalType* IPIV, double* B, const OrdinalType ldb, OrdinalType* info) const
  {
    DGESV_F77(&n, &nrhs, A, &lda, IPIV, B, &ldb, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GEEQU(const OrdinalType m, const OrdinalType n, const double* A, const OrdinalType lda, double* R, double* C, double* rowcond, double* colcond, double* amax, OrdinalType* info) const
  {
    DGEEQU_F77(&m, &n, A, &lda, R, C, rowcond, colcond, amax, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const double* A, const OrdinalType lda, const double* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    DGERFS_F77(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, B, &ldb, X, &ldx, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType,double>::GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, double* A, const OrdinalType lda, double* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, double* R, double* C, double* B, const OrdinalType ldb, double* X, const OrdinalType ldx, double* rcond, double* FERR, double* BERR, double* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    DGESVX_F77(CHAR_MACRO(FACT), CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, AF, &ldaf, IPIV, CHAR_MACRO(EQUED), R, C, B, &ldb, X, &ldx, rcond, FERR, BERR, WORK, IWORK, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, double* A, const OrdinalType lda, double* TAU, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DGEHRD_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, double* H, const OrdinalType ldh, double* WR, double* WI, double* Z, const OrdinalType ldz, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DHSEQR_F77(CHAR_MACRO(JOB), CHAR_MACRO(COMPZ), &n, &ilo, &ihi, H, &ldh, WR, WI, Z, &ldz, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEES(const char JOBVS, const char SORT, OrdinalType* SELECT, const OrdinalType n, double* A, const OrdinalType lda, OrdinalType* sdim, double* WR, double* WI, double* VS, const OrdinalType ldvs, double* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const    
  {
    DGEES_F77(CHAR_MACRO(JOBVS), CHAR_MACRO(SORT), SELECT, &n, A, &lda, sdim, WR, WI, VS, &ldvs, WORK, &lwork, BWORK, info);
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, double* A, const OrdinalType lda, double* WR, double* WI, double* VL, const OrdinalType ldvl, double* VR, const OrdinalType ldvr, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DGEEV_F77(CHAR_MACRO(JOBVL), CHAR_MACRO(SORT), &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, WORK, &lwork, info);
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, double* A, const OrdinalType lda, const double* TAU, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DORGHR_F77(&n, &ilo, &ihi, A, &lda, TAU, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const double* A, const OrdinalType lda, const double* TAU, double* C, const OrdinalType ldc, double* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    DORMHR_F77(CHAR_MACRO(SIDE), CHAR_MACRO(TRANS), &m, &n, &ilo, &ihi, A, &lda, TAU, C, &ldc, WORK, &lwork, info);
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::TREVC(const char SIDE, const char HOWMNY, OrdinalType* SELECT, const OrdinalType n, const double* T, const OrdinalType ldt, double* VL, const OrdinalType ldvl, double* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, double* WORK, OrdinalType* info) const
  {
    if(HOWMNY=='S')
      {
	*info = -3; // We do not support 'S' since it requires a logical array (yuck!)
      }
    else
      {
	DTREVC_F77(CHAR_MACRO(SIDE), CHAR_MACRO(HOWMNY), SELECT, &n, T, &ldt, VL, &ldvl, VR, &ldvr, &mm, m, WORK, info);
      }
  }
  
  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::TREXC(const char COMPQ, const OrdinalType n, double* T, const OrdinalType ldt, double* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, double* WORK, OrdinalType* info) const
  {
    DTREXC_F77(CHAR_MACRO(COMPQ), &n, T, &ldt, Q, &ldq, &ifst, &ilst, WORK, info);
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LARND( const OrdinalType idist, OrdinalType* seed ) const
  {
    return(DLARND_F77(&idist, seed));
  }

  template<typename OrdinalType>
  void LAPACK<OrdinalType, double>::LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, double* v ) const
  {
    DLARNV_F77(&idist, seed, &n, v);
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LAMCH(const char CMACH) const
  {
    return(DLAMCH_F77(CHAR_MACRO(CMACH)));
  }

  template<typename OrdinalType>
  double LAPACK<OrdinalType, double>::LAPY2(const double x, const double y) const
  {
    return DLAPY2_F77(&x, &y);
  }

  // END DOUBLE PARTIAL SPECIALIZATION IMPLEMENTATION //
#endif

} // namespace Teuchos

#endif // _TEUCHOS_LAPACK_HPP_
