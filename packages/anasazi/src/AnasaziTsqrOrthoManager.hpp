// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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

/// \file AnasaziTsqrOrthoManager.hpp
/// \brief Orthogonalization manager based on Tall Skinny QR (TSQR)

#ifndef __AnasaziTsqrOrthoManager_hpp
#define __AnasaziTsqrOrthoManager_hpp

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziTsqrAdaptor.hpp"
//#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOrthoManager.hpp"
#include "Teuchos_LAPACK.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Anasazi {

  /// \class TsqrOrthoError
  /// \brief TsqrOrthoManager error
  class TsqrOrthoError : public OrthoError
  {
  public: 
    TsqrOrthoError (const std::string& what_arg) : 
      OrthoError (what_arg) {}
  };

  /// \class TsqrOrthoFault
  /// \brief Orthogonalization fault
  ///
  /// Stewart (SISC 2008) gives a Block Gram-Schmidt (BGS) with
  /// reorthogonalization algorithm.  An "orthogonalization fault" is
  /// what happens when the second BGS pass does not succeed (which is
  /// possible in BGS, but not possible in (non-block) Gram-Schmidt if
  /// you use Stewart's randomization procedure).  Stewart gives an
  /// algorithm for recovering from an orthogonalization fault, but
  /// the algorithm is expensive: it involves careful
  /// reorthogonalization with non-block Gram-Schmidt.  We choose
  /// instead to report the orthogonalization fault and let users
  /// recover from it.
  ///
  /// \note This is not a (subclass of) TsqrOrthoError, because the
  /// latter is a logic or runtime bug, whereas a TsqrOrthoFault is a
  /// property of the input and admits recovery.
  class TsqrOrthoFault : public OrthoError
  {
  public: 
    TsqrOrthoFault (const std::string& what_arg) : 
      OrthoError (what_arg) {}
  };


  class NonNullOperatorError : public OrthoError
  {
  public:
    NonNullOperatorError () : 
      OrthoError ("Sorry, TsqrOrthoManager doesn\'t work with a non-null Op "
		  "argument.  I know this is bad class design, but it will "
		  "have to do for now, since TsqrOrthoManager has to inherit "
		  "from MatOrthoManager in order to work in Anasazi\'s "
		  "eigensolvers.  If you want to solve this problem yourself, "
		  "the thing to do is to have TsqrOrthoManager degrade to "
		  "SVQBOrthoManager when this->getOp() != Teuchos::null.")
    {}
  };

  /// \class TsqrOrthoManager
  /// \brief MatOrthoManager implementation using TSQR
  /// 
  /// TsqrOrthoManager is a subclass of MatOrthoManager that uses a
  /// combination of Tall Skinny QR (TSQR) and Block Gram-Schmidt
  /// (BGS) to orthogonalize the columns of groups of multivectors.
  ///
  /// The Block Gram-Schmidt procedure used here is inspired by that
  /// of G. W. Stewart ("Block Gram-Schmidt Orthogonalization", SISC
  /// vol 31 #1 pp. 761--775, 2008), except that we use TSQR+SVD
  /// instead of standard Gram-Schmidt with orthogonalization to
  /// handle the current block.  "Orthogonalization faults" may still
  /// happen, but we do not handle them by default.  Rather, we make
  /// one BGS pass, do TSQR+SVD, check the resulting column norms, and
  /// make a second BGS pass (+ TSQR+SVD) if necessary.  If we then
  /// detect an orthogonalization fault, we throw TsqrOrthoFault.
  ///
  /// \note Only use Teuchos::null for the Op, since TSQR can only
  ///   orthogonalize with respect to the Euclidean norm.
  template< class ScalarType, class MV >
  class TsqrOrthoManagerImpl // public OrthoManager< ScalarType, MV >
  {
  private:
    typedef typename Teuchos::ScalarTraits< ScalarType >::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits< ScalarType >    SCT;
    typedef Teuchos::ScalarTraits< MagnitudeType > SCTM;
    typedef MultiVecTraits< ScalarType, MV >       MVT;

    typedef typename Anasazi::TsqrAdaptor< ScalarType, MV >::adaptor_type tsqr_adaptor_type;
    typedef Teuchos::RCP< tsqr_adaptor_type > tsqr_adaptor_ptr;

  public:

    /// \brief Constructor
    ///
    /// \param tsqrParams [in] Configuration parameters for TSQR.  See
    ///   TSQR documentation for how to set those.  They depend on
    ///   which multivector (MV) class you are using (since each MV
    ///   class maps to a specific TSQR implementation).
    ///
    /// \param Op [in] Inner product.  Don't set to anything not 
    ///   Teuchos::null, otherwise an exception will be thrown.
    ///   Also, don't call setOp().
    ///
    TsqrOrthoManagerImpl (const Teuchos::ParameterList& tsqrParams) :
      tsqrParams_ (tsqrParams),
      tsqrAdaptor_ (Teuchos::null),
      Q_ (Teuchos::null),
      eps_ (SCT::eps()), // ScalarTraits< ScalarType >::eps() returns MagnitudeType
      blockReorthogThreshold_ (MagnitudeType(1) / MagnitudeType(2)),
      relativeRankTolerance_ (MagnitudeType(100)*SCT::eps()),
      throwOnReorthogFault_ (true)
    {}

    // void 
    // setOp (Teuchos::RCP< const OP > Op)
    // {
    //   static_cast< MatOrthoManager< ScalarType, MV, OP >* const >(this)->setOp (Op);
    //   if (getOp() != Teuchos::null) 
    // 	throw NonNullOperatorError();
    // }

    typedef Teuchos::RCP< MV >       mv_ptr;
    typedef Teuchos::RCP< const MV > const_mv_ptr;
    typedef Teuchos::Array< const_mv_ptr >                       const_prev_mvs_type;
    typedef Teuchos::SerialDenseMatrix< int, ScalarType >        serial_matrix_type;
    typedef Teuchos::RCP< serial_matrix_type >                   serial_matrix_ptr;
    typedef Teuchos::Array< Teuchos::RCP< serial_matrix_type > > prev_coeffs_type;

    void 
    innerProd (const MV& X, 
	       const MV& Y, 
	       serial_matrix_type& Z) const
    {
      MVT::MvTransMv (SCT::one(), X, Y, Z);
    }

    void
    norm (const MV& X, 
	  std::vector< MagnitudeType >& normvec) const
    {
      const int ncols = MVT::GetNumberVecs (X);

      // mfh 20 Jul 2010: MatOrthoManager::normMat() computes norms
      // one column at a time, which results in too much
      // communication.  Instead, we compute the whole inner product,
      // which takes more computation but is the only option
      // available, given the current MVT interface.
      serial_matrix_type XTX (ncols, ncols); 
      innerProd (X, X, XTX);

      // Record the results.  Make sure normvec is long enough.
      if (normvec.size() < ncols)
	normvec.resize (ncols);
      for (int k = 0; k < ncols; ++k)
	normvec[k] = SCTM::squareroot (XTX(k,k));
    }

    /// \brief Compute \f$C := Q^* X\f$ and \f$X := X - Q C\f$.
    ///
    /// \warning This method does not do reorthogonalization.  This is
    /// because the reorthogonalization test requires normalization as
    /// well.
    void 
    project (MV& X, 
	     const_prev_mvs_type Q,
	     prev_coeffs_type C = Teuchos::tuple (serial_matrix_ptr (Teuchos::null))) const
    {
      int nrows_X, ncols_X, num_Q_blocks, ncols_Q_total;
      checkProjectionDims (nrows_X, ncols_X, num_Q_blocks, ncols_Q_total, X, Q);
      // Test for quick exit: any dimension of X is zero, or there are
      // zero Q blocks, or the total number of columns of the Q blocks
      // is zero.
      if (nrows_X == 0 || ncols_X == 0 || num_Q_blocks == 0 || ncols_Q_total == 0)
	return;

      // If we don't have enough C, expanding it creates null references.
      // If we have too many, resizing just throws away the later ones.
      // If we have exactly as many as we have Q, this call has no effect.
      C.resize (num_Q_blocks);
      for (int i = 0; i < num_Q_blocks; ++i) 
	{
	  const int ncols_Q_i = MVT::GetNumberVecs (*Q[i]);
	  // Create a new C[i] if necessary.
	  if (C[i] == Teuchos::null)
	    C[i] = Teuchos::rcp (new serial_matrix_type (ncols_Q_i, ncols_X));
	}
      rawProject (X, Q, C);
    }


    int 
    normalize (MV& X,
	       serial_matrix_ptr B = Teuchos::null) const
    {
      // Internal data used by this method require a specific MV
      // object for initialization (e.g., to get a Map / communicator,
      // and to initialize scratch space).  Thus, we delay (hence
      // "lazy") initialization until we get an X.
      lazyInit (X);
      
      // MVT returns int for these two quantities, even though
      // local_ordinal_type of the MV may be some other type.
      const int nrows = MVT::GetVecLength (X);
      const int ncols = MVT::GetNumberVecs (X);

      // TSQR's rank-revealing part doesn't work unless B is provided.
      // If B is not provided, allocate a temporary B for use in TSQR.
      // If it is provided, adjust dimensions as necessary.
      const bool B_is_null_on_input = (B == Teuchos::null);
      if (B_is_null_on_input)
	B = Teuchos::rcp (new serial_matrix_type (ncols, ncols));
      else
	B->reshape (ncols, ncols);
      //
      // Compute rank-revealing decomposition (in this case, TSQR of X
      // followed by SVD of the R factor and appropriate updating of
      // the resulting Q factor) of X.  X is modified in place, and Q
      // contains the results.
      //
      // Compute TSQR and SVD of X.  Resulting orthogonal vectors go
      // into Q_, and coefficients (not necessarily upper triangular)
      // go into B.
      int rank;
      try {
	typedef typename tsqr_adaptor_type::factor_output_type 
	  factor_output_type;
	factor_output_type factorOutput = tsqrAdaptor_->factor (X, *B);
	tsqrAdaptor_->explicitQ (X, factorOutput, Q_);
	rank = tsqrAdaptor_->revealRank (Q_, B, relativeRankTolerance());
      } catch (std::exception& e) {
	throw OrthoError (e.what());
      }
      // Now we should copy Q_ back into X, but don't do it yet: if we
      // want to fill the last ncols-rank columns with random data, we
      // should do so in Q_, because it's fresh in the cache.
      if (false)
	{
	  // If X did not have full (numerical rank), augment the last
	  // ncols-rank columns of X with random data.
	  const int ncolsToFill = ncols - rank;
	  if (ncolsToFill > 0)
	    {
	      // ind: Indices of columns of X to fill with random data.
	      std::vector< int > fillIndices (ncolsToFill);
	      for (int j = 0; j < ncolsToFill; ++j)
		fillIndices[j] = j + rank;

	      mv_ptr Q_null = MVT::CloneViewNonConst (Q_, fillIndices);
	      MVT::MvRandom (*Q_null);
	      Q_null = Teuchos::null;
	    }
	}
      // MultiVecTraits (MVT) doesn't have a "deep copy from one MV
      // into an existing MV in place" method, but it does have two
      // methods which may be used to this effect:
      // 
      // 1. MvAddMv() (compute X := 1*Q_ + 0*X)
      //
      // 2. SetBlock() (Copy from A to mv by setting the index vector
      //    to [0, 1, ..., GetNumberVecs(mv)-1])
      //
      // MVT doesn't state the aliasing rules for MvAddMv(), but I'm
      // guessing it's OK for B and mv to alias one another (that is
      // the way AXPY works in the BLAS).
      MVT::MvAddMv (ScalarType(1), Q_, ScalarType(0), X, X);

      // Don't deallocate B, even if the user provided Teuchos::null
      // as the B input (or the default B input took over).  The RCP's
      // destructor should work just fine.
      return rank;
    }


    int 
    projectAndNormalize (MV &X,
			 const_prev_mvs_type Q,
			 prev_coeffs_type C = Teuchos::tuple (serial_matrix_ptr (Teuchos::null)),
			 serial_matrix_ptr B = Teuchos::null) const
    {
      // if (_hasOp || _Op != Teuchos::null)
      // 	throw NonNullOperatorError();
      // else if (MX != Teuchos::null)
      // 	throw OrthoError("Sorry, TsqrOrthoManager::projectAndNormalizeMat() "
      // 			 "doesn\'t work with MX non-null yet");

      // Fetch dimensions of X and Q.
      int nrows_X, ncols_X, num_Q_blocks, ncols_Q_total;
      checkProjectionDims (nrows_X, ncols_X, num_Q_blocks, ncols_Q_total, X, Q);
      // Test for quick exit: any dimension of X is zero, or there are
      // zero Q blocks, or the total number of columns of the Q blocks
      // is zero.
      if (nrows_X == 0 || ncols_X == 0 || num_Q_blocks == 0 || ncols_Q_total == 0)
	return;

      // If we don't have enough C, expanding it creates null references.
      // If we have too many, resizing just throws away the later ones.
      // If we have exactly as many as we have Q, this call has no effect.
      C.resize (num_Q_blocks);
      prev_coeffs_type newC (num_Q_blocks);
      for (int i = 0; i < num_Q_blocks; ++i) 
	{
	  const int ncols_Q_i = MVT::GetNumberVecs (*Q[i]);
	  // Create a new C[i] if necessary.
	  if (C[i] == Teuchos::null)
	    C[i] = Teuchos::rcp (new serial_matrix_type (ncols_Q_i, ncols_X));
	  // Fill C[i] with zeros.
	  C[i]->putScalar (ScalarType(0));
	  // Create newC[i] as a clone of C[i].  (All that really
	  // matters is that is has the same dimensions and is filled
	  // with zeros; cloning C[i] accomplishes both in this case.)
	  newC[i] = Teuchos::rcp (new serial_matrix_type (*C[i]));
	}

      // Keep track of the column norms of X, both before and after
      // each orthogonalization pass.
      std::vector< MagnitudeType > normsBeforeFirstPass (ncols_X, MagnitudeType(0));
      std::vector< MagnitudeType > normsAfterFirstPass (ncols_X, MagnitudeType(0));
      std::vector< MagnitudeType > normsAfterSecondPass (ncols_X, MagnitudeType(0));
      MVT::MvNorm (X, normsBeforeFirstPass);

      // First BGS pass.  "Modified Gram-Schmidt" version of Block
      // Gram-Schmidt:
      //
      // \li \f$C^{\text{new}}_i := Q_i^* \cdot X\f$
      // \li \f$X := X - Q_i \cdot C^{\text{new}}_i\f$
      rawProject (X, Q, newC);

      // Update the C matrices:
      //
      // \li \f$C_i := C_i + C^{\text{new}}_i\f$
      for (int i = 0; i < num_Q_blocks; ++i)
	*C[i] += *newC[i];

      // Normalize the matrix X.
      int rank = normalize (X, B);
      
      // Compute post-first-pass (pre-normalization) norms, using B.
      // normalize() doesn't guarantee in general that B is upper
      // triangular, so we compute norms using the entire column of B.
      Teuchos::BLAS< int, ScalarType > blas;
      for (int j = 0; j < ncols_X; ++j)
	{
	  const ScalarType* const B_j = &(*B)(0,j);
	  // Teuchos::BLAS returns a MagnitudeType result on
	  // ScalarType inputs.
	  normsAfterFirstPass[j] = blas.NRM2 (ncols_X, B_j, 1);
	}
      // Test whether any of the norms dropped below the
      // reorthogonalization threshold.
      bool reorthog = false;
      for (int j = 0; j < ncols_X; ++j)
	if (normsBeforeFirstPass[j] / normsAfterFirstPass[j] <= blockReorthogThreshold())
	  {
	    reorthog = true; 
	    break;
	  }

      // Perform another BGS pass if necessary.  "Twice is enough"
      // (Kahan's theorem) for a Krylov method, unless (using
      // Stewart's term) there is an "orthogonalization fault"
      // (indicated by reorthogFault).
      bool reorthogFault = false;
      // Indices of X at which there was an orthogonalization fault.
      std::vector< int > faultIndices (0);
      if (reorthog)
	{
	  // Block Gram-Schmidt (again):
	  //
	  // \li \f$C^{\text{new}} = Q^* X\f$
	  // \li \f$X := X - Q C^{\text{new}}\f$
	  // \li \f$C := C + C^{\text{new}}\f$
	  rawProject (X, Q, newC);
	  for (int i = 0; i < num_Q_blocks; ++i)
	    *C[i] += *newC[i];

	  // Normalize the matrix X.
	  serial_matrix_ptr B_new (new serial_matrix_type (ncols_X, ncols_X));
	  rank = normalize (X, B_new);
	  *B += *B_new;

	  // Compute post-second-pass (pre-normalization) norms, using
	  // B.  normalize() doesn't guarantee in general that B is
	  // upper triangular, so we compute norms using the entire
	  // column of B.
	  for (int j = 0; j < ncols_X; ++j)
	    {
	      // FIXME Should we use B_new or B here?
	      const ScalarType* const B_j = &(*B)(0,j);
	      // Teuchos::BLAS returns a MagnitudeType result on
	      // ScalarType inputs.
	      normsAfterSecondPass[j] = blas.NRM2 (ncols_X, B_j, 1);
	    }
	  // Test whether any of the norms dropped below the
	  // reorthogonalization threshold.  If so, it's an
	  // orthogonalization fault, which requires expensive
	  // recovery.
	  reorthogFault = false;
	  for (int j = 0; j < ncols_X; ++j)
	    if (normsAfterSecondPass[j] / normsAfterFirstPass[j] <= blockReorthogThreshold())
	      {
		reorthogFault = true; 
		faultIndices.push_back (j);
	      }
	} // if (reorthog) // reorthogonalization pass

      if (reorthogFault)
	{
	  if (throwOnReorthogFault_)
	    {
	      using std::endl;
	      std::ostringstream os;
	      os << "Orthogonalization fault at the following column(s) of X:" << endl;
	      os << "Column\tNorm decrease factor" << endl;
	      for (int k = 0; k < faultIndices.size(); ++k)
		{
		  const int index = faultIndices[k];
		  const MagnitudeType decreaseFactor = 
		    normsAfterSecondPass[index] / normsAfterFirstPass[index];
		  os << index << "\t" << decreaseFactor << endl;
		}
	      throw TsqrOrthoFault (os.str());
	    }
	  else // Slowly reorthogonalize, one vector at a time, the offending vectors of X.
	    {
	      throw std::logic_error ("Not implemented yet");

	      // for (int k = 0; k < faultIndices.size(); ++k)
	      // 	{
	      // 	  // Get a view of the current column of X.
	      // 	  std::vector<int> currentIndex (1, faultIndices[k]);
	      // 	  Teuchos::RCP< MV > X_j = MVT::CloneViewNonConst (X, currentIndex);

	      // 	  // Reorthogonalize X_j against all columns of each Q[i].
	      // 	}
	    }
	}
      return rank;
    }

    /// \brief Return \f$ \| I - X^* \cdot X \|_F \f$
    ///
    /// Return the Frobenius norm of I - X^* X, which is an absolute
    /// measure of the orthogonality of the columns of X.
    MagnitudeType 
    orthonormError (const MV &X) const
    {
      const ScalarType ONE = SCT::one();
      const int ncols = MVT::GetNumberVecs(X);
      serial_matrix_type XTX (ncols, ncols);
      innerProd (X, X, XTX);
      for (int k = 0; k < ncols; ++k)
	XTX(k,k) -= ONE;
      return XTX.normFrobenius();
    }

    MagnitudeType 
    orthogError (const MV &X1, 
		 const MV &X2) const
    {
      const int ncols_X1 = MVT::GetNumberVecs (X1);
      const int ncols_X2 = MVT::GetNumberVecs (X2);
      serial_matrix_type X1_T_X2 (ncols_X1, ncols_X2);
      innerProd (X1, X2, X1_T_X2);
      return X1_T_X2.normFrobenius();
    }

    /// Relative tolerance for triggering a block reorthogonalization.
    /// If any column norm in a block decreases by this amount, then
    /// we reorthogonalize.
    MagnitudeType blockReorthogThreshold() const { return blockReorthogThreshold_; }
    /// Relative tolerance for determining (via the SVD) whether a
    /// block is of full numerical rank.
    MagnitudeType relativeRankTolerance() const { return relativeRankTolerance_; }

  private:
    /// 
    /// Parameters for initializing TSQR
    Teuchos::ParameterList tsqrParams_;
    ///
    /// Interface between Anasazi and TSQR (the Tall Skinny QR
    /// factorization).
    tsqr_adaptor_ptr tsqrAdaptor_;
    ///
    /// Scratch space for TSQR
    mv_ptr Q_;
    /// 
    // Machine precision for ScalarType
    MagnitudeType eps_;
    ///
    /// Reorthogonalization threshold (relative) in Block Gram-Schmidt.
    MagnitudeType blockReorthogThreshold_;
    ///
    /// Relative tolerance for measuring the numerical rank of a matrix.
    MagnitudeType relativeRankTolerance_;
    
    /// Whether to throw an exception when an orthogonalization fault
    /// occurs.  Recovery is possible, but expensive.
    bool throwOnReorthogFault_;

    /// \brief Initialize the TSQR adaptor and scratch space
    ///
    /// Initialize the TSQR adaptor and scratch space for TSQR.  Both
    /// require a specific MV object, so we have to delay their
    /// initialization until we get an X input (for normalize(), since
    /// only that method uses the TSQR adaptor and the scratch space).
    /// (Hence, "lazy," for delayed initialization.)
    void
    lazyInit (const MV& X)
    {
      // The TSQR adaptor object requires a specific MV object for
      // initialization.  As long as subsequent MV objects use the
      // same communicator (e.g., the same Teuchos::Comm<int>), we
      // don't need to reinitialize the adaptor.
      //
      // FIXME (mfh 15 Jul 2010) If tsqrAdaptor_ has already been
      // initialized, check to make sure that X has the same
      // communicator as the multivector previously used to initialize
      // tsqrAdaptor_.
      if (tsqrAdaptor_ == Teuchos::null)
	tsqrAdaptor_ = Teuchos::rcp (X, tsqrParams_);

      const int nrows = MVT::GetVecLength (X);
      const int ncols = MVT::GetNumberVecs (X);

      // Q_ is temporary workspace.  It must have the same dimensions
      // as X.  If not, we have to reallocate.  We also have to
      // allocate (not "re-") if we haven't allocated Q_ before.  (We
      // can't allocate Q_ until we have some X, so we need a
      // multivector as the "prototype.")
      if (Q_ == Teuchos::null || 
	  nrows != MVT::GetVecLength(X) || 
	  ncols != MVT::GetVecLength(X))
	// Allocate Q_ to have the same dimensions as X.  Contents of
	// X are not copied.
	Q_ = MVT::Clone (X, ncols);
    }

    void
    checkProjectionDims (int& nrows_X, 
			 int& ncols_X, 
			 int& num_Q_blocks,
			 int& ncols_Q_total,
			 const MV& X, 
			 const_prev_mvs_type Q) const
    {
      // Test for quick exit
      num_Q_blocks = Q.length();
      if (num_Q_blocks == 0)
	return;
      nrows_X = MVT::GetVecLength (X);
      ncols_X = MVT::GetNumberVecs (X);
      if (nrows_X == 0 || ncols_X == 0)
	return;

      //
      // Make sure that for each i, the dimensions of X and Q[i] are
      // compatible.
      //
      ncols_Q_total = 0; // total over all Q blocks
      for (int i = 0; i < num_Q_blocks; ++i)
	{
	  const int nrows_Q = MVT::GetVecLength (*Q[i]);
	  TEST_FOR_EXCEPTION( (nrows_Q != nrows_X), 
			      std::invalid_argument,
			      "Anasazi::TsqrOrthoManager::project(): "
			      "Size of X not consistant with size of Q" );
	  ncols_Q_total += MVT::GetNumberVecs (*Q[i]);
	}
    }

    void
    rawProject (MV& X, 
		const_prev_mvs_type Q,
		prev_coeffs_type C = Teuchos::tuple (serial_matrix_ptr (Teuchos::null))) const
    {
      int nrows_X, ncols_X, num_Q_blocks, ncols_Q_total;
      checkProjectionDims (nrows_X, ncols_X, num_Q_blocks, ncols_Q_total, X, Q);
      // Test for quick exit: any dimension of X is zero, or there are
      // zero Q blocks, or the total number of columns of the Q blocks
      // is zero.
      if (nrows_X == 0 || ncols_X == 0 || num_Q_blocks == 0 || ncols_Q_total == 0)
	return;

      // "Modified Gram-Schmidt" version of Block Gram-Schmidt.
      for (int i = 0; i < num_Q_blocks; ++i)
	{
	  innerProd (*Q[i], X, *C[i]);
	  MVT::MvTimesMatAddMv (ScalarType(-1), *Q[i], *C[i], ScalarType(1), X);
	}
    }
  };



  template< class ScalarType, class MV >
  class TsqrOrthoManager : public OrthoManager< ScalarType, MV > {
  public:
    TsqrOrthoManager (const Teuchos::ParameterList& tsqrParams) :
      impl_ (tsqrParams) 
    {}

    virtual ~TsqrOrthoManager() {}

    virtual void 
    innerProd (const MV &X, 
	       const MV &Y, 
	       Teuchos::SerialDenseMatrix<int,ScalarType>& Z) const
    {
      return impl_.innerProd (X, Y, Z);
    }

    virtual void 
    norm (const MV& X, 
	  std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec) const
    {
      return impl_.norm (X, normvec);
    }

    virtual void 
    project (MV &X, 
	     Teuchos::Array<Teuchos::RCP<const MV> > Q,
	     Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C 
	     = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null))) const
    {
      return impl_.project (X, Q, C);
    }

    virtual int 
    normalize (MV &X, 
	       Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null) const
    {
      return impl_.normalize (X, B);
    }

    virtual int 
    projectAndNormalize (MV &X, 
			 Teuchos::Array<Teuchos::RCP<const MV> > Q,
			 Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C 
			 = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
			 Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null) const
    {
      return impl_.projectAndNormalize (X, Q, C, B);
    }

    virtual typename Teuchos::ScalarTraits< ScalarType >::magnitudeType 
    orthonormError (const MV &X) const
    {
      return impl_.orthonormError (X);
    }

    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthogError (const MV &X1, 
		 const MV &X2) const 
    {
      return impl_.orthogError (X1, X2);
    }

  private:
    TsqrOrthoManagerImpl< ScalarType, MV > impl_;
  };

  // template< class ScalarType, class MV >
  // class TsqrMatOrthoManager : public MatOrthoManager< ScalarType, MV > {
  // public:
  //   TsqrMatOrthoManager (const Teuchos::ParameterList& tsqrParams,
  // 			 Teuchos::RCP<const OP> Op = Teuchos::null) :
  //     MatOrthoManager< ScalarType, MV > (Op),
  //     impl_ (tsqrParams)
  //   {
  //     if (getOp() != Teuchos::null)
  // 	throw NonNullOperatorError();
  //   }




  // private:
  //   TsqrOrthoManagerImpl impl_;
  // };


} // namespace Anasazi

#endif // __AnasaziTsqrOrthoManager_hpp

