// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_NATIVE_NORM_STATUS_TEST_HPP
#define BELOS_NATIVE_NORM_STATUS_TEST_HPP

#include "Belos_ToleranceBasedStatusTestBase.hpp"

namespace Belos {

///
/** Defines a status test based on the relative norm of the unscaled
 * unpreconditioned residual.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NativeNormStatusTest : public ToleranceBasedStatusTestBase<Scalar> {
public:

	///
	typedef Teuchos::ScalarTraits<Scalar>  ST;

	/// Set the stream that output will be sent to
	STANDARD_COMPOSITION_MEMBERS( std::ostream, out )

	/// Set the leading string that will be printed at the beginning of each new line of output.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, leadingOutputStr )

	/** @name Constructor/initializers/accessors */
	//@{

	///
	/** Construct to uninitialized.
	 */
	NativeNormStatusTest();

	///
	/** Calls <tt>initialize(tol)</tt>
	 */
	NativeNormStatusTest(
		const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        tol
		);

	///
	/** Calls <tt>initialize(totalNumRhs,tols)</tt>
	 */
	NativeNormStatusTest(
		const int                                                          totalNumRhs
		,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType       tols[]
		);

	//@}

protected:

	/** @name Overridden from AttachStatusTestBase */
	//@{
	///
	void protectedCheckStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
		);
	//@}

private:

	mutable std::vector<typename ST::magnitudeType>  R_native_norms_;
	mutable std::vector<int>                         currRhsIndexes_;

};

// //////////////////////////////////
// Implementation

// Constructor/initializers/accessors

template<class Scalar>
NativeNormStatusTest<Scalar>::NativeNormStatusTest()
{}

template<class Scalar>
NativeNormStatusTest<Scalar>::NativeNormStatusTest(
	const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        tol
	)
{
	this->initialize(tol);
}

template<class Scalar>
NativeNormStatusTest<Scalar>::NativeNormStatusTest(
	const int                                                          totalNumRhs
	,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType       tols[]
	)
{
	this->initialize(totalNumRhs,tols);
}

// Overridden from AttachStatusTestBase

template<class Scalar>
void NativeNormStatusTest<Scalar>::protectedCheckStatus(
	const BasicIterationState<Scalar>         &bis
	,const int                                currBlockSize
	,const int                                currNumRhs
	,EStatusType                              status[]
	)
{
	const std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &tols = this->tols();
#ifdef _DEBUG
	TEST_FOR_EXCEPT( tols.size() > 1 && ( bis.getProblem().getTotalNumRhs() != static_cast<int>(tols.size()) ) );
#endif
	if(static_cast<int>(R_native_norms_.size()) < currBlockSize) R_native_norms_.resize(currBlockSize);
	bis.getCurrNativeResiduals( currBlockSize, &R_native_norms_[0] );
	if(tols.size()==1 && currNumRhs > 1) {
		// One tolerance for all RHSs.
		const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol = tols[0];
		for( int k = 0; k < currNumRhs; ++k ) {
			if( R_native_norms_[k] <= tol )  status[k] = STATUS_CONVERGED;
			else                             status[k] = STATUS_UNCONVERGED;
		}
	}
	else {
		// Must match a tolerance to each specific RHS
		if( static_cast<int>(currRhsIndexes_.size()) < currNumRhs ) currRhsIndexes_.resize(currNumRhs);
		bis.getProblem().getCurrRhsIndexes( currNumRhs, &currRhsIndexes_[0] );
		for( int k = 0; k < currNumRhs; ++k ) {
			const int origRhsIndex = currRhsIndexes_[k];
			if( R_native_norms_[k] <= tols[origRhsIndex-1] )  status[k] = STATUS_CONVERGED;
			else                                              status[k] = STATUS_UNCONVERGED;
		}
	}
}

} // end Belos namespace

#endif // BELOS_NATIVE_NORM_STATUS_TEST_HPP
