/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_IDENTITY_SOLVER_DEF_HPP
#define IFPACK2_IDENTITY_SOLVER_DEF_HPP

#include "Ifpack2_IdentitySolver_decl.hpp"
#include "Ifpack2_Condest.hpp"

namespace Ifpack2 {

template<class MatrixType>
IdentitySolver<MatrixType>::IdentitySolver(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& A)
 : isInitialized_(false),
   isComputed_(false),
   matrix_(A),
   numInitialize_(0),
   numCompute_(0),
   numApply_(0),
   condEst_ (-Teuchos::ScalarTraits<magnitude_type>::one ())
{
}

template<class MatrixType>
IdentitySolver<MatrixType>::~IdentitySolver()
{
}

template<class MatrixType>
void IdentitySolver<MatrixType>::setParameters(const Teuchos::ParameterList& /*params*/)
{
}

template<class MatrixType>
void IdentitySolver<MatrixType>::initialize()
{
  if (isInitialized_) return;

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
void IdentitySolver<MatrixType>::compute()
{
  if (! isInitialized_) {
    initialize ();
  }
  isComputed_ = false;
  if (matrix_.is_null ()) {
    isComputed_ = true;
    return;
  }

  isComputed_ = true;
  ++numCompute_;
}

template<class MatrixType>
void IdentitySolver<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp /*mode*/,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayRCP;

  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::IdentitySolver::apply() ERROR, compute() hasn't been called yet.");

  ++numApply_;
  //copy X in to Y
  ArrayRCP<const scalar_type> const xData = X.getData(0);
  ArrayRCP<scalar_type> yData = Y.getDataNonConst(0);
  for (size_t i=0; i< xData.size(); ++i)
    yData[i] = xData[i];
}

template<class MatrixType>
typename IdentitySolver<MatrixType>::magnitude_type IdentitySolver<MatrixType>::computeCondEst(CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix)
{
  const magnitude_type minusOne = Teuchos::ScalarTraits<magnitude_type>::one ();

  return minusOne;
}

template <class MatrixType>
int IdentitySolver<MatrixType>::getNumInitialize() const {
  return(numInitialize_);
}

template <class MatrixType>
int IdentitySolver<MatrixType>::getNumCompute() const {
  return(numCompute_);
}

template <class MatrixType>
int IdentitySolver<MatrixType>::getNumApply() const {
  return(numApply_);
}

template <class MatrixType>
double IdentitySolver<MatrixType>::getInitializeTime() const {
  return(initializeTime_);
}

template<class MatrixType>
double IdentitySolver<MatrixType>::getComputeTime() const {
  return(computeTime_);
}

template<class MatrixType>
double IdentitySolver<MatrixType>::getApplyTime() const {
  return(applyTime_);
}

template <class MatrixType>
std::string IdentitySolver<MatrixType>::description() const
{
  return std::string("Ifpack2::IdentitySolver");
}

template <class MatrixType>
void IdentitySolver<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  if (verbLevel != Teuchos::VERB_NONE) {
    out << this->description() << std::endl;
    out << "  numApply: " << numApply_ << std::endl;
  }
}

template <class MatrixType>
Teuchos::RCP<const typename IdentitySolver<MatrixType>::map_type> IdentitySolver<MatrixType>::getDomainMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null (), std::runtime_error, "Ifpack2::IdentitySolver::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return matrix_->getDomainMap ();
}

template <class MatrixType>
Teuchos::RCP<const typename IdentitySolver<MatrixType>::map_type> IdentitySolver<MatrixType>::getRangeMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null (), std::runtime_error, "Ifpack2::IdentitySolver::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return matrix_->getRangeMap ();
}

template<class MatrixType>
void IdentitySolver<MatrixType>::setMatrix (const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null () && A->getComm ()->getSize () == 1 &&
    A->getNodeNumRows () != A->getNodeNumCols (),
    std::runtime_error, "Ifpack2::IdentitySolver::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getNodeNumRows () << " rows and "
    << A->getNodeNumCols () << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  isInitialized_ = false;
  isComputed_ = false;
  matrix_ = A;
}

}//namespace Ifpack2

#endif
