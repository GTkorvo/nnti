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

#ifndef IFPACK2_KRYLOV_DEF_HPP
#define IFPACK2_KRYLOV_DEF_HPP

namespace Ifpack2 {

template <class MatrixType>
Krylov<MatrixType>::
Krylov (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  // Default values
  IterationType_ (1),
  Iterations_ (5),
  ResidualTolerance_ (0.001), // FIXME (mfh 17 Jan 2014) Make a function of STS::eps()
  BlockSize_ (1),
  ZeroStartingSolution_ (true),
  PreconditionerType_ (1),
  // General
  Condest_ (-STM::one ()),
  IsInitialized_ (false),
  IsComputed_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumMyRows_ (-1), // FIXME (mfh 17 Jan 2014) What does this mean???
  NumGlobalNonzeros_ (0)
{}


template <class MatrixType>
Krylov<MatrixType>::~Krylov() {}


template <class MatrixType>
void Krylov<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null () && A->getComm ()->getSize () == 1 &&
    A->getNodeNumRows () != A->getNodeNumCols (),
    std::runtime_error, "Ifpack2::Krylov::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getNodeNumRows () << " rows and "
    << A->getNodeNumCols () << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  IsInitialized_ = false;
  IsComputed_ = false;
  Condest_ = -STM::one ();

  A_ = A;
}


template <class MatrixType>
void Krylov<MatrixType>::setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Read in parameters
  Ifpack2::getParameter(params, "krylov: iteration type",IterationType_);
  Ifpack2::getParameter(params, "krylov: number of iterations",Iterations_);
  Ifpack2::getParameter(params, "krylov: residual tolerance",ResidualTolerance_);
  Ifpack2::getParameter(params, "krylov: block size",BlockSize_);
  Ifpack2::getParameter(params, "krylov: zero starting solution",ZeroStartingSolution_);
  Ifpack2::getParameter(params, "krylov: preconditioner type",PreconditionerType_);
  params_=params;
  // Separate preconditioner parameters into another list
  //
  // FIXME (mfh 17 Jan 2014) Inner preconditioner's parameters should
  // be a sublist, not part of the main list!!!
  if (PreconditionerType_ == 1) {
    precParams_.set ("relaxation: sweeps",
                     params_.get ("relaxation: sweeps", 1));
    precParams_.set ("relaxation: damping factor",
                     params_.get("relaxation: damping factor", (scalar_type) 1.0));
    precParams_.set ("relaxation: min diagonal value",
                     params_.get ("relaxation: min diagonal value", STS::one ()));
    precParams_.set ("relaxation: zero starting solution",
                     params_.get ("relaxation: zero starting solution", true));
    precParams_.set ("relaxation: backward mode",
                     params_.get ("relaxation: backward mode", false));
  }
  // FIXME (mfh 17 Jan 2014) AdditiveSchwarz's ParameterList no longer
  // takes parameters for its subdomain solver!  You have to pass them
  // into a sublist.
  if (PreconditionerType_ == 2 || PreconditionerType_ == 3) {
    // FIXME (mfh 17 Jan 2014) should be an integer, given how ILUT
    // works!  Furthermore, this parameter does not mean what you
    // think it means.
    precParams_.set ("fact: ilut level-of-fill",
                     params_.get ("fact: ilut level-of-fill", (double) 1.0));
    // FIXME (mfh 17 Jan 2014) scalar_type or magnitude_type? not
    // sure, but double is definitely wrong.
    precParams_.set ("fact: absolute threshold",
                     params_.get ("fact: absolute threshold", (double) 0.0));
    // FIXME (mfh 17 Jan 2014) scalar_type or magnitude_type? not
    // sure, but double is definitely wrong.
    precParams_.set ("fact: relative threshold",
                     params_.get("fact: relative threshold", (double) 1.0));
    // FIXME (mfh 17 Jan 2014) scalar_type or magnitude_type? not
    // sure, but double is definitely wrong.
    precParams_.set ("fact: relax value",
                     params_.get ("fact: relax value", (double) 0.0));
  }
  if (PreconditionerType_ == 3) {
    precParams_.set ("schwarz: compute condest",
                     params_.get ("schwarz: compute condest",true));
    precParams_.set ("schwarz: combine mode",
                     params_.get ("schwarz: combine mode", "Zero"));
    // FIXME (mfh 17 Jan 2014) AdditiveSchwarz only allows setting
    // this to true if Xpetra and Zoltan2 are enabled!  Otherwise it
    // throws an exception if this is true.
    precParams_.set ("schwarz: use reordering",
                     params_.get ("schwarz: use reordering", true));
    precParams_.set ("schwarz: filter singletons",
                     params_.get ("schwarz: filter singletons", false));
    precParams_.set ("schwarz: overlap level",
                     params_.get ("schwarz: overlap level", (int) 0));
  }
}


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Krylov<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Krylov::getComm: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType>::getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Krylov::getDomainMap: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Krylov::getRangeMap: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool Krylov<MatrixType>::hasTransposeApply () const {
  // FIXME (mfh 17 Jan 2014) apply() does not currently work with mode
  // != NO_TRANS, so it's correct to return false here.
  return false;
}


template <class MatrixType>
int Krylov<MatrixType>::getNumInitialize () const {
  return NumInitialize_;
}


template <class MatrixType>
int Krylov<MatrixType>::getNumCompute () const {
  return NumCompute_;
}


template <class MatrixType>
int Krylov<MatrixType>::getNumApply () const {
  return NumApply_;
}


template <class MatrixType>
double Krylov<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template <class MatrixType>
double Krylov<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template <class MatrixType>
double Krylov<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}


template <class MatrixType>
typename Krylov<MatrixType>::magnitude_type
Krylov<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) { // cannot compute right now
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}


template <class MatrixType>
void Krylov<MatrixType>::initialize ()
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> TMV;
  typedef Tpetra::Operator<scalar_type, local_ordinal_type,
                           global_ordinal_type, node_type> TOP;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Krylov::initialize: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // check only in serial
  //
  // FIXME (mfh 17 Jan 2014) Why do we need to check this?  Belos'
  // LSQR doesn't require a square matrix.
  TEUCHOS_TEST_FOR_EXCEPTION(
    getComm ()->getSize () == 1 && A_->getNodeNumRows () != A_->getNodeNumCols (),
    std::runtime_error, "Ifpack2::Krylov::initialize: If the input matrix A's "
    "communicator has only one process, then the matrix must be square.  "
    "Instead, the matrix has " << A_->getNodeNumRows () << " row(s) and "
    << A_->getNodeNumCols () << " column(s).");

  // clear any previous allocation
  IsInitialized_ = false;
  IsComputed_ = false;

  Teuchos::Time timer ("initialize");
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (timer);

    NumMyRows_ = A_->getNodeNumRows ();

    // Belos parameter list
    RCP<ParameterList> belosList = rcp (new ParameterList ("GMRES"));
    belosList->set ("Maximum Iterations", Iterations_);
    belosList->set ("Convergence Tolerance", ResidualTolerance_);

    // FIXME (17 Jan 2014) This whole "preconditioner type" thing is not
    // how we want Krylov to initialize its inner preconditioner.
    // Krylov should be initialized like AdditiveSchwarz: the Factory
    // should create it, in order to avoid circular dependencies.

    if (PreconditionerType_ == 0) {
      // no preconditioner
    }
    else if (PreconditionerType_==1) {
      ifpack2_prec_=rcp (new Relaxation<MatrixType> (A_));
    }
    else if (PreconditionerType_==2) {
      ifpack2_prec_=rcp (new ILUT<MatrixType> (A_));
    }
    else if (PreconditionerType_==3) {
      ifpack2_prec_ = rcp (new AdditiveSchwarz<MatrixType, ILUT<MatrixType> > (A_));
    }
    else if (PreconditionerType_==4) {
      ifpack2_prec_ = rcp (new Chebyshev<MatrixType> (A_));
    }
    if (PreconditionerType_>0) {
      ifpack2_prec_->initialize();
      ifpack2_prec_->setParameters(precParams_);
    }
    belosProblem_ = rcp (new Belos::LinearProblem<scalar_type,TMV,TOP> ());
    belosProblem_->setOperator (A_);

    // FIXME (mfh 17 Jan 2014) It's also nonsense to use enums to pick
    // the Belos solver type.  Users should just name the Belos solver
    // type and give a sublist of parameters.  This could then go
    // straight into the Belos::SolverFactory to create the solver.

    if (IterationType_ == 1) {
      belosSolver_ =
        rcp (new Belos::BlockGmresSolMgr<scalar_type,TMV,TOP> (belosProblem_, belosList));
    }
    else {
      belosSolver_ =
        rcp (new Belos::BlockCGSolMgr<scalar_type,TMV,TOP> (belosProblem_, belosList));
    }
  }
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
void Krylov<MatrixType>::compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Krylov::compute: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // Don't time the initialize(); that gets timed separately.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::Time timer ("compute");
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (timer);
    if (PreconditionerType_ > 0) {
      ifpack2_prec_->compute ();
      belosProblem_->setLeftPrec (ifpack2_prec_);
    }
  }
  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
void Krylov<MatrixType>::
apply (const Tpetra::MultiVector<typename MatrixType::scalar_type,
       typename MatrixType::local_ordinal_type,
       typename MatrixType::global_ordinal_type,
       typename MatrixType::node_type>& X,
       Tpetra::MultiVector<typename MatrixType::scalar_type,
                           typename MatrixType::local_ordinal_type,
                           typename MatrixType::global_ordinal_type,
                           typename MatrixType::node_type>& Y,
       Teuchos::ETransp mode,
       typename MatrixType::scalar_type alpha,
       typename MatrixType::scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::Krylov::apply: You must call compute() before you may call apply().");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::Krylov::apply: The MultiVector inputs X and Y do not have the "
    "same number of columns.  X.getNumVectors() = " << X.getNumVectors ()
    << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");

  // Catch unimplemented cases: alpha != 1, beta != 0, mode != NO_TRANS.
  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha != STS::one (), std::logic_error,
    "Ifpack2::Krylov::apply: alpha != 1 has not been implemented.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    beta != STS::zero (), std::logic_error,
    "Ifpack2::Krylov::apply: zero != 0 has not been implemented.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    mode != Teuchos::NO_TRANS, std::logic_error,
    "Ifpack2::Krylov::apply: mode != Teuchos::NO_TRANS has not been implemented.");

  Teuchos::Time timer ("apply");
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (timer);

    // If X and Y are pointing to the same memory location,
    // we need to create an auxiliary vector, Xcopy
    RCP<const MV> Xcopy;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      Xcopy = rcp (new MV (X));
    } else {
      Xcopy = rcpFromRef (X);
    }

    RCP<MV> Ycopy = rcpFromRef (Y);
    if (ZeroStartingSolution_) {
      Ycopy->putScalar (STS::zero ());
    }

    // Set left and right hand sides for Belos
    belosProblem_->setProblem (Ycopy, Xcopy);
    belosSolver_->solve (); // solve the linear system
  }
  ++NumApply_;
  ApplyTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
std::string Krylov<MatrixType>::description () const
{
  std::ostringstream out;

  out << "\"Ifpack2::Krylov\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: " << this->getObjectLabel () << ", ";
  }
  out << "Initialized: " << (isInitialized () ? "true" : "false")
      << ", Computed: " << (isComputed () ? "true" : "false")
      << ", Global number of rows: " << A_->getGlobalNumRows ()
      << ", Global number of columns: " << A_->getGlobalNumCols ()
      << "}";
  return out.str ();
}


template <class MatrixType>
void Krylov<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  if (vl != VERB_NONE) {
    // describe() always starts with a tab by convention.
    Teuchos::OSTab tab0 (out);
    out << "\"Ifpack2::Krylov\":";

    Teuchos::OSTab tab1 (out);
    if (this->getObjectLabel () != "") {
      out << "Label: " << this->getObjectLabel () << endl;
    }
    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl
        << "Computed: " << (isComputed () ? "true" : "false") << endl
        << "Global number of rows: " << A_->getGlobalNumRows () << endl
        << "Global number of columns: " << A_->getGlobalNumCols () << endl
        << "Matrix:";
    if (A_.is_null ()) {
      out << " null" << endl;
    } else {
      A_->describe (out, vl);
    }
  }
}

} // namespace Ifpack2

#endif /* IFPACK2_KRYLOV_DEF_HPP */
