/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_DETAILS_DENSESOLVER_DECL_HPP
#define IFPACK2_DETAILS_DENSESOLVER_DECL_HPP

/// \file Ifpack2_DenseContainer_decl.hpp
/// \brief Ifpack2::DenseContainer class declaration

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"


namespace Ifpack2 {
namespace Details {

/// \struct LapackSupportsScalar
/// \brief Trait for whether LAPACK supports the given scalar type.
/// \tparam ScalarType Scalar type to test.
///
/// This is an implementation detail of DenseSolver.  It might be
/// useful to promote this to a commonly used utility, but for now,
/// I'll leave it here.
template<class ScalarType>
struct LapackSupportsScalar {
  //! Whether LAPACK supports \c ScalarType.
  static const bool value = false;
};

template<>
struct LapackSupportsScalar<float> {
  static const bool value = true;
};

template<>
struct LapackSupportsScalar<double> {
  static const bool value = true;
};

// FIXME (mfh 15 Nov 2013) Should we write IFPACK2_HAVE_COMPLEX ?
#ifdef TEUCHOS_HAVE_COMPLEX
template<>
struct LapackSupportsScalar<std::complex<float> > {
  static const bool value = true;
};

template<>
struct LapackSupportsScalar<std::complex<double> > {
  static const bool value = true;
};
#endif // TEUCHOS_HAVE_COMPLEX


/// \class DenseSolver
/// \brief "Preconditioner" that uses LAPACK's dense LU.
/// \tparam MatrixType Specialization of Tpetra::RowMatrix.
/// \tparam stub Whether this is a stub implementation.  The default
///   is false.  If true, then this class does nothing and its
///   constructor throws an exception.  You should always use the
///   default value.  This template parameter is an implementation
///   detail; it prevents this class from using LAPACK for types that
///   LAPACK does not support.
///
/// \warning This class only works for the four types supported by
///   LAPACK.  Using any other type will result in a run-time error.
template<class MatrixType,
         const bool stub = ! LapackSupportsScalar<typename MatrixType::scalar_type>::value>
class DenseSolver :
    public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                   typename MatrixType::local_ordinal_type,
                                   typename MatrixType::global_ordinal_type,
                                   typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{};

//! Partial specialization for stub=false (the default).
template<class MatrixType>
class DenseSolver<MatrixType, false> :
    public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                   typename MatrixType::local_ordinal_type,
                                   typename MatrixType::global_ordinal_type,
                                   typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{
public:
  //! \name Public typedefs
  //@{

  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization
  /// (preferred) or a Tpetra::CrsMatrix specialization.
  typedef MatrixType matrix_type;

  //! The type of entries in the input (global) matrix.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input (global) matrix.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input (global) matrix.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Kokkos Node type of the input (global) matrix.
  typedef typename MatrixType::node_type node_type;

  //! The type of the absolute value (magnitude) of a \c scalar_type.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Specialization of Tpetra::RowMatrix used by this class.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  //! Specialization of Tpetra::Map used by this class.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  //@}
  //! \name Constructor and destructor
  //@{

  /// \brief Constructor.
  ///
  /// \param matrix [in] The original input matrix.
  DenseSolver (const Teuchos::RCP<const row_matrix_type>& matrix);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~DenseSolver ();

  //@}
  //! Implementation of Tpetra::Operator
  //@{


  /// \brief The domain Map of this operator.
  ///
  /// The domain Map describes the distribution of valid input vectors
  /// X to the apply() method.
  Teuchos::RCP<const map_type> getDomainMap () const;

  /// \brief The range Map of this operator.
  ///
  /// The range Map describes the distribution of valid output vectors
  /// Y to the apply() method.
  Teuchos::RCP<const map_type> getRangeMap () const;

  /// \brief Apply the preconditioner to X, putting the result in Y.
  ///
  /// If the result of applying this preconditioner to a vector X is
  /// \f$M^{-1} \cdot X$, then this method computes \f$\beta Y + \alpha M^{-1} \cdot X\f$.
  /// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //@}

  //! Set the solvers's parameters.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Set up the graph structure of this preconditioner.
  ///
  /// If the graph structure of the constructor's input matrix has
  /// changed, or if you have not yet called initialize(), you must
  /// call initialize() before you may call compute() or apply().
  ///
  /// Thus, initialize() corresponds to the "symbolic factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  void initialize ();

  //! True if the preconditioner has been successfully initialized, else false.
  bool isInitialized () const;

  /// \brief Set up the numerical values in this preconditioner.
  ///
  /// If the values of the constructor's input matrix have changed, or
  /// if you have not yet called compute(), you must call compute()
  /// before you may call apply().
  ///
  /// Thus, compute() corresponds to the "numeric factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  void compute ();

  //! True if the preconditioner has been successfully computed, else false.
  bool isComputed () const;

  /// \brief Compute the condition number estimate and return its value.
  ///
  /// \warning This method is DEPRECATED.  It was inherited from
  ///   Ifpack, and Ifpack never clearly stated what this method
  ///   computes.  Furthermore, Ifpack's method just estimates the
  ///   condition number of the matrix A, and ignores the
  ///   preconditioner -- which is probably not what users thought it
  ///   did.  If there is sufficient interest, we might reintroduce
  ///   this method with a different meaning and a better algorithm.
  virtual magnitude_type TEUCHOS_DEPRECATED
  computeCondEst (CondestType CT = Ifpack2::Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const row_matrix_type>& Matrix = Teuchos::null);

  /// \brief Return the computed condition number estimate, or -1 if not computed.
  ///
  /// \warning This method is DEPRECATED.  See warning for computeCondEst().
  virtual magnitude_type TEUCHOS_DEPRECATED getCondEst() const;

  //! The input matrix given to the constructor.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! Change the matrix to precondition.
  void setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //! The number of calls to initialize().
  int getNumInitialize () const;

  //! The number of calls to compute().
  int getNumCompute () const;

  //! The number of calls to apply().
  int getNumApply() const;

  //! The total time (in seconds) spent in initialize().
  double getInitializeTime() const;

  //! The total time (in seconds) spent in compute().
  double getComputeTime() const;

  //! The total time (in seconds) spent in apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  std::string description () const;

  //! Print the object with some verbosity level to the given FancyOStream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}
private:

  //! Implementation of describe() for "local" (per process) data.
  void
  describeLocal (Teuchos::FancyOStream& out,
                 const Teuchos::EVerbosityLevel verbLevel) const;

  //! Reset stored data.
  void reset ();

  /// \brief Extract the local dense matrix from the local sparse matrix.
  ///
  /// \param A_local_dense [out] The local dense matrix.
  /// \param A_local [out] The local sparse matrix; the result of
  ///   applying a LocalFilter to the original input matrix A.
  ///   (If A is already "local," then this need not necessarily
  ///   be the result of a LocalFilter.)
  static void
  extract (Teuchos::SerialDenseMatrix<int, scalar_type>& A_local_dense,
           const row_matrix_type& A_local);

  /// \brief Factor the local dense matrix A in place.
  ///
  /// \param A [in/out] On input: The local dense matrix to factor.
  ///   On output: The resulting LU factorization.
  ///
  /// \param ipiv [out] The pivot array from the LU factorization
  ///   (with partial pivoting).
  ///
  /// Call this after calling extract().
  static void
  factor (Teuchos::SerialDenseMatrix<int, scalar_type>& A,
          const Teuchos::ArrayView<int>& ipiv);

  //! Specialization of Tpetra::MultiVector used in implementation.
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;

  //! Specialization of Tpetra::Import used in implementation.
  typedef Tpetra::Import<local_ordinal_type,
                         global_ordinal_type, node_type> import_type;

  //! Specialization of Tpetra::Export used in implementation.
  typedef Tpetra::Export<local_ordinal_type,
                         global_ordinal_type, node_type> export_type;

  //! Specialization of Teuchos::ScalarTraits used in implementation.
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  /// \brief Post-permutation, post-view version of apply().
  ///
  /// apply() first does any necessary subset permutation and view
  /// creation (or copying data), then calls this method to solve the
  /// linear system with the diagonal block.
  ///
  /// \param X [in] Subset permutation of the input X of apply().
  /// \param Y [in] Subset permutation of the input/output Y of apply().
  void
  applyImpl (const MV& X,
             MV& Y,
             const Teuchos::ETransp mode,
             const scalar_type alpha,
             const scalar_type beta) const;

  //! The (original) input matrix.
  Teuchos::RCP<const row_matrix_type> A_;

  //! The result of a LocalFilter on A_.
  Teuchos::RCP<const row_matrix_type> A_local_;

  //! The local (dense) matrix, which compute() extracts and factor() factors in place.
  Teuchos::SerialDenseMatrix<int, scalar_type> A_local_dense_;

  //! Permutation array from LAPACK (GETRF).
  Teuchos::Array<int> ipiv_;

  //! The total time (in seconds) spent in initialize().
  double initializeTime_;

  //! The total time (in seconds) spent in compute().
  double computeTime_;

  //! The total time (in seconds) spent in apply().
  mutable double applyTime_;

  //! Total number of successful initialize() calls.
  int numInitialize_;

  //! Total number of successful compute() calls.
  int numCompute_;

  //! Total number of successful apply() calls.
  mutable int numApply_;

  //! If \c true, the container has been successfully initialized.
  bool isInitialized_;

  //! If \c true, the container has been successfully computed.
  bool isComputed_;
};


//! Partial specialization for stub=true.
template<class MatrixType>
class DenseSolver<MatrixType, true> :
    public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                   typename MatrixType::local_ordinal_type,
                                   typename MatrixType::global_ordinal_type,
                                   typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{
public:
  //! \name Public typedefs
  //@{

  /// \brief The first template parameter of this class.
  ///
  /// This must be either a Tpetra::RowMatrix specialization
  /// (preferred) or a Tpetra::CrsMatrix specialization.
  typedef MatrixType matrix_type;

  //! The type of entries in the input (global) matrix.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input (global) matrix.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input (global) matrix.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Kokkos Node type of the input (global) matrix.
  typedef typename MatrixType::node_type node_type;

  //! The type of the absolute value (magnitude) of a \c scalar_type.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Specialization of Tpetra::RowMatrix used by this class.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  //! Specialization of Tpetra::Map used by this class.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  //@}
  //! \name Constructor and destructor
  //@{

  /// \brief Constructor.
  ///
  /// \param matrix [in] The original input matrix.
  DenseSolver (const Teuchos::RCP<const row_matrix_type>& matrix);

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~DenseSolver ();

  //@}
  //! Implementation of Tpetra::Operator
  //@{

  /// \brief The domain Map of this operator.
  ///
  /// The domain Map describes the distribution of valid input vectors
  /// X to the apply() method.
  Teuchos::RCP<const map_type> getDomainMap () const;

  /// \brief The range Map of this operator.
  ///
  /// The range Map describes the distribution of valid output vectors
  /// Y to the apply() method.
  Teuchos::RCP<const map_type> getRangeMap () const;

  /// \brief Apply the preconditioner to X, putting the result in Y.
  ///
  /// If the result of applying this preconditioner to a vector X is
  /// \f$M^{-1} \cdot X$, then this method computes \f$\beta Y + \alpha M^{-1} \cdot X\f$.
  /// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //@}

  //! Set the solvers's parameters.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Set up the graph structure of this preconditioner.
  ///
  /// If the graph structure of the constructor's input matrix has
  /// changed, or if you have not yet called initialize(), you must
  /// call initialize() before you may call compute() or apply().
  ///
  /// Thus, initialize() corresponds to the "symbolic factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  void initialize ();

  //! True if the preconditioner has been successfully initialized, else false.
  bool isInitialized () const;

  /// \brief Set up the numerical values in this preconditioner.
  ///
  /// If the values of the constructor's input matrix have changed, or
  /// if you have not yet called compute(), you must call compute()
  /// before you may call apply().
  ///
  /// Thus, compute() corresponds to the "numeric factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  void compute ();

  //! True if the preconditioner has been successfully computed, else false.
  bool isComputed () const;

  //! Compute the condition number estimate and return its value.
  magnitude_type TEUCHOS_DEPRECATED
  computeCondEst (CondestType CT = Ifpack2::Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const row_matrix_type>& Matrix = Teuchos::null);

  //! Return the computed condition number estimate, or -0 if not computed.
  magnitude_type TEUCHOS_DEPRECATED getCondEst () const;

  //! The input matrix given to the constructor.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! Change the matrix to precondition.
  void setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //! The number of calls to initialize().
  int getNumInitialize () const;

  //! The number of calls to compute().
  int getNumCompute () const;

  //! The number of calls to apply().
  int getNumApply() const;

  //! The total time (in seconds) spent in initialize().
  double getInitializeTime() const;

  //! The total time (in seconds) spent in compute().
  double getComputeTime() const;

  //! The total time (in seconds) spent in apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  std::string description () const;

  //! Print the object with some verbosity level to the given FancyOStream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}
private:
  //! Specialization of Tpetra::MultiVector used in implementation.
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_DENSESOLVER_DECL_HPP
