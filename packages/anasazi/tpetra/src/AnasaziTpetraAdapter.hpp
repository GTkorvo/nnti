// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_TPETRA_ADAPTER_HPP
#define ANASAZI_TPETRA_ADAPTER_HPP

#include <Kokkos_NodeTrace.hpp>

/// \file AnasaziTpetraAdapter.hpp
/// \brief Adaptor between Anasazi and Tpetra::{MultiVector,Operator}
///
/// Specializations of Anasazi multi-vector and operator traits
/// classes for the Tpetra MultiVector and Operator classes.
///
/// \note TODO: Here, we assume the solver, multivector and operator
/// are all templated on the same scalar.  This assumption should be
/// relaxed, e.g., for implementing multiprecision solvers.

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#include <AnasaziConfigDefs.hpp>
#include <AnasaziTypes.hpp>
#include <AnasaziMultiVecTraits.hpp>
#include <AnasaziOperatorTraits.hpp>
#include <Kokkos_NodeAPIConfigDefs.hpp>

#ifdef HAVE_ANASAZI_TSQR
#  include <Tpetra_TsqrAdaptor.hpp>
#endif // HAVE_ANASAZI_TSQR


namespace Anasazi {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Tpetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Anasazi::MultiVecTraits class using the Tpetra::MultiVector class.

    This interface will ensure that any Tpetra::MultiVector will be accepted by the Anasazi
    templated solvers.  */
  template<class Scalar, class LO, class GO, class Node>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > {
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  public:
    /// \brief Create a new multivector with \c numvecs columns.
    ///
    /// The returned Tpetra::MultiVector has the same Tpetra::Map
    /// (distribution over one or more parallel processes) as \c X.
    /// Its entries are not initialized and have undefined values.
    static Teuchos::RCP<MV> Clone (const MV& X, const int numvecs) {
      return Teuchos::rcp (new MV (X.getMap (), numvecs, false));
    }

    static Teuchos::RCP<MV> CloneCopy (const MV& X)
    {
      // Make a deep copy of X.  The one-argument copy constructor
      // does a shallow copy by default; the second argument tells it
      // to do a deep copy.
      Teuchos::RCP<MV> X_copy (new MV (X, Teuchos::Copy));
      // Make Tpetra::MultiVector use the new view semantics.  This is
      // a no-op for the Kokkos refactor version of Tpetra; it only
      // does something for the "classic" version of Tpetra.  This
      // shouldn't matter because Belos only handles MV through RCP
      // and through this interface anyway, but it doesn't hurt to set
      // it and make sure that it works.
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<MV>
    CloneCopy (const MV& mv, const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Anasazi::MultiVecTraits::CloneCopy(mv,index)";
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::runtime_error, fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= mv.getNumVectors (),
        std::runtime_error,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << mv.getNumVectors () << " of the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in MultiVector::subCopy, so we
      // don't have to check here.
      Teuchos::RCP<MV> X_copy = mv.subCopy (columns ());
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<MV>
    CloneCopy (const MV& mv, const Teuchos::Range1D& index)
    {
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 &&
        index.ubound() < GetNumberVecs(mv);
      if (! validRange) { // invalid range; generate error message
        std::ostringstream os;
        os << "Anasazi::MultiVecTraits::CloneCopy(mv,index=["
           << index.lbound() << "," << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= GetNumberVecs(mv), std::invalid_argument,
          os.str() << "Index range exceeds number of vectors "
          << mv.getNumVectors() << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<MV> X_copy = mv.subCopy (index);
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<MV>
    CloneViewNonConst (MV& mv, const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Anasazi::MultiVecTraits::CloneViewNonConst(mv,index)";
      const size_t numVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= numVecs,
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << numVecs << " in the input MultiVector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in
      // MultiVector::subViewNonConst, so we don't have to check here.
      Teuchos::RCP<MV> X_view = mv.subViewNonConst (columns ());
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                       const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneViewNonConst(mv,index=["
           << index.lbound() << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative inde{x,ices}, which is "
          "not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numCols, std::invalid_argument,
          os.str() << "Index range exceeds number of vectors " << numCols
          << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<MV> X_view = mv.subViewNonConst (index);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<const MV>
    CloneView (const MV& mv, const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits<Scalar, "
        "Tpetra::MultiVector<...> >::CloneView(mv,index)";
      const size_t numVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= numVecs,
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << numVecs << " in the input MultiVector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in MultiVector::subView, so we
      // don't have to check here.
      Teuchos::RCP<const MV> X_view = mv.subView (columns);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<const MV>
    CloneView (const MV& mv, const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange) {
        std::ostringstream os;
        os << "Anasazi::MultiVecTraits::CloneView(mv, index=["
           << index.lbound () << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numCols, std::invalid_argument,
          os.str() << "Index range exceeds number of vectors " << numCols
          << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<const MV> X_view = mv.subView (index);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static int GetVecLength (const MV& mv) {
      return static_cast<int> (mv.getGlobalLength ());
    }

    static int GetNumberVecs( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { return mv.getNumVectors(); }

    static void
    MvTimesMatAddMv (Scalar alpha,
                     const MV& A,
                     const Teuchos::SerialDenseMatrix<int, Scalar>& B,
                     Scalar beta,
                     MV& mv)
    {
      using Teuchos::ArrayView;
      using Teuchos::Comm;
      using Teuchos::rcpFromRef;
      typedef Tpetra::Map<LO, GO, Node> map_type;

#ifdef HAVE_ANASAZI_TPETRA_TIMERS
      const std::string timerName ("Anasazi::MVT::MvTimesMatAddMv");
      Teuchos::RCP<Teuchos::Time> timer =
        Teuchos::TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = Teuchos::TimeMonitor::getNewCounter (timerName);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        timer.is_null (), std::logic_error,
        "Anasazi::MultiVecTraits::MvTimesMatAddMv: "
        "Failed to look up timer \"" << timerName << "\".  "
        "Please report this bug to the Belos developers.");

      // This starts the timer.  It will be stopped on scope exit.
      Teuchos::TimeMonitor timeMon (*timer);
#endif

      // Check if B is 1-by-1, in which case we can just call update()
      if (B.numRows () == 1 && B.numCols () == 1) {
        mv.update (alpha*B(0,0), A, beta);
        return;
      }

      // Create local map
      Teuchos::SerialComm<int> serialComm;
      map_type LocalMap (B.numRows (), A.getMap ()->getIndexBase (),
                         rcpFromRef<const Comm<int> > (serialComm),
                         Tpetra::LocallyReplicated, A.getMap ()->getNode ());
      // encapsulate Teuchos::SerialDenseMatrix data in ArrayView
      ArrayView<const Scalar> Bvalues (B.values (), B.stride () * B.numCols ());
      // create locally replicated MultiVector with a copy of this data
      MV B_mv (rcpFromRef (LocalMap), Bvalues, B.stride (), B.numCols ());
      mv.multiply (Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    /// \brief <tt>mv := alpha*A + beta*B</tt>
    ///
    /// The Tpetra specialization of this method ignores and
    /// completely overwrites any NaN or Inf entries in A.  Thus, it
    /// does <i>not</i> mean the same thing as <tt>mv := 0*mv +
    /// alpha*A + beta*B</tt> in IEEE 754 floating-point arithmetic.
    /// (Remember that NaN*0 = NaN.)
    static void
    MvAddMv (Scalar alpha,
             const MV& A,
             Scalar beta,
             const MV& B,
             MV& mv)
    {
      mv.update (alpha, A, beta, B, Teuchos::ScalarTraits<Scalar>::zero ());
    }

    static void MvScale (MV& mv, Scalar alpha) {
      mv.scale (alpha);
    }

    static void MvScale (MV& mv, const std::vector<Scalar>& alphas) {
      mv.scale (alphas);
    }

    static void
    MvTransMv (Scalar alpha,
               const MV& A,
               const MV& B,
               Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
      using Tpetra::LocallyReplicated;
      using Teuchos::Comm;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      using Teuchos::SerialComm;
      using Teuchos::TimeMonitor;
      typedef Tpetra::Map<LO,GO,Node> map_type;

#ifdef HAVE_ANASAZI_TPETRA_TIMERS
      const std::string timerName ("Anasazi::MVT::MvTransMv");
      RCP<Teuchos::Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        timer.is_null (), std::logic_error, "Anasazi::MvTransMv: "
        "Failed to look up timer \"" << timerName << "\".  "
        "Please report this bug to the Belos developers.");

      // This starts the timer.  It will be stopped on scope exit.
      TimeMonitor timeMon (*timer);
#endif // HAVE_ANASAZI_TPETRA_TIMERS

      // Form alpha * A^H * B, then copy into the SerialDenseMatrix.
      // We will create a multivector C_mv from a a local map.  This
      // map has a serial comm, the purpose being to short-circuit the
      // MultiVector::reduce() call at the end of
      // MultiVector::multiply().  Otherwise, the reduced multivector
      // data would be copied back to the GPU, only to turn around and
      // have to get it back here.  This saves us a round trip for
      // this data.
      const int numRowsC = C.numRows ();
      const int numColsC = C.numCols ();
      const int strideC  = C.stride ();

      // Check if numRowsC == numColsC == 1, in which case we can call dot()
      if (numRowsC == 1 && numColsC == 1) {
        A.dot (B, Teuchos::ArrayView<Scalar> (C.values (), 1));
        return;
      }

      RCP<Comm<int> > serialComm (new SerialComm<int> ());
      // create local map with serial comm
      RCP<const map_type> LocalMap =
        rcp (new map_type (numRowsC, 0, serialComm, LocallyReplicated,
                           A.getMap ()->getNode ()));
      // create local multivector to hold the result
      const bool INIT_TO_ZERO = true;
      MV C_mv (LocalMap, numColsC, INIT_TO_ZERO);

      // multiply result into local multivector
      C_mv.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, alpha, A, B,
                     Teuchos::ScalarTraits<Scalar>::zero ());
      // get comm
      RCP<const Comm<int> > pcomm = A.getMap ()->getComm ();
      // create arrayview encapsulating the Teuchos::SerialDenseMatrix
      Teuchos::ArrayView<Scalar> C_view (C.values (), strideC * numColsC);
      if (pcomm->getSize() == 1) {
        // No accumulation to do; simply extract the multivector data
        // into C.  Extract a copy of the result into the array view
        // (and therefore, the SerialDenseMatrix).
        C_mv.get1dCopy (C_view, strideC);
      }
      else {
        // get a const host view of the data in C_mv
        Teuchos::ArrayRCP<const Scalar> C_mv_view = C_mv.get1dView();
        if (strideC == numRowsC) {
          // sum-all into C
          reduceAll<int, Scalar> (*pcomm, REDUCE_SUM, numColsC*numRowsC,
                                  C_mv_view.getRawPtr (), C_view.getRawPtr ());
        }
        else {
          // sum-all into temp, copy into C
          Teuchos::Array<Scalar> destBuff (numColsC * numRowsC);
          reduceAll<int,Scalar> (*pcomm, REDUCE_SUM, numColsC*numRowsC,
                                 C_mv_view.getRawPtr (), destBuff.getRawPtr ());
          for (int j=0; j < numColsC; ++j) {
            for (int i=0; i < numRowsC; ++i) {
              C_view[strideC*j+i] = destBuff[numRowsC*j+i];
            }
          }
        }
      }
    }

    //! For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
    static void
    MvDot (const MV& A, const MV& B, std::vector<Scalar> &dots)
    {
      const size_t numVecs = A.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        numVecs != B.getNumVectors (), std::invalid_argument,
        "Anasazi::MultiVecTraits::MvDot(A,B,dots): "
        "A and B must have the same number of columns.  "
        "A has " << numVecs << " column(s), "
        "but B has " << B.getNumVectors () << " column(s).");
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        dots.size() < numVecs, std::invalid_argument,
        "Anasazi::MultiVecTraits::MvDot(A,B,dots): "
        "The output array 'dots' must have room for all dot products.  "
        "A and B each have " << numVecs << " column(s), "
        "but 'dots' only has " << dots.size() << " entry(/ies).");
#endif // HAVE_TPETRA_DEBUG

      Teuchos::ArrayView<Scalar> av (dots);
      A.dot (B, av (0, numVecs));
    }

    //! For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    static void
    MvNorm (const MV& mv,
            std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.getNumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(normvec);
      mv.norm2(av(0,mv.getNumVectors()));
    }

    static void
    SetBlock (const MV& A, const std::vector<int>& index, MV& mv)
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      const size_t inNumVecs = A.getNumVectors ();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        inNumVecs < static_cast<size_t> (index.size ()), std::invalid_argument,
        "Anasazi::MultiVecTraits::SetBlock(A,index,mv): 'index' argument must "
        "have no more entries as the number of columns in the input MultiVector"
        " A.  A.getNumVectors() = " << inNumVecs << " < index.size () = "
        << index.size () << ".");
#endif // HAVE_TPETRA_DEBUG
      RCP<MV> mvsub = CloneViewNonConst (mv, index);
      if (inNumVecs > static_cast<size_t> (index.size ())) {
        RCP<const MV> Asub = A.subView (Range1D (0, index.size () - 1));
        Tpetra::deep_copy (*mvsub, *Asub);
      } else {
        Tpetra::deep_copy (*mvsub, A);
      }
    }

    static void
    SetBlock (const MV& A, const Teuchos::Range1D& index, MV& mv)
    {
      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt =
        static_cast<size_t> (Teuchos::OrdinalTraits<int>::max ());
      const bool overflow =
        maxInt < A.getNumVectors () && maxInt < mv.getNumVectors ();
      if (overflow) {
        std::ostringstream os;
        os << "Anasazi::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < A.getNumVectors (), std::range_error, os.str () << "Number "
          "of columns (size_t) in the input MultiVector 'A' overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < mv.getNumVectors (), std::range_error, os.str () << "Number "
          "of columns (size_t) in the output MultiVector 'mv' overflows int.");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors ());
      const int numColsMv = static_cast<int> (mv.getNumVectors ());
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex =
        index.lbound () >= 0 && index.ubound () < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size () <= numColsA;

      if (! validIndex || ! validSource) {
        std::ostringstream os;
        os << "Anasazi::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Range lower bound must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numColsMv, std::invalid_argument,
          os.str() << "Range upper bound must be less than the number of "
          "columns " << numColsA << " in the 'mv' output argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() > numColsA, std::invalid_argument,
          os.str() << "Range must have no more elements than the number of "
          "columns " << numColsA << " in the 'A' input argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Should never get here!");
      }

      // View of the relevant column(s) of the target multivector mv.
      // We avoid view creation overhead by only creating a view if
      // the index range is different than [0, (# columns in mv) - 1].
      Teuchos::RCP<MV> mv_view;
      if (index.lbound () == 0 && index.ubound () + 1 == numColsMv) {
        mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      } else {
        mv_view = CloneViewNonConst (mv, index);
      }

      // View of the relevant column(s) of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first index.size() columns of A.
      Teuchos::RCP<const MV> A_view;
      if (index.size () == numColsA) {
        A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      } else {
        A_view = CloneView (A, Teuchos::Range1D (0, index.size () - 1));
      }

      Tpetra::deep_copy (*mv_view, *A_view);
    }

    static void
    Assign (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      const char errPrefix[] = "Anasazi::MultiVecTraits::Assign(A, mv): ";

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt =
        static_cast<size_t> (Teuchos::OrdinalTraits<int>::max ());
      const bool overflow =
        maxInt < A.getNumVectors () && maxInt < mv.getNumVectors ();
      if (overflow) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < A.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the input multivector 'A' "
          "(a size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < mv.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the output multivector 'mv' "
          "(a size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Should never get here!");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors ());
      const int numColsMv = static_cast<int> (mv.getNumVectors ());
      if (numColsA > numColsMv) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          numColsA > numColsMv, std::invalid_argument,
          errPrefix << "Input multivector 'A' has " << numColsA << " columns, "
          "but output multivector 'mv' has only " << numColsMv << " columns.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
      }
      if (numColsA == numColsMv) {
        Tpetra::deep_copy (mv, A);
      } else {
        Teuchos::RCP<MV> mv_view =
          CloneViewNonConst (mv, Teuchos::Range1D (0, numColsA - 1));
        Tpetra::deep_copy (*mv_view, A);
      }
    }

    static void MvRandom (MV& mv) {
      mv.randomize ();
    }

    static void
    MvInit (MV& mv, const Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero ())
    {
      mv.putScalar (alpha);
    }

    static void MvPrint (const MV& mv, std::ostream& os) {
      mv.print (os);
    }

#ifdef HAVE_ANASAZI_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Tpetra::MultiVector
    typedef Tpetra::TsqrAdaptor<Tpetra::MultiVector<Scalar, LO, GO, Node> > tsqr_adaptor_type;
#endif // HAVE_ANASAZI_TSQR
  };

  //! Partial specialization of OperatorTraits for Tpetra objects.
  template <class Scalar, class LO, class GO, class Node>
  class OperatorTraits<Scalar,
                       Tpetra::MultiVector<Scalar,LO,GO,Node>,
                       Tpetra::Operator<Scalar,LO,GO,Node> >
  {
  public:
    static void
    Apply (const Tpetra::Operator<Scalar,LO,GO,Node>& Op,
           const Tpetra::MultiVector<Scalar,LO,GO,Node>& X,
           Tpetra::MultiVector<Scalar,LO,GO,Node>& Y)
    {
      Op.apply (X, Y, Teuchos::NO_TRANS);
    }
  };

} // end of Anasazi namespace

#endif
// end of file ANASAZI_TPETRA_ADAPTER_HPP
