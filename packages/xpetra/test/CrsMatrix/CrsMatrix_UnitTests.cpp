// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * BlockedCrsMatrix_UnitTests.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_ConfigDefs.hpp>

#ifdef HAVE_XPETRA_EPETRAEXT
// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>

//#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

namespace {

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Xpetra::DefaultPlatform;
  using Xpetra::Matrix;
  using Xpetra::CrsMatrix;
#ifdef HAVE_XPETRA_TPETRA
  using Xpetra::TpetraCrsMatrix; //TMP
#endif
  using Xpetra::Map;

  using Xpetra::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;



  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  /////////////////////////////////////////////////////

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
                  "test-mpi", "test-serial", &testMpi,
                  "Test MPI (if available) or force test of serial.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
    clp.setOption(
                  "error-tol-slack", &errorTolSlack,
                  "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Apply, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA

    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        matrix->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    matrix->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);

    vec->putScalar(1.0);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec_sol =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(matrix->getRangeMap());

    vec_sol->putScalar(0.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);

    vec_sol->putScalar(2.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, -0.5);

    TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Epetra_ReplaceLocalValues, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA

    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
    //Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        matrix->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    matrix->fillComplete();
    matrix->resumeFill();

    Teuchos::Array<GO> indout(1,0);
    Teuchos::Array<Scalar> valout(1,5.0);
    matrix->replaceLocalValues(0, indout.view(0,indout.size()), valout.view(0,valout.size()));
    matrix->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);

    vec->putScalar(1.0);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec_sol =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(matrix->getRangeMap());

    vec_sol->putScalar(0.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vectest =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
    vectest->putScalar(1.0);
    Teuchos::ArrayRCP<Scalar> vectestData = vectest->getDataNonConst(0);
    vectestData[0] = 5.0;

    vec_sol->update(-1.0,*vectest,1.0);

    TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);
#endif
  }

  // just a copy of the Epetra_ReplaceLocalValues test for Tpetra
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Tpetra_ReplaceLocalValues, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_TPETRA
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    typedef Xpetra::Map<LO, GO, Node> map_type;
    typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
    typedef Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node> crs_matrix_factory_type;
    typedef Xpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
    typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> vec_factory_type;
    typedef Xpetra::Vector<Scalar, LO, GO, Node> vec_type;
    typedef typename Teuchos::Array<LO>::size_type size_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;

    out << "Tpetra replaceLocalValues test" << endl;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm ();

    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
    //Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    out << "Create Map and matrix" << endl;

    // generate problem
    LO nEle = 63;
    RCP<const map_type> map = map_factory_type::Build(lib, nEle, 0, comm);

    RCP<crs_matrix_type> matrix = crs_matrix_factory_type::Build (map, 10);
    const LO NumMyElements = map->getNodeNumElements ();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    // Make the matrix the identity matrix.
    out << "Fill matrix by calling insertGlobalValues" << endl;
    for (LO i = 0; i < NumMyElements; ++i) {
      matrix->insertGlobalValues (MyGlobalElements[i],
                                  Teuchos::tuple<GO>(MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar> (1.0));
    }

    out << "Call fillComplete and resumeFill on matrix" << endl;
    matrix->fillComplete();
    matrix->resumeFill();

    // Change the 0,0 local entry, on each process, to be 5.0.
    out << "Modify entries of the matrix using replaceLocalValues, "
      "and test the result before calling fillComplete" << endl;
    Teuchos::Array<LO> indout (1, 0);
    Teuchos::Array<Scalar> valout(1, 5.0);

    // Every process should have a local row index 0.
    TEST_ASSERT( map->isNodeLocalElement (0) );
    TEST_ASSERT( ! matrix->getColMap ().is_null () );

    if (map->isNodeLocalElement (0) && ! matrix->getColMap ().is_null ()) {
      bool validLocalColumnIndices = true;
      for (size_type k = 0; k < indout.size (); ++k) {
        if (! matrix->getColMap ()->isNodeLocalElement (indout[k])) {
          validLocalColumnIndices = false;
          break;
        }
      }
      // Every process should have a local column index 0.
      TEST_ASSERT( validLocalColumnIndices );
      if (validLocalColumnIndices) {
        // Make sure that we are changing the first diagonal entry on
        // this process.  We determine whether a matrix is diagonal
        // using global indices.
        TEST_ASSERT( matrix->getColMap ()->getGlobalElement (indout[0]) ==
                     map->getGlobalElement (0) );
        // Replace the local (0,0) entry with valout[0].  We know from
        // the above test that the local (0,0) entry is the first
        // diagonal entry on the calling process.
        matrix->replaceLocalValues (0, indout.view (0, indout.size ()),
                                    valout.view (0, valout.size ()));
      }

      // Make sure that replaceLocalValues worked, by getting the
      // values in the local row 0.
      const size_t numEnt = matrix->getNumEntriesInLocalRow (0);
      TEST_EQUALITY_CONST( numEnt, static_cast<size_t> (1) );

      if (numEnt == static_cast<size_t> (1)) {
        Teuchos::Array<LO> ind (numEnt);
        Teuchos::Array<Scalar> val (numEnt);
        size_t numEntOut = 0;
        matrix->getLocalRowCopy (0, ind (), val (), numEntOut);
        TEST_EQUALITY( numEnt, numEntOut );

        if (numEntOut == static_cast<size_t> (1)) {
          TEST_EQUALITY( ind[0], 0 );
          TEST_EQUALITY( val[0], 5.0 );
        }
      }
    }

    // Make sure that all processes got this far.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    out << "Call fillComplete on matrix for the second time" << endl;
    matrix->fillComplete ();

    out << "Test the result of replaceLocalValues after fillComplete" << endl;
    if (map->isNodeLocalElement (0)) {
      // Make sure that replaceLocalValues worked, by getting the
      // values in the local row 0.
      const size_t numEnt = matrix->getNumEntriesInLocalRow (0);
      TEST_EQUALITY_CONST( numEnt, static_cast<size_t> (1) );

      if (numEnt == static_cast<size_t> (1)) {
        Teuchos::Array<LO> ind (numEnt);
        Teuchos::Array<Scalar> val (numEnt);
        size_t numEntOut = 0;
        matrix->getLocalRowCopy (0, ind (), val (), numEntOut);
        TEST_EQUALITY( numEnt, numEntOut );

        if (numEntOut == static_cast<size_t> (1)) {
          TEST_EQUALITY( ind[0], 0 );
          TEST_EQUALITY( val[0], 5.0 );
        }
      }
    }

    RCP<vec_type> vec = vec_factory_type::Build (map);
    vec->putScalar (1.0);
    out << "Test that vec->putScalar(1.0) filled vec with ones" << endl;
    {
      const MT N = static_cast<MT> (vec->getGlobalLength ());
      const MT expectedNorm2 = STM::squareroot (N);
      const MT actualNorm2 = vec->norm2 ();
      TEST_EQUALITY( actualNorm2, expectedNorm2 );
    }

    RCP<const map_type> rangeMap = matrix->getRangeMap ();
    TEST_ASSERT( ! rangeMap.is_null () );
    RCP<vec_type> vec_sol = vec_factory_type::Build (rangeMap);
    vec_sol->putScalar (0.0);
    out << "Test that vec_sol->putScalar(0.0) filled vec with zeros" << endl;
    {
      const MT expectedNorm2 = STM::zero ();
      const MT actualNorm2 = vec_sol->norm2 ();
      TEST_EQUALITY( actualNorm2, expectedNorm2 );
    }

    // Compute vec_sol := matrix*vec.  The result _should_ be a vector
    // of ones everywhere, except for the entry at local index zero
    // (on every process), which should be 5.0.
    matrix->apply (*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);
    if (rangeMap->getNodeNumElements () > 0) {
      // Test this both for a const view and for a nonconst view.
      // This may also be a test for {T,X}petra::MultiVector::getData
      // and {T,X}petra::MultiVector::getDataNonConst.

      // Create the const view.
      Teuchos::ArrayRCP<const Scalar> outData = vec_sol->getData (0);
      TEST_ASSERT( outData.size () == rangeMap->getNodeNumElements () );
      if (outData.size () == rangeMap->getNodeNumElements () &&
          outData.size () > static_cast<size_type> (0)) {
        TEST_EQUALITY( outData[0], 5.0 );
      }
      if (rangeMap->getNodeNumElements () > static_cast<size_t> (1)) {
        bool allOnes = true;
        for (size_type k = 1; k < rangeMap->getNodeNumElements (); ++k) {
          if (! outData[k] == 1.0) {
            allOnes = false;
          }
        }
        TEST_ASSERT( allOnes );
      }

      // Invalidate the const view, before creating a nonconst view.
      outData = Teuchos::null;
      // Create the nonconst view.
      Teuchos::ArrayRCP<Scalar> outDataNonConst = vec_sol->getDataNonConst (0);
      TEST_ASSERT( outDataNonConst.size () == rangeMap->getNodeNumElements () );
      if (outDataNonConst.size () == rangeMap->getNodeNumElements () &&
          outDataNonConst.size () > static_cast<size_type> (0)) {
        TEST_EQUALITY( outDataNonConst[0], 5.0 );
      }
      if (rangeMap->getNodeNumElements () > static_cast<size_t> (1)) {
        bool allOnes = true;
        for (size_type k = 1; k < rangeMap->getNodeNumElements (); ++k) {
          if (! outDataNonConst[k] == 1.0) {
            allOnes = false;
          }
        }
        TEST_ASSERT( allOnes );
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (gblSuccess == 1) {
      out << "Vector result is correct" << endl;
    }


    RCP<vec_type> vectest = vec_factory_type::Build (map);
    vectest->putScalar (1.0);
    Teuchos::ArrayRCP<Scalar> vectestData = vectest->getDataNonConst(0);
    vectestData[0] = 5.0;

    vec_sol->update(-1.0,*vectest,1.0);

    TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);

    // Make sure that all processes got this far.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && (gblSuccess == 1);
    TEST_EQUALITY_CONST( gblSuccess, 1 );
#endif
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TpetraDeepCopy, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": (CrsMatrix, TpetraDeepCopy) test" << endl;
      cerr << os.str ();
    }

    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    // Create a Map, which will be the row, domain, and range Map of the matrix A.
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Creating matrix" << endl;
      cerr << os.str ();
    }

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build (map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Filling matrix" << endl;
      cerr << os.str ();
    }

    // Make A the identity matrix.
    for (LO i = 0; i < NumMyElements; ++i) {
      A->insertGlobalValues (MyGlobalElements[i],
                            Teuchos::tuple<GO> (MyGlobalElements[i]),
                            Teuchos::tuple<Scalar> (1.0));
    }

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Calling fillComplete on matrix A" << endl;
      cerr << os.str ();
    }

    A->fillComplete ();

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Building an input Vector" << endl;
      cerr << os.str ();
    }

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (A->getRangeMap ());
    v->setSeed (8675309);
    v->randomize (true);

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Testing that Xpetra::Vector::operator= does a deep copy" << endl;
      cerr << os.str ();
    }

    // Remember the norm of v, to make sure that neither apply() call changes it.
    const magnitude_type v_norm = v->norm2 ();

    // Keep a copy of v, to test that neither apply() call changes it.
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vcopy =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (map);
    // Xpetra's operator= does a deep copy, like Epetra, but unlike
    // Tpetra (as of early 2014).
    *vcopy = *v;

    // Make sure that vcopy and v have the same norm.  It's OK for the
    // norms to be slightly different, due to nondeterminism in
    // parallel collectives.
    const magnitude_type vcopy_norm = vcopy->norm2 ();

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": v->norm2() = " << v_norm
         << ", vcopy->norm2() = " << vcopy_norm << endl;
      cerr << os.str ();
    }

    const magnitude_type norm_tol =
      static_cast<magnitude_type> (map->getGlobalNumElements ()) * STM::eps ();
    TEUCHOS_TEST_COMPARE(STM::magnitude (v_norm - vcopy_norm), <, norm_tol, out, success);

    // Make sure that if you change vcopy, v doesn't change.
    // That is, vcopy must be a true deep copy of v.
    {
      Teuchos::ArrayRCP<Scalar> vcopy_data = vcopy->getDataNonConst (0);
      if (NumMyElements != 0) {
        vcopy_data[0] += static_cast<magnitude_type> (10000.0);
      }
      // Destroy the view, so that the changes get written back to the Vector.
      vcopy_data = Teuchos::null;

      // Adding 10000 to an entry had better change the 2-norm by at least sqrt(10000) = 100.
      const magnitude_type norm_tol2 = static_cast<magnitude_type> (100.0);
      TEUCHOS_TEST_COMPARE(STM::magnitude (vcopy_norm - vcopy->norm2 ()), >, norm_tol2, out, success);

      // Restore the original vcopy, by doing a deep copy again.
      // Xpetra's operator= does a deep copy, like Epetra, but unlike
      // Tpetra (as of early 2014).
      *vcopy = *v;

      // Make sure the original copy got restored.
      TEUCHOS_TEST_COMPARE(STM::magnitude (vcopy_norm - vcopy->norm2 ()), <, norm_tol, out, success);
    }

    // r and rcopy are distinct Vectors with the same Map, namely the
    // range Map of A.  All the Vectors v, r, and rcopy are distinct.
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (A->getRangeMap ());
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy =
      Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (A->getRangeMap ());

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Applying matrix A" << endl;
      cerr << os.str ();
    }

    // r := A * v.
    A->apply (*v, *r, Teuchos::NO_TRANS, STS::one (), STS::zero ());

    // Since A is the identity matrix, after the above line finishes,
    // r and v should be exactly equal.  This should be true even in
    // finite-precision arithmetic.  Test this here.
    {
      RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (A->getRangeMap ());

      // diff := 1*v + (-1)*r.
      diff->update (STS::one (), *v, -STS::one (), *r, STS::zero ());
      Teuchos::Array<magnitude_type> norms (1);
      diff->norm2 (norms ());
      // The norm of v - r must be _exactly_ zero.
      TEST_EQUALITY(norms[0], STM::zero ());
    }

    // Make sure that the above apply() call didn't change v, by
    // testing against vcopy.
    {
      RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (A->getRangeMap ());

      // diff := 1*v + (-1)*vcopy.
      diff->update (STS::one (), *v, -STS::one (), *vcopy, STS::zero ());
      Teuchos::Array<magnitude_type> norms (1);
      diff->norm2 (norms ());
      // The norm of v - vcopy must be _exactly_ zero.
      TEST_EQUALITY(norms[0], STM::zero ());
    }
    // TODO Make sure norm of v didn't change


    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Deep-copying matrix" << endl;
      cerr << os.str ();
    }

    using Teuchos::rcp_static_cast;
    // NOTE (mfh 24 Apr 2014): This invokes the
    // Xpetra::TpetraCrsMatrix copy constructor, on line 329 of
    // Xpetra_TpetraCrsMatrix.hpp (as of 24 Apr 2014).  That in turn
    // calls Tpetra::CrsMatrix::clone().
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > Acopy =
      rcp (new Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node> (* (rcp_static_cast<Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node> > (A))));

    // Make sure that A and Acopy have the same gross properties.  For
    // example, they must be both fill complete and locally indexed,
    // and their four Maps must match.

    const bool bothFillComplete = A->isFillComplete () && Acopy->isFillComplete ();
    TEST_EQUALITY_CONST(bothFillComplete, true);

    const bool bothLocallyIndexed = A->isLocallyIndexed () && Acopy->isLocallyIndexed ();
    TEST_EQUALITY_CONST(bothLocallyIndexed, true);

    const bool bothNotGloballyIndexed = ! A->isGloballyIndexed () && ! Acopy->isGloballyIndexed ();
    TEST_EQUALITY_CONST(bothNotGloballyIndexed, true);

    const bool rowMapsMatch = A->getRowMap ()->isSameAs (* (Acopy->getRowMap ()));
    TEST_EQUALITY_CONST(rowMapsMatch, true);

    const bool colMapsMatch = A->getColMap ()->isSameAs (* (Acopy->getColMap ()));
    TEST_EQUALITY_CONST(colMapsMatch, true);

    const bool domainMapsMatch = A->getDomainMap ()->isSameAs (* (Acopy->getDomainMap ()));
    TEST_EQUALITY_CONST(domainMapsMatch, true);

    const bool rangeMapsMatch = A->getRangeMap ()->isSameAs (* (Acopy->getRangeMap ()));
    TEST_EQUALITY_CONST(rangeMapsMatch, true);

    TEST_EQUALITY(A->getGlobalNumRows (), Acopy->getGlobalNumRows ());
    TEST_EQUALITY(A->getGlobalNumCols (), Acopy->getGlobalNumCols ());
    TEST_EQUALITY(A->getGlobalNumEntries (), Acopy->getGlobalNumEntries ());
    TEST_EQUALITY(A->getGlobalNumDiags (), Acopy->getGlobalNumDiags ());
    TEST_EQUALITY(A->getGlobalMaxNumRowEntries (), Acopy->getGlobalMaxNumRowEntries ());

    // FIXME (mfh 24 Apr 2014) Need to test separately on each MPI
    // process and do an all-reduce to check if all got it right.
    TEST_EQUALITY(A->getNodeNumRows (), Acopy->getNodeNumRows ());
    TEST_EQUALITY(A->getNodeNumCols (), Acopy->getNodeNumCols ());
    TEST_EQUALITY(A->getNodeNumEntries (), Acopy->getNodeNumEntries ());
    TEST_EQUALITY(A->getNodeNumDiags (), Acopy->getNodeNumDiags ());
    TEST_EQUALITY(A->getNodeMaxNumRowEntries (), Acopy->getNodeMaxNumRowEntries ());

    // Acopy and A should be identically the same.  We can verify this
    // in two ways.  First, we can directly compare the rows of both
    // matrices on each process.  Second, we can repeat the apply()
    // operation with Acopy and verify that it produces the same
    // result.  We will take both approaches here.

    // This test only makes sense if the row Maps of A and Acopy
    // match.  Otherwise, some row indices that are valid for one
    // matrix might not be valid for the other.
    if (rowMapsMatch) {
      typedef typename Teuchos::ArrayView<const GO>::size_type size_type;
      // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
      // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
      //
      // Teuchos::Array<GO> A_ginds (A->getNodeMaxNumRowEntries ());
      Teuchos::Array<LO> A_linds (A->getNodeMaxNumRowEntries ());
      Teuchos::Array<Scalar> A_vals (A->getNodeMaxNumRowEntries ());

      // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
      // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
      //
      // Teuchos::Array<GO> Acopy_ginds (Acopy->getNodeMaxNumRowEntries ());
      Teuchos::Array<LO> Acopy_linds (Acopy->getNodeMaxNumRowEntries ());
      Teuchos::Array<Scalar> Acopy_vals (Acopy->getNodeMaxNumRowEntries ());

      for (size_type k = 0; k < static_cast<size_type> (NumMyElements); ++k) {
        const LO lrow = static_cast<LO> (k);
        //const GO grow = MyGlobalElements[k];
        size_t A_numEnt = 0;
        size_t Acopy_numEnt = 0;

        // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
        // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
        //
        // A->getGlobalRowCopy (grow, A_ginds (), A_vals (), A_numEnt);
        // Acopy->getGlobalRowCopy (grow, Acopy_ginds (), Acopy_vals (), Acopy_numEnt);
        // TEST_COMPARE_ARRAYS( A_ginds (0, A_numEnt), Acopy_ginds (0, Acopy_numEnt) );
        // TEST_COMPARE_ARRAYS( A_vals (0, A_numEnt), Acopy_vals (0, Acopy_numEnt) );

        TEST_EQUALITY(A->getNumEntriesInLocalRow (lrow), Acopy->getNumEntriesInLocalRow (lrow));

        // mfh 24 Apr 2014: Apparently, Xpetra::CrsMatrix implements
        // neither getGlobalRowCopy nor getNumEntriesInGlobalRow.
        //
        // TEST_EQUALITY(A->getNumEntriesInGlobalRow (grow), Acopy->getNumEntriesInGlobalRow (grow));

        A->getLocalRowCopy (lrow, A_linds (), A_vals (), A_numEnt);
        Acopy->getLocalRowCopy (lrow, Acopy_linds (), Acopy_vals (), Acopy_numEnt);

        TEST_COMPARE_ARRAYS( A_linds (0, A_numEnt), Acopy_linds (0, Acopy_numEnt) );
        TEST_COMPARE_ARRAYS( A_vals (0, A_numEnt), Acopy_vals (0, Acopy_numEnt) );
      }
    }

    // Make sure that A doesn't exist anymore.
    A = Teuchos::null;

    {
      using std::cerr;
      using std::endl;

      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Process " << myRank << ": Applying matrix Acopy" << endl;
      cerr << os.str ();
    }

    Acopy->apply (*v, *rcopy, Teuchos::NO_TRANS, STS::one (), STS::zero ());

    {
      // Since A is the identity matrix and Acopy is a deep copy of A,
      // after this line finishes, rcopy and v should be exactly
      // equal.  This should be true even in finite-precision
      // arithmetic.  Test this here.
      RCP<Xpetra::Vector<Scalar, LO, GO, Node> > diff =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build (Acopy->getRangeMap ());

      // diff := 1*v + (-1)*rcopy.
      diff->update (STS::one (), *v, -STS::one (), *rcopy, STS::zero ());
      Teuchos::Array<magnitude_type> norms (1);
      diff->norm2 (norms ());
      // The norm of v - r must be _exactly_ zero.
      TEST_EQUALITY(norms[0], STM::zero ());
    }

    // Repeat the above test in a different way.
    Teuchos::ArrayRCP<const Scalar> rdata = r->getData (0);
    Teuchos::ArrayRCP<const Scalar> rdatacopy = rcopy->getData (0);
    magnitude_type s = STM::zero ();
    for (LO i = 0; i < NumMyElements; ++i) {
      s += Teuchos::ScalarTraits<magnitude_type>::magnitude (rdata[i] - rdatacopy[i]);
    }
    TEUCHOS_TEST_COMPARE(s, <, 1e-16, out, success);
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EpetraDeepCopy, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA

    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > A =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        A->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    A->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > v = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);
    v->setSeed(8675309);
    v->randomize(true);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > r = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());
    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > rcopy = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(A->getRangeMap());

    A->apply(*v, *r, Teuchos::NO_TRANS, 1.0, 0.0);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > Acopy(new Xpetra::EpetraCrsMatrix(*(Teuchos::rcp_static_cast<Xpetra::EpetraCrsMatrix>(A))));
    A = Teuchos::null;

    Acopy->apply(*v, *rcopy, Teuchos::NO_TRANS, 1.0, 0.0);

    Teuchos::ArrayRCP<Scalar> rdata = r->getDataNonConst(0), rdatacopy = rcopy->getDataNonConst(0);
    Scalar s = Teuchos::ScalarTraits<Scalar>::zero ();
    for (LO i = 0; i < NumMyElements; i++) {
      s += Teuchos::ScalarTraits<Scalar>::magnitude (rdata[i] - rdatacopy[i]);
    }
    TEUCHOS_TEST_COMPARE(s, <, 1e-16, out, success);
#endif
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Apply, SC, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Epetra_ReplaceLocalValues, SC, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Tpetra_ReplaceLocalValues, SC, LO, GO, Node )
#define UNIT_TEST_GROUP_ORDINAL1( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TpetraDeepCopy, SC, LO, GO, Node )
#define UNIT_TEST_GROUP_ORDINAL2( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EpetraDeepCopy, SC, LO, GO, Node )

  typedef KokkosClassic::DefaultNode::DefaultNodeType DefaultNodeType;

  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)
  UNIT_TEST_GROUP_ORDINAL1(double, int, int, DefaultNodeType)
  UNIT_TEST_GROUP_ORDINAL2(double, int, int, DefaultNodeType)
}

