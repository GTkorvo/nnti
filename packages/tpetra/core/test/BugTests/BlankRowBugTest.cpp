/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_UnitTestHarness.hpp>

namespace {

  using Teuchos::RCP;
  using Teuchos::Comm;
  using Tpetra::tuple;

  TEUCHOS_UNIT_TEST(CrsMatrix, BlankRowImport)
  {
    typedef Tpetra::Map<int,int>                    Map;
    typedef Tpetra::CrsMatrix<double,int,int> CrsMatrix;
    typedef Tpetra::Import<int,int>              Import;
    RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    // We run this test explicitly in MPI mode with 2 processors as described in
    // the CMakeLists.txt file. This is just asserting that fact.
    TEST_EQUALITY_CONST(comm->getSize(), 2)

    const int rank = comm->getRank();
    RCP<const Map> destRowMap, sourceRowMap;
    if (rank ==0) {
      sourceRowMap = Tpetra::createNonContigMap<int,int>( tuple<int>(0), comm );
    } else {
      sourceRowMap = Tpetra::createNonContigMap<int,int>( tuple<int>(1), comm );
    }
    destRowMap   = Tpetra::createNonContigMap<int,int>( tuple<int>(0,1), comm );

    RCP<CrsMatrix> srcMat = Tpetra::createCrsMatrix<double>(sourceRowMap);
    if (rank == 0) {
      srcMat->insertGlobalValues(0, tuple<int>(0), tuple<double>(1.0) );
    }
    srcMat->fillComplete();
    /*
       srcMat = [1 ] // proc 0
                [  ] // proc 1
     */
    if (rank == 0) {
      TEST_EQUALITY_CONST( srcMat->getNumEntriesInGlobalRow(0), 1 );
    } else {
      TEST_EQUALITY_CONST( srcMat->getNumEntriesInGlobalRow(1), 0 );
    }

    RCP<CrsMatrix> dstMat = Tpetra::createCrsMatrix<double>(destRowMap);
    RCP<const Import> importer = Tpetra::createImport(sourceRowMap, destRowMap);
    // global row 1 in srcMat is empty: this is a null communication to dstMat
    dstMat->doImport(*srcMat, *importer, Tpetra::INSERT);
    /*
       dstMat_p0 = [1 ]
                   [  ]
       dstMat_p1 = [1 ]
                   [  ]
    */
    TEST_EQUALITY_CONST( dstMat->getNumEntriesInGlobalRow(0), 1 );
    TEST_EQUALITY_CONST( dstMat->getNumEntriesInGlobalRow(1), 0 );
  }

}


