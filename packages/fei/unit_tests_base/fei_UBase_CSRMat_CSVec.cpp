/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "fei_CSVec.hpp"
#include "fei_CSRMat.hpp"

#include "fei_iostream.hpp"

#include <cmath>
#include <limits>

namespace {

TEUCHOS_UNIT_TEST(CSRMatCSVec, FillableVec_1)
{
  fei::FillableVec fv;

  if (fv.hasEntry(0)) {
    throw std::runtime_error("FillableVec failed 1");
  }

  bool test_passed = true;
  try {
    fv.getEntry(0);
    test_passed = false;
  }
  catch(...) {}

  if (test_passed == false) {
    throw std::runtime_error("FillableVec failed 2");
  }

  fv.addEntry(0, 0.0);
  fv.addEntry(1, 1.0);
  fv.putEntry(2, 2.0);
  fv.addEntry(2, 2.0);

  test_passed = true;
  try {
    double coef = fv.getEntry(2);
    const double fei_eps = std::numeric_limits<double>::epsilon();
    if (std::abs(coef - 4.0) > fei_eps) test_passed = false;
  }
  catch(...) {test_passed = false;}

  if (test_passed == false) {
    throw std::runtime_error("FillableVec failed 3");
  }

  if (!fv.hasEntry(1)) {
    throw std::runtime_error("FillableVec failed 4");
  }

  if (fv.size() != 3) {
    throw std::runtime_error("FillableVec failed 5");
  }
}

TEUCHOS_UNIT_TEST(CSRMatCSVec, FillableMat_1)
{
  fei::FillableMat fm;

  if (fm.hasRow(0)) {
    throw std::runtime_error("FillableMat failed 1");
  }

  bool test_passed = true;
  try {
    fm.getRow(0);
    test_passed = false;
  }
  catch(...) {}

  if (test_passed == false) {
    throw std::runtime_error("FillableMat failed 2");
  }

  fm.sumInCoef(0, 0, 0.0);
  fm.sumInCoef(1, 1, 1.0);
  fm.putCoef(2, 2, 2.0);
  fm.sumInCoef(2, 2, 2.0);

  test_passed = true;
  try {
    const fei::FillableVec* row = fm.getRow(2);
    if (row->size() != 1) test_passed = false;
  }
  catch(...) {test_passed = false;}

  if (test_passed == false) {
    throw std::runtime_error("FillableMat failed 3");
  }

  if (!fm.hasRow(1)) {
    throw std::runtime_error("FillableMat failed 4");
  }

  if (fm.getNumRows() != 3) {
    throw std::runtime_error("FillableMat failed 5");
  }
}

TEUCHOS_UNIT_TEST(CSRMatCSVec, multiply_CSRMat_CSVec)
{
  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(0, 1, 0.1);
  fm.putCoef(1, 0, 1.0);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(2, 0, 2.0);
  fm.putCoef(2, 1, 2.1);
  fm.putCoef(2, 2, 2.2);

  fei::FillableVec fv;

  fv.putEntry(0, 1.0);
  fv.putEntry(1, 2.0);
  fv.putEntry(2, 3.0);

  fei::CSVec x(fv), y;
  fei::CSRMat A(fm);

  fei::multiply_CSRMat_CSVec(A, x, y);

  if (y.size() != 3) {
    throw std::runtime_error("CSRMat y=Ax test 1 failed.");
  }

  std::vector<int>& y_ind = y.indices();
  std::vector<double>& y_coef = y.coefs();

  if (y_ind[1] != 1) {
    throw std::runtime_error("CSRMat y=Ax test 2 failed.");
  }

  if (std::abs(y_coef[1] - 3.2) > 1.e-13) {
    throw std::runtime_error("CSRMat y=Ax test 3 failed.");
  }

  if (std::abs(y_coef[2] - 12.8) > 1.e-13) {
    throw std::runtime_error("CSRMat y=Ax test 4 failed.");
  }

  fei::multiply_trans_CSRMat_CSVec(A, x, y);

  if (y.size() != 3) {
    throw std::runtime_error("CSRMat y=A^Tx test 1 failed.");
  }

  std::vector<int>& y_ind2 = y.indices();
  std::vector<double>& y_coef2 = y.coefs();

  if (y_ind2[1] != 1) {
    throw std::runtime_error("CSRMat y=A^Tx test 2 failed.");
  }

  if (std::abs(y_coef2[1] - 8.6) > 1.e-13) {
    throw std::runtime_error("CSRMat y=A^Tx test 3 failed.");
  }

  if (std::abs(y_coef2[2] - 6.6) > 1.e-13) {
    throw std::runtime_error("CSRMat y=A^Tx test 4 failed.");
  }
}

TEUCHOS_UNIT_TEST(CSRMatCSVec, multiply_CSRMat_CSRMat)
{
  fei::FillableMat fa, fb;

  fa.putCoef(0, 0, 1.0);
  fa.putCoef(0, 1, 1.0);
  fa.putCoef(1, 1, 2.0);
  fa.putCoef(2, 1, 3.0);
  fa.putCoef(2, 2, 3.0);

  fb.putCoef(0, 0, 1.0);
  fb.putCoef(0, 1, 1.0);
  fb.putCoef(1, 1, 2.0);
  fb.putCoef(2, 1, 3.0);
  fb.putCoef(2, 2, 3.0);

  fei::CSRMat A(fa), B(fb), C;

  fei::multiply_CSRMat_CSRMat(A, B, C);

  if (C.getNumRows() != 3) {
    throw std::runtime_error("CSRMAT C=AB test 1 failed.");
  }

  std::vector<int>& cols = C.getGraph().packedColumnIndices;
  std::vector<double>& coefs = C.getPackedCoefs();

  if (cols.size() != 5) {
    throw std::runtime_error("CSRMAT C=AB test 2 failed.");
  }

  if (cols[1] != 1) {
    throw std::runtime_error("CSRMAT C=AB test 3 failed.");
  }

  if (cols[3] != 1) {
    throw std::runtime_error("CSRMAT C=AB test 4 failed.");
  }

  if (std::abs(coefs[1] - 3.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test 5 failed.");
  }

  if (std::abs(coefs[3] - 15.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test 6 failed.");
  }

  if (std::abs(coefs[4] - 9.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test 7 failed.");
  }

  fei::FillableMat fa2, fb2;

  fa2.putCoef(1, 1, 1.0);
  fa2.putCoef(1, 2, 1.0);
  fa2.putCoef(2, 1, 1.0);
  fa2.putCoef(2, 2, 1.0);
  fa2.putCoef(3, 1, 1.0);
  fa2.putCoef(3, 2, 1.0);
  fa2.putCoef(3, 3, 1.0);

  fb2.putCoef(1, 1, 1.0);
  fb2.putCoef(1, 2, 1.0);
  fb2.putCoef(2, 1, 1.0);
  fb2.putCoef(2, 2, 1.0);
  fb2.putCoef(3, 1, 1.0);
  fb2.putCoef(3, 2, 1.0);
  fb2.putCoef(3, 3, 1.0);

  fei::CSRMat A2(fa2), B2(fb2), C2;

  fei::multiply_CSRMat_CSRMat(A2, B2, C2);

  if (C.getNumRows() != 3) {
    throw std::runtime_error("CSRMAT C=AB test 7 failed.");
  }

  std::vector<int>& cols2 = C2.getGraph().packedColumnIndices;
  std::vector<double>& coefs2 = C2.getPackedCoefs();

  if (cols2.size() != 7) {
    throw std::runtime_error("CSRMAT C=AB test 8 failed.");
  }

  if (cols2[3] != 2) {
    throw std::runtime_error("CSRMAT C=AB test 9 failed.");
  }

  if (std::abs(coefs2[3] - 2.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 10 failed.");
  }

  if (std::abs(coefs2[4] - 3.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 11 failed.");
  }

  if (std::abs(coefs2[6] - 1.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=AB test 12 failed.");
  }

  fa2.clear();

  if (fa2.getNumRows() != 0 || fa2.begin() != fa2.end()) {
    throw std::runtime_error("FillableMat::clear() test failed.");
  }

  fei::multiply_trans_CSRMat_CSRMat(A, B, C);

  if (C.getNumRows() != 3) {
    throw std::runtime_error("CSRMAT C=A^TB test 1 failed.");
  }

  std::vector<int>& tcols = C.getGraph().packedColumnIndices;
  std::vector<double>& tcoefs = C.getPackedCoefs();

  if (tcols.size() != 7) {
    throw std::runtime_error("CSRMAT C=A^TB test 2 failed.");
  }

  if (tcols[2] != 0) {
    throw std::runtime_error("CSRMAT C=A^TB test 3 failed.");
  }

  if (std::abs(tcoefs[3] - 14.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=A^TB test 4 failed.");
  }

  if (std::abs(tcoefs[4] - 9.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=A^TB test 5 failed.");
  }

  if (std::abs(tcoefs[6] - 9.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=A^TB test 6 failed.");
  }
}

TEUCHOS_UNIT_TEST(CSRMatCSVec, multiply_CSRMat_CSRMat_sparse)
{
  fei::FillableMat fa, fb;

  fa.putCoef(0, 0, 1.0);
  fa.putCoef(0, 1, 1.0);
  fa.putCoef(2, 1, 3.0);
  fa.putCoef(2, 2, 3.0);

  fb.putCoef(1, 1, 2.0);
  fb.putCoef(2, 1, 2.0);
  fb.putCoef(3, 2, 3.0);

  fei::CSRMat A(fa), B(fb), C;

  fei::multiply_CSRMat_CSRMat(A, B, C);

  if (C.getNumRows() != 2) {
    throw std::runtime_error("CSRMAT C=AB test s1 failed.");
  }

  std::vector<int>& cols = C.getGraph().packedColumnIndices;
  std::vector<double>& coefs = C.getPackedCoefs();

  if (cols.size() != 2) {
    throw std::runtime_error("CSRMAT C=AB test s2 failed.");
  }

  if (cols[1] != 1) {
    throw std::runtime_error("CSRMAT C=AB test s3 failed.");
  }

  if (std::abs(coefs[0] - 2.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test s4 failed.");
  }

  if (std::abs(coefs[1] - 12.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test s5 failed.");
  }

  fei::multiply_CSRMat_CSRMat(B, A, C);

  if (C.getNumRows() != 1) {
    throw std::runtime_error("CSRMAT C=AB test s2.1 failed.");
  }

  std::vector<int>& cols2 = C.getGraph().packedColumnIndices;
  std::vector<double>& coefs2 = C.getPackedCoefs();

  if (cols2.size() != 2) {
    throw std::runtime_error("CSRMAT C=AB test s2.2 failed.");
  }

  if (cols2[1] != 2) {
    throw std::runtime_error("CSRMAT C=AB test s2.3 failed.");
  }

  if (std::abs(coefs2[0] - 9.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test s2.4 failed.");
  }

  if (std::abs(coefs2[1] - 9.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=AB test s2.5 failed.");
  }

  fei::FillableMat fa2, fb2;

  fa2.putCoef(1, 1, 2.0);
  fa2.putCoef(2, 1, 2.0);
  fa2.putCoef(3, 2, 3.0);

  fb2.putCoef(0, 0, 1.0);
  fb2.putCoef(0, 1, 1.0);
  fb2.putCoef(2, 1, 3.0);
  fb2.putCoef(2, 2, 3.0);

  fei::CSRMat A2(fa2), B2(fb2), C2;

  fei::multiply_trans_CSRMat_CSRMat(A2, B2, C2);

  if (C2.getNumRows() != 1) {
    throw std::runtime_error("CSRMAT C=A^TB test s7 failed.");
  }

  std::vector<int>& cols3 = C2.getGraph().packedColumnIndices;
  std::vector<double>& coefs3 = C2.getPackedCoefs();

  if (cols3.size() != 2) {
    throw std::runtime_error("CSRMAT C=A^TB test s8 failed.");
  }

  if (cols3[1] != 2) {
    throw std::runtime_error("CSRMAT C=A^TB test s9 failed.");
  }

  if (std::abs(coefs3[0] - 6.0) > 1.e-13) {
    throw std::runtime_error("CSRMAT C=A^TB test s10 failed.");
  }

  if (std::abs(coefs3[1] - 6.0) > 1.e-14) {
    throw std::runtime_error("CSRMAT C=A^TB test s11 failed.");
  }

  fa2.clear();

  if (fa2.getNumRows() != 0 || fa2.begin() != fa2.end()) {
    throw std::runtime_error("FillableMat::clear() test failed.");
  }
}

TEUCHOS_UNIT_TEST(CSRMatCSVec, csvec_add_entry)
{
  std::vector<int> ind(3);
  std::vector<double> coef(3);

  ind[0] = 2; ind[1] = 5; ind[2] = 8;
  coef[0] = 2.0; coef[1] = 5.0; coef[2] = 8.0;

  fei::CSVec csv;

  for(int i=ind.size()-1; i>=0; --i) {
    fei::add_entry(csv, ind[i], coef[i]);
  }

  if (csv.indices() != ind) {
    throw std::runtime_error("add_entry(CSVec... failed 1.");
  }

  if (csv.coefs() != coef) {
    throw std::runtime_error("add_entry(CSVec... failed 2.");
  }

  coef[1] = 7.0;
  fei::put_entry(csv, ind[1], 7.0);

  TEUCHOS_TEST_EQUALITY(csv.indices() == ind, true, out, success);
  TEUCHOS_TEST_EQUALITY(csv.coefs() == coef, true, out, success);
}

TEUCHOS_UNIT_TEST(CSRMatCSVec, constructors)
{
  fei::FillableVec fv;

  fv.putEntry(0, 0.0);
  fv.putEntry(1, 1.0);
  fv.addEntry(2, 2.0);

  fei::CSVec csv(fv);

  std::vector<int>& inds = csv.indices();
  std::vector<double>& coefs = csv.coefs();

  if (inds.size() != fv.size()) {
    throw std::runtime_error("CSVec ctor test failed.");
  }

  fei::FillableVec::iterator iter = fv.begin(), iter_end = fv.end();
  unsigned i=0;
  for(; iter != iter_end; ++iter, ++i) {
    TEUCHOS_TEST_EQUALITY(inds[i], iter->first, out, success);
    TEUCHOS_TEST_EQUALITY(coefs[i], iter->second, out, success);
  }

  fei::FillableMat fm;

  fm.sumInCoef(0, 0, 0.0);
  fm.sumInCoef(1, 1, 1.0);
  fm.sumInCoef(2, 2, 2.0);
  fm.sumInCoef(3, 3, 3.0);

  fei::CSRMat csrm(fm);

  if (csrm.getNumRows() != 4) {
    throw std::runtime_error("CSRMat ctor test failed.");
  }

  if (csrm.getGraph().packedColumnIndices.size() != 4) {
    throw std::runtime_error("CSRMat ctor test 2 failed.");
  }

  if (csrm.getPackedCoefs()[2] != 2.0) {
    throw std::runtime_error("CSRMat ctor test 3 failed.");
  }
}

}//namespace <anonymous>
