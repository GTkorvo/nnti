/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
 /* @HEADER@ */


#ifndef TSF_MATRIXMATRIXTESTER_HPP
#define TSF_MATRIXMATRIXTESTER_HPP

#include "TSFLinearOperatorDecl.hpp"
#include "TSFEpetraMatrixMatrixProduct.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "TSFLinearCombinationImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFSimpleComposedOpImpl.hpp"
#include "TSFSimpleScaledOpImpl.hpp"
#include "TSFSimpleAddedOpImpl.hpp"
#include "TSFSimpleDiagonalOpImpl.hpp"
#include "TSFRandomSparseMatrixBuilderImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;
using Thyra::TestSpecifier;

namespace TSFExtended
{

/** */
template <class Scalar>
class MatrixMatrixTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  MatrixMatrixTester(const LinearOperator<Scalar>& A,
    const LinearOperator<Scalar>& B,
    const TestSpecifier<Scalar>& prodSpec,
    const TestSpecifier<Scalar>& diagProdSpec);

  /** */
  bool runAllTests() const ;

  /** */
  bool prodTest() const ;

  /** */
  bool diagProdTest() const ;


private:

  LinearOperator<Scalar> A_;

  LinearOperator<Scalar> B_;

  TestSpecifier<Scalar> prodSpec_;

  TestSpecifier<Scalar> diagProdSpec_;

};

template <class Scalar> 
inline MatrixMatrixTester<Scalar>
::MatrixMatrixTester(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B,
  const TestSpecifier<Scalar>& prodSpec,
  const TestSpecifier<Scalar>& diagProdSpec)
  : TesterBase<Scalar>(), 
    A_(A),
    B_(B),
    prodSpec_(prodSpec),
    diagProdSpec_(diagProdSpec)
{;}

template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = prodTest() && pass;
  pass = diagProdTest() && pass;

  return pass;
}

template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::prodTest() const 
{
  if (prodSpec_.doTest())
  {
    Out::root() << "running matrix-matrix multiply test..." << endl;
    LinearOperator<Scalar> composed = A_ * B_;
    LinearOperator<Scalar> multiplied = epetraMatrixMatrixProduct(A_, B_);

    Vector<Scalar> x = B_.domain().createMember();
    randomizeVec(x);
    Out::root() << "computing implicit product y1 = (A*B)*x..." << endl;
    Vector<Scalar> y1 = composed*x;
    Out::root() << "computing explicit product y2 = (A*B)*x..." << endl;
    Vector<Scalar> y2 = multiplied*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << endl;
    return checkTest(prodSpec_, err, "matrix-matrix multiply");
  }
  Out::root() << "skipping matrix-matrix multiply test..." << endl;
  return true;
}


template <class Scalar> 
inline bool MatrixMatrixTester<Scalar>
::diagProdTest() const 
{
  if (diagProdSpec_.doTest())
  {
    Out::root() << "running diagonal matrix multiplication test..." << endl;

    Vector<Scalar> x = A_.domain().createMember();
    randomizeVec(x);

    Vector<Scalar> d = A_.range().createMember();
    randomizeVec(d);
        
    LinearOperator<Scalar> D = diagonalOperator(d);
    LinearOperator<Scalar> DA = epetraLeftScale(d, A_);

    Out::root() << "computing implicit y1 = D*A*x..." << endl;
    Vector<Scalar> y1 = D*A_*x;
    Out::root() << "computing explicit y2 = D*A*x..." << endl;
    Vector<Scalar> y2 = DA*x;

    ScalarMag err = (y1 - y2).norm2();

    Out::root() << "|y1-y2| = " << err << endl;
    return checkTest(diagProdSpec_, err, "diagonal matrix multiplication");
  }
  Out::root() << "skipping diagonal matrix-matrix test..." << endl;
  return true;
}

  
  
}
#endif
