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


#ifndef TSF_COMPOUNDTESTER_HPP
#define TSF_COMPOUNDTESTER_HPP

#include "TSFLinearOperator.hpp"
#include "TSFScaledOperator.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;
using Thyra::TestSpecifier;

namespace TSFExtended
{

  /** */
  template <class Scalar>
  class CompoundTester
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** */
    CompoundTester(const LinearOperator<Scalar>& A,
                   const LinearOperator<Scalar>& B,
                   const TestSpecifier<Scalar>& sumSpec,
                   const TestSpecifier<Scalar>& composedSpec,
                   const TestSpecifier<Scalar>& scaledSpec);

    /** */
    bool runAllTests() const ;

    /** */
    bool sumTest() const ;

    /** */
    bool composedTest() const ;

    /** */
    bool scaledTest() const ;

  private:

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;
    
    LinearOperator<Scalar> A_;

    LinearOperator<Scalar> B_;

    TestSpecifier<Scalar> sumSpec_;

    TestSpecifier<Scalar> composedSpec_;

    TestSpecifier<Scalar> scaledSpec_;

  };

  template <class Scalar> 
  inline CompoundTester<Scalar>
  ::CompoundTester(const LinearOperator<Scalar>& A,
                   const LinearOperator<Scalar>& B,
                   const TestSpecifier<Scalar>& sumSpec,
                   const TestSpecifier<Scalar>& composedSpec,
                   const TestSpecifier<Scalar>& scaledSpec)
    : A_(A),
      B_(B),
      sumSpec_(sumSpec),
      composedSpec_(composedSpec),
      scaledSpec_(scaledSpec)
  {;}

  template <class Scalar> 
  inline bool CompoundTester<Scalar>
  ::runAllTests() const
  {
    bool pass = true;

    pass = sumTest() && pass;
    pass = composedTest() && pass;
    pass = scaledTest() && pass;

    return pass;
  }

  template <class Scalar> 
  inline void CompoundTester<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Thyra::randomize(Scalar(-ST::one()),Scalar(+ST::one()),x.ptr().get());
    
  }

  template <class Scalar> 
  inline bool CompoundTester<Scalar>
  ::sumTest() const 
  {
    if (sumSpec_.doTest())
      {
        cerr << "running operator addition test..." << endl;
        LinearOperator<Scalar> sum = A_ + B_;

        Vector<Scalar> x = A_.domain().createMember();
        randomizeVec(x);
        cerr << "computing y1 = (A+B)*x..." << endl;
        Vector<Scalar> y1 = sum*x;
        cerr << "computing y2 = A*x + B*x..." << endl;
        Vector<Scalar> y2 = A_*x + B_*x;
    
        ScalarMag err = (y1 - y2).norm2();

        cerr << "|y1-y2| = " << err << endl;
        if (err > sumSpec_.errorTol())
          {
            cerr << "operator addition test FAILED: tol = " 
                 << sumSpec_.errorTol() << endl;
            return false;
          }
        else if (err > sumSpec_.warningTol())
          {
            cerr << "WARNING: operator addition test could not beat tol = " 
                 << sumSpec_.warningTol() << endl;
          }
      }
    else
      {
        cerr << "skipping operator addition test..." << endl;
      }
    cerr << "operator addition test PASSED: tol = " 
         << sumSpec_.errorTol() << endl;
    return true;
  }

  template <class Scalar> 
  inline bool CompoundTester<Scalar>
  ::composedTest() const 
  {
    if (composedSpec_.doTest())
      {
        cerr << "running operator composition test..." << endl;
        LinearOperator<Scalar> composed = A_ * B_;

        Vector<Scalar> x = B_.domain().createMember();
        randomizeVec(x);
        cerr << "computing y1 = (A*B)*x..." << endl;
        Vector<Scalar> y1 = composed*x;
        cerr << "computing y2 = B*x..." << endl;
        Vector<Scalar> y2 = B_*x;
        cerr << "computing y3 = A*y2..." << endl;
        Vector<Scalar> y3 = A_*y2;

        ScalarMag err = (y1 - y3).norm2();

        cerr << "|y1-y3| = " << err << endl;
        if (err > composedSpec_.errorTol())
          {
            cerr << "operator composition test FAILED: tol = " 
                 << composedSpec_.errorTol() << endl;
            return false;
          }
        else if (err > composedSpec_.warningTol())
          {
            cerr << "WARNING: operator composition test could not beat tol = " 
                 << composedSpec_.warningTol() << endl;
          }
      }
    else
      {
        cerr << "skipping operator composition test..." << endl;
      }
    cerr << "operator composition test PASSED: tol = " 
         << composedSpec_.errorTol() << endl;
    return true;
  }


  template <class Scalar> 
  inline bool CompoundTester<Scalar>
  ::scaledTest() const 
  {
    if (scaledSpec_.doTest())
      {
        cerr << "running operator scaling test..." << endl;
        Scalar alpha = sqrt(2.0);
        LinearOperator<Scalar> scaled = new ScaledOperator<Scalar>(A_, alpha);

        Vector<Scalar> x = A_.domain().createMember();
        randomizeVec(x);
        cerr << "computing y1 = (alpha*A)*x..." << endl;
        Vector<Scalar> y1 = scaled*x;
        cerr << "computing y2 = A*x..." << endl;
        Vector<Scalar> y2 = A_*x;
        cerr << "computing y3 = alpha*y2..." << endl;
        Vector<Scalar> y3 = alpha*y2;

        ScalarMag err = (y1 - y3).norm2();

        cerr << "|y1-y3| = " << err << endl;
        if (err > scaledSpec_.errorTol())
          {
            cerr << "operator scaling test FAILED: tol = " 
                 << scaledSpec_.errorTol() << endl;
            return false;
          }
        else if (err > scaledSpec_.warningTol())
          {
            cerr << "WARNING: operator scaling test could not beat tol = " 
                 << scaledSpec_.warningTol() << endl;
          }
      }
    else
      {
        cerr << "skipping operator scaling test..." << endl;
      }
    cerr << "operator scaling test PASSED: tol = " 
         << scaledSpec_.errorTol() << endl;
    return true;
  }

  
  
}
#endif
