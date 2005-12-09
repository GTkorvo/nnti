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


#ifndef TSF_VECTORTESTER_HPP
#define TSF_VECTORTESTER_HPP

#include "TSFVector.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;
using Thyra::TestSpecifier;

namespace TSFExtended
{

  /** 
   * Run comparisons between element-wise calculations and Vector member
   * functions.
   */
  template <class Scalar>
  class VectorTester
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** */
    VectorTester(const VectorSpace<Scalar>& space,
                 const TestSpecifier<Scalar>& spec);

    /** */
    bool runAllTests() const ;

    /** */
    bool sumTest() const ;

    /** */
    bool setElementTest() const ;

    /** */
    bool dotStarTest() const ;

    /** */
    bool dotSlashTest() const ;

    /** */
    bool scalarMultTest() const ;

    /** */
    bool overloadedUpdateTest() const ;

  private:

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;

    TestSpecifier<Scalar> spec_;

    VectorSpace<Scalar> space_;

  };

  template <class Scalar> 
  inline VectorTester<Scalar>
  ::VectorTester(const VectorSpace<Scalar>& space,
                 const TestSpecifier<Scalar>& spec)
    : spec_(spec), space_(space)
  {;}

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::runAllTests() const
  {
    bool pass = true;

    pass = sumTest() && pass;
    pass = setElementTest() && pass;
    pass = dotStarTest() && pass;
    pass = dotSlashTest() && pass;
    pass = scalarMultTest() && pass;
    pass = overloadedUpdateTest() && pass;

    return pass;
  }

  template <class Scalar> 
  inline void VectorTester<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Thyra::randomize(Scalar(-ST::one()),Scalar(+ST::one()),x.ptr().get());
    
  }

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::sumTest() const 
  {
    if (spec_.doTest())
      {
        cerr << "running vector addition test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);

        /* do the operation with member functions */
        x = a + b ;

        /* do the operation with member functions */

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            double b_i = b.getElement(i);
            y.setElement(i, a_i + b_i );
          }
	
        double err = (x-y).normInf();

        cerr << "|sum error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector sum test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector sum test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector addition test..." << endl;
      }
    cerr << "vector addition test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::setElementTest() const 
  {
    if (spec_.doTest())
      {
        cerr << "running setElement test..." << endl;

        Vector<Scalar> a = space_.createMember();
	
        /* we will load a vector with a_i = i, and then do
         * the sum of all elements. If done correctly, the sum will equal 
         * N*(N+1)*(2N+1)/6.
         */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            a.setElement(i, i);
          }
        cerr << "a = " << endl << a << endl;
        Vector<double> b = a.copy();
        b = b.dotStar(a);

        double sum = 0.0;
        for (int i=low; i<high; i++)
          {
            cerr << i << " " << a.getElement(i) << " " << i*a.getElement(i)
                 << endl;
            sum += i * a.getElement(i);
          }
	
        double thyraSum = Thyra::sum(*(b.ptr()));
        cerr << "thyra sum = " << sum << endl;
        cerr << "thyra sum = " << thyraSum << endl;

        double err = ::fabs(sum - thyraSum);

        cerr << "|setElement error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector setElement test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector setElement test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector setElement test..." << endl;
      }
    cerr << "vector setElement test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::dotStarTest() const 
  {
    if (spec_.doTest())
      {
        cerr << "running vector dotStar test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = a.dotStar(b);

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            double b_i = b.getElement(i);
            y.setElement(i, a_i * b_i);
          }
	
        double err = (x-y).normInf();

        cerr << "|dotStar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector dotStar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector dotStar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector dotStar test..." << endl;
      }
    cerr << "vector dotStar test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }


  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::dotSlashTest() const 
  {
    if (spec_.doTest())
      {
        cerr << "running vector dotSlash test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = a.dotSlash(b);

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            double b_i = b.getElement(i);
            y.setElement(i, a_i / b_i);
          }
	
        double err = (x-y).normInf();

        cerr << "|dotSlash error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector dotSlash test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector dotSlash test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector dotSlash test..." << endl;
      }
    cerr << "vector dotSlash test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  
  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::scalarMultTest() const 
  {
    if (spec_.doTest())
      {
        cerr << "running vector scalarMult test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);


        /* do the operation with member functions */
        x = 3.14*a;

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            y.setElement(i, 3.14*a_i);
          }
	
        double err = (x-y).normInf();

        cerr << "|scalarMult error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector scalarMult test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector scalarMult test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector scalarMult test..." << endl;
      }
    cerr << "vector scalarMult test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }
 
  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::overloadedUpdateTest() const 
  {
    if (spec_.doTest())
      {
        cerr << "running vector overloadedUpdate test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = 3.14*a + 1.4*b;

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            double b_i = b.getElement(i);
            y.setElement(i, 3.14*a_i + 1.4*b_i);
          }
	
        double err = (x-y).normInf();

        cerr << "|overloadedUpdate error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector overloadedUpdate test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector overloadedUpdate test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector overloadedUpdate test..." << endl;
      }
    cerr << "vector overloadedUpdate test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  
  
}
#endif
