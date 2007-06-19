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
#include "TSFProductVectorSpace.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Thyra_SUNDIALS_Ops.hpp"
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
                 const TestSpecifier<Scalar>& spec,
                 const Teuchos::MPIComm& comm = Teuchos::MPIComm::world());

    /** */
    bool runAllTests() const ;

    /** */
    bool sumTest() const ;

    /** */
    bool dotStarTest() const ;

    /** */
    bool dotSlashTest() const ;

    /** */
    bool scalarMultTest() const ;

    /** */
    bool overloadedUpdateTest() const ;

    /** */
    bool reciprocalTest() const ;

    /** */
    bool minQuotientTest() const ;

    /** */
    bool addScalarTest() const ;

    /** */
    bool compareToScalarTest() const ;

    /** */
    bool constraintMaskTest() const ;

    /** */
    bool indexTest() const ;

  private:

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;

    TestSpecifier<Scalar> spec_;

    VectorSpace<Scalar> space_;

    Teuchos::MPIComm comm_;

  };

  template <class Scalar> 
  inline VectorTester<Scalar>
  ::VectorTester(const VectorSpace<Scalar>& space,
                 const TestSpecifier<Scalar>& spec,
                 const Teuchos::MPIComm& comm)
    : spec_(spec), space_(space), comm_(comm)
  {;}

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::runAllTests() const
  {
    bool pass = true;

    pass = sumTest() && pass;
    pass = dotStarTest() && pass;
    pass = dotSlashTest() && pass;
    pass = scalarMultTest() && pass;
    pass = overloadedUpdateTest() && pass;
    pass = reciprocalTest() && pass;
    pass = minQuotientTest() && pass;
    pass = constraintMaskTest() && pass;
    pass = compareToScalarTest() && pass;
    pass = indexTest() && pass;

    return pass;
  }

  template <class Scalar> 
  inline void VectorTester<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;

    /* do the operation elementwise */
    SequentialIterator<Scalar> i;
    for (i=space_.begin(); i != space_.end(); i++)
      {
        x[i] = 2.0*(drand48()-0.5);
      }    
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
        x.zero();
        y.zero();
        randomizeVec(a);
        randomizeVec(b);

        /* do the operation elementwise */
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = a[i] + b[i];
          }

        /* do the operation with member functions */
        x = a + b ;

        cout << "x=" << endl << x << endl;
        cout << "y=" << endl << y << endl;
	
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
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = a[i] * b[i];
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
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = a[i] / b[i];
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
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = 3.14 * a[i];
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
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = 3.14*a[i] + 1.4*b[i];
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

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::reciprocalTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        cerr << "running vector reciprocal test..." << endl;

        Vector<Scalar> a = space_.createMember();
        randomizeVec(a);

        Vector<Scalar> y = space_.createMember();

        /* load the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        int denomsAreOK = true;
	for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            if (a_i != Teuchos::ScalarTraits<Scalar>::zero()) 
              {
                y.setElement(i, 1.0/a_i );
              }
            else
              {
                denomsAreOK=false;
                y.setElement(i, a_i);
              }
          }
        
        Vector<Scalar> x = space_.createMember();
        int tDenomsAreOK = Thyra::VInvTest(*(a.ptr()), x.ptr().get());
        double err = (x - y).norm2();

#ifdef HAVE_MPI
        int localDenomsAreOK = denomsAreOK;
        MPI_Allreduce( (void*) &localDenomsAreOK, (void*) &denomsAreOK, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
#endif

        cerr << "|reciprocal error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector reciprocal test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (tDenomsAreOK != denomsAreOK)
          {
            cerr << "vector reciprocal test FAILED: trilinosDenomsOK="
                 << tDenomsAreOK << ", denomsOK=" << denomsAreOK << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector reciprocal test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        cerr << "skipping vector reciprocal test..." << endl;
      }
    cerr << "vector reciprocal test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::minQuotientTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        cerr << "running vector minQuotient test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();

        randomizeVec(a);
        randomizeVec(b);

        /* perform the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        double minQLocal = Teuchos::ScalarTraits<Scalar>::rmax();
	for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            double b_i = b.getElement(i);
            if (b_i != Teuchos::ScalarTraits<Scalar>::zero())
              {
                double q = a_i/b_i;
                if (q < minQLocal) minQLocal = q;
              }
          }

        double minQ = minQLocal;
        comm_.allReduce((void*) &minQLocal, (void*) &minQ, 1, Teuchos::MPIComm::DOUBLE,
                        Teuchos::MPIComm::MIN);
	

        double tMinQ = Thyra::VMinQuotient(*(a.ptr()), *(b.ptr()));
        cerr << "trilinos minQ = " << tMinQ << endl;
        cerr << "elemwise minQ = " << minQ << endl;
        double err = fabs(minQ - tMinQ);
        
        cerr << "min quotient error=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "min quotient test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: min quotient test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        cerr << "skipping min quotient test..." << endl;
      }
    cerr << "min quotient test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }



  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::constraintMaskTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        cerr << "running vector constraintMask test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> c = space_.createMember();
        randomizeVec(a);

        Vector<Scalar> y = space_.createMember();
        Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

        /* load the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        int allFeasible = true;
	for (int i=low; i<high; i++)
          {
            int feasible = true;
            double a_i = a.getElement(i);
            switch(i%4)
              {
              case 0:
                c.setElement(i, -2.0);
                feasible = a_i < zero;
                break;
              case 1:
                c.setElement(i, -1.0);
                feasible = a_i <= zero;
                break;
              case 2:
                c.setElement(i, 1.0);
                feasible = a_i > zero;
                break;
              case 3:
                c.setElement(i, 2.0);
                feasible = a_i >= zero;
                break;
              default:
                TEST_FOR_EXCEPTION(true, logic_error, "impossible!");
              }
            y.setElement(i, (Scalar) !feasible);
            allFeasible = allFeasible && feasible;
          }
	
        Vector<Scalar> m = space_.createMember();
        int tAllFeasible = Thyra::VConstrMask(*(a.ptr()), *(c.ptr()), m.ptr().get());
        double err = (m - y).norm2();

#ifdef HAVE_MPI
        int localAllFeas = allFeasible;
        cerr << "local all feas=" << localAllFeas << endl;
        MPI_Allreduce( (void*) &localAllFeas, (void*) &allFeasible, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
        cerr << "globalal all feas=" << allFeasible << endl;
#endif

        cerr << "|constraintMask error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector constraintMask test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (allFeasible != tAllFeasible)
          {
            cerr << "vector constraintMask test FAILED: trilinosFeas="
                 << tAllFeasible << ", feas=" << allFeasible << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector constraintMask test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        cerr << "skipping vector constraintMask test..." << endl;
      }
    cerr << "vector constraintMask test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::compareToScalarTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        cerr << "running vector compare-to-scalar test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);

        /* do the operation with member functions */
        Scalar s = 0.5;
        Thyra::VCompare(s, *(a.ptr()), x.ptr().get());

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

	for (int i=low; i<high; i++)
          {
            double a_i = a.getElement(i);
            y.setElement(i, fabs(a_i) >= s );
          }
	
        double err = (x-y).normInf();

        cerr << "|compare-to-scalar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector compare-to-scalar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector compare-to-scalar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector compare-to-scalar test..." << endl;
      }
    cerr << "vector compare-to-scalar test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }


  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::indexTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        cerr << "running vector index test..." << endl;
        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        /* do the operation with member functions */
        Scalar s = 0.5;
        Thyra::VCompare(s, *(a.ptr()), x.ptr().get());

        /* do the operation elementwise */
        int low = space_.lowestLocallyOwnedIndex();
        int high = low + space_.numLocalElements();

        for (int i=low; i<high; i++)
          {
            //double a_i = a.getElement(i);
	    double a_i = a[i];
	    y[i] =  fabs(a_i) >= s;
            //y.setElement(i, fabs(a_i) >= s );
          }
	
        double err = (x-y).normInf();

        cerr << "|index error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cerr << "vector index test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cerr << "WARNING: vector index test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cerr << "skipping vector index test..." << endl;
      }
    cerr << "vector index test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  


}
#endif
