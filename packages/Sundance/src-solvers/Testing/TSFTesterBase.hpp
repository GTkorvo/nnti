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


#ifndef TSF_TESTERBASE_HPP
#define TSF_TESTERBASE_HPP

#include "TSFLinearOperator.hpp"
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
  class TesterBase : public ObjectWithVerbosity<TesterBase<Scalar> >
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** */
    TesterBase(){;}

    /** */
    virtual ~TesterBase(){;}

    /** */
    virtual bool runAllTests() const = 0 ;


    /** */
    bool checkTest(const TestSpecifier<Scalar>& spec,
                    const ScalarMag& err, 
                    const string& testName) const ;

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;

  };

  template <class Scalar> 
  inline void TesterBase<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Thyra::randomize(Scalar(-ST::one()),Scalar(+ST::one()),x.ptr().get());
    
  }

  template <class Scalar> 
  inline bool TesterBase<Scalar>
  ::checkTest(const TestSpecifier<Scalar>& spec,
              const ScalarMag& err, 
              const string& testName) const 
  {
    bool rtn = true;
    if (err > spec.errorTol())
      {
        cerr << testName << " test FAILED: err=" << err << ", tol = " 
             << spec.errorTol() << endl;
        rtn = false;
      }
    else if (err > spec.warningTol())
      {
        cerr << "WARNING: " << testName << " test err="
             << err << " could not beat tol = " 
             << spec.warningTol() << endl;
      }
    else
      {
        cerr << "test " << testName << " PASSED with tol=" << spec.errorTol() << endl;
      }
    return rtn;
  }
  
}
#endif
