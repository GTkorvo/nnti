// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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
// ***********************************************************************
// @HEADER

#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAdjointLinearOpWithSolve.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"

#include "OperatorSolveHelpers.hpp"


namespace Thyra {


/** \brief Main test driver that tests the LinearOpWithSolveBase interface and
 * supporting software.
 */
template <class Scalar>
bool run_linear_op_with_solve_tests(
  const int n,
  const typename ScalarTraits<Scalar>::magnitudeType maxRelErr,
  const bool showAllTests,
  const bool dumpAll,
  FancyOStream &out
  )
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  using Teuchos::as;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Ordinal Index;

  out
    << "\n***"
    << "\n*** Entering run_linear_op_with_solve_tests<"<<ST::name()<<">(...) ..."
    << "\n***\n";

  bool success = true;

  out << "\nA) Creat a serial vector space of dimension "<<n<<" ...\n";
  const RCP<const Thyra::VectorSpaceBase<Scalar> >
    vs = Thyra::defaultSpmdVectorSpace<Scalar>(n);

  out << "\nB) Create a nonsingular MV object ...\n";
  const RCP<const MultiVectorBase<Scalar> > M =
    createNonsingularMultiVector(vs);

  out << "\nC) Create DefaultSerialDenseLinearOpWithSolve object M_lows from M ...\n";

  RCP<Thyra::LinearOpWithSolveBase<Scalar> >
    M_lows_nonconst = Thyra::linearOpWithSolve<Scalar>(
      *Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>(),
      M );

  RCP<const Thyra::LinearOpWithSolveBase<Scalar> > 
    M_lows = M_lows_nonconst;

  out << "\nD) Test the LOB interface of M_lows ...\n";
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.set_all_error_tol(maxRelErr);
  linearOpTester.set_all_warning_tol(1e-2*maxRelErr);
  linearOpTester.show_all_tests(showAllTests);
  linearOpTester.dump_all(dumpAll);
  {
    OSTab tab(out);
    const bool result = linearOpTester.check(*M_lows, &out);
    if(!result) success = false;
  }

  out << "\nE) Test the LOWSB interface of M_lows ...\n";
  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  {
    RCP<ParameterList> pl = parameterList();
    pl->set("All Solve Tol", maxRelErr);
    pl->set("All Slack Error Tol", 1e+1*maxRelErr);
    pl->set("All Slack Warning Tol", maxRelErr);
    pl->set("Show All Tests", showAllTests);
    pl->set("Dump All", dumpAll);
    linearOpWithSolveTester.setParameterList(pl);
  }
  {
    OSTab tab(out);
    const bool result = linearOpWithSolveTester.check(*M_lows, &out);
    if(!result) success = false;
  }

  out << "\nF) Create a DefaultInverseLinearOp object invM from M_lows and test the LOB interface ...\n";
  RCP<const Thyra::LinearOpBase<Scalar> > invM = inverse(M_lows);
  {
    OSTab tab(out);
    const bool result = linearOpTester.check(*invM, &out);
    if(!result) success = false;
  }

  out << "\nG) Test DefaultAdjointLinearOpWithSolve ...\n";
  {
    OSTab tab(out);

    out << "\nG.1) Create and test DefaultAdjointLinearOpWithSolve object M_lows_adj wrapping const M_lows ...\n";
    RCP<const Thyra::LinearOpWithSolveBase<Scalar> >
      M_lows_adj = adjointLows(M_lows);

    {
      OSTab tab2(out);
      
      out << "\nG.1.a) Test that we can extract the underlying const M_lows ...\n";
      {
        OSTab tab3(out);
        TEST_EQUALITY( M_lows,
          rcp_dynamic_cast<const Thyra::DefaultAdjointLinearOpWithSolve<Scalar> >(M_lows_adj,true)->getOp() );
      }
      
      out << "\nG.1.b) Testing LOB interface of DefaultAdjointLinearOpWithSolve object M_lows_adj ...\n";
      {
        OSTab tab3(out);
        const bool result = linearOpTester.check(*M_lows_adj, &out);
        if(!result) success = false;
      }
      
      out << "\nG.1.c) Testing LOWSB interface of DefaultAdjointLinearOpWithSolve object M_lows_adj ...\n";
      {
        OSTab tab3(out);
        const bool result = linearOpWithSolveTester.check(*M_lows_adj, &out);
        if(!result) success = false;
      }
      
      out << "\nG.1.d) Testing that M_lows_adj is the adjoint of M (M_adj) ...\n";
      const RCP<const LinearOpBase<Scalar> > M_adj = Thyra::adjoint<Scalar>(M);
      {
        OSTab tab3(out);
        const bool result = linearOpTester.compare(*M_lows_adj, *M_adj, &out);
        if(!result) success = false;
      }

    }

    out << "\nG.2) Create and test DefaultAdjointLinearOpWithSolve object M_lows_adj_nonconst wrapping non-const M_lows ...\n";
    RCP<Thyra::LinearOpWithSolveBase<Scalar> >
      M_lows_adj_nonconst = nonconstAdjointLows<Scalar>(M_lows_nonconst);

    {
      OSTab tab3(out);
      
      out << "\nG.2.a) Test that we can extract the underlying non-const and const M_lows ...\n";
      {
        OSTab tab4(out);
        TEST_EQUALITY( M_lows,
          rcp_dynamic_cast<Thyra::DefaultAdjointLinearOpWithSolve<Scalar> >(M_lows_adj_nonconst,true)->getOp() );
        TEST_EQUALITY( M_lows_nonconst,
          rcp_dynamic_cast<Thyra::DefaultAdjointLinearOpWithSolve<Scalar> >(M_lows_adj_nonconst,true)->getNonconstOp() );
      }
      
      out << "\nG.2.b) Only testing LOB interface of DefaultAdjointLinearOpWithSolve object M_lows_adj_nonconst ...\n";
      {
        OSTab tab4(out);
        const bool result = linearOpTester.check(*M_lows_adj_nonconst, &out);
        if(!result) success = false;
      }
    }
    
  }
  
  return success;

}


} // namespace Thyra


int main( int argc, char* argv[] )
{

  using Teuchos::CommandLineProcessor;
  using Teuchos::ScalarTraits;
  using Teuchos::as;

  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from the command-line
    //


    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int n = 4;
    clp.setOption( "n", &n, "Size of the system." );

    double epsScale = 2e+2;
    clp.setOption( "eps-scale", &epsScale,
      "Constant (greater than 1) to scale eps by in error tests." );

    bool showAllTests = false;
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests,
      "Determines if detailed tests are shown or not." );

    bool dumpAll = false;
    clp.setOption( "dump-all", "no-dump", &dumpAll,
      "Determines if quantities are dumped or not." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Run the tests
    //

#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !Thyra::run_linear_op_with_solve_tests<float>(
          n, as<float>(epsScale*ScalarTraits<float>::eps()), showAllTests, dumpAll, *out)
      ) success = false;
#endif
    if( !Thyra::run_linear_op_with_solve_tests<double>(
          n, as<double>(epsScale*ScalarTraits<double>::eps()), showAllTests, dumpAll, *out)
      ) success = false;
#if defined(HAVE_THYRA_COMPLEX) && defined(HAVE_THYRA_TEUCHOS_BLASFLOAT)
    if( !Thyra::run_linear_op_with_solve_tests<std::complex<float> >(
          n, as<float>(epsScale*ScalarTraits<float>::eps()), showAllTests, dumpAll, *out)
      ) success = false;
#endif
#if defined(HAVE_THYRA_COMPLEX)
    if( !Thyra::run_linear_op_with_solve_tests<std::complex<double> >(
          n, as<double>(epsScale*ScalarTraits<double>::eps()), showAllTests, dumpAll, *out)
      ) success = false;
#endif

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success);

  if(success)
    *out << "\nAll of the tests seem to have run successfully!\n";
  else
    *out << "\nOh no! at least one of the test failed!\n";	
  
  return success ? 0 : 1;

}
