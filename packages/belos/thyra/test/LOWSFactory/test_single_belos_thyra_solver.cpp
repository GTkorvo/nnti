#include "test_single_belos_thyra_solver.hpp"

#ifndef __sun

#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_BELOS_IFPACK
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif

#endif // __sun

bool Thyra::test_single_belos_thyra_solver(
  const std::string                       matrixFile
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *solveParamList
  ,Teuchos::FancyOStream                  *out_arg
  )
{
  using Teuchos::rcp;
  using Teuchos::OSTab;
  bool result, success = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = Teuchos::rcp(out_arg,false);

  try {

#ifndef __sun

    if(out.get()) {
      *out << "\n***"
           << "\n*** Testing Thyra::BelosLinearOpWithSolveFactory (and Thyra::BelosLinearOpWithSolve)"
           << "\n***\n"
           << "\nEchoing input options:"
           << "\n  matrixFile             = " << matrixFile
           << "\n  testTranspose          = " << testTranspose
           << "\n  numRandomVectors       = " << numRandomVectors
           << "\n  maxFwdError            = " << maxFwdError
           << "\n  maxResid               = " << maxResid
           << "\n  showAllTests           = " << showAllTests
           << "\n  dumpAll                = " << dumpAll
           << std::endl;
    }

    if(out.get()) *out << "\nA) Reading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";
  
    Epetra_SerialComm comm;
    Teuchos::RefCountPtr<Epetra_CrsMatrix> epetra_A;
    EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

    Teuchos::RefCountPtr<LinearOpBase<double> > A = Teuchos::rcp(new EpetraLinearOp(epetra_A));

    if(out.get() && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A,Teuchos::VERB_EXTREME);

    if(out.get()) *out << "\nB) Creating a BelosLinearOpWithSolveFactory object opFactory ...\n";

    Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
      opFactory;
    if(1) {
      Teuchos::RefCountPtr<BelosLinearOpWithSolveFactory<double> >
        belosOpFactory = Teuchos::rcp(new BelosLinearOpWithSolveFactory<double>());
      opFactory = belosOpFactory;
    }
    opFactory->setParameterList(Teuchos::rcp(solveParamList,false));

    if(out.get()) *out << "\nC) Creating a BelosLinearOpWithSolve object nsA from A ...\n";

    Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
      nsA = opFactory->createOp();

    opFactory->initializeOp( A, &*nsA );

    if(out.get()) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

    LinearOpTester<double> linearOpTester;
    linearOpTester.check_adjoint(testTranspose);
    linearOpTester.num_random_vectors(numRandomVectors);
    linearOpTester.set_all_error_tol(maxFwdError);
    linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
    linearOpTester.show_all_tests(showAllTests);
    linearOpTester.dump_all(dumpAll);
    Thyra::seed_randomize<double>(0);
    result = linearOpTester.check(*nsA,OSTab(out).getOStream().get());
    if(!result) success = false;

    if(out.get()) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
    LinearOpWithSolveTester<double> linearOpWithSolveTester;
    linearOpWithSolveTester.turn_off_all_tests();
    linearOpWithSolveTester.check_forward_default(true);
    linearOpWithSolveTester.check_forward_residual(true);
    if(testTranspose) {
      linearOpWithSolveTester.check_adjoint_default(true);
      linearOpWithSolveTester.check_adjoint_residual(true);
    }
    else {
      linearOpWithSolveTester.check_adjoint_default(false);
      linearOpWithSolveTester.check_adjoint_residual(false);
    }
    linearOpWithSolveTester.set_all_solve_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_error_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_warning_tol(1e+1*maxResid);
    linearOpWithSolveTester.forward_default_residual_error_tol(2*maxResid);
    linearOpWithSolveTester.forward_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.adjoint_default_residual_error_tol(2*maxResid);
    linearOpWithSolveTester.adjoint_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.show_all_tests(showAllTests);
    linearOpWithSolveTester.dump_all(dumpAll);
    Thyra::seed_randomize<double>(0);
    result = linearOpWithSolveTester.check(*nsA,OSTab(out).getOStream().get());
    if(!result) success = false;

    
#else // __sun
		
		if(out.get()) *out << "\nTest failed since is was not even compiled since __sun was defined!\n";
		success = false;

#endif // __sun

  }
	catch( const std::exception &excpt ) {
    if(out.get()) *out << std::flush;
		std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
		success = false;
	}
   
  return success;
    
}
