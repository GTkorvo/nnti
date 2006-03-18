// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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

#include "test_single_amesos_thyra_solver.hpp"

#ifndef __sun

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"

#endif // __sun

bool Thyra::test_single_amesos_thyra_solver(
  const std::string                       matrixFile
  ,Teuchos::ParameterList                 *amesosLOWSFPL
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxError
  ,const double                           maxResid
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::FancyOStream                  *out_arg
  )
{
  using Teuchos::OSTab;

  bool result, success = true;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::rcp(out_arg,false);

#ifndef __sun

  if(out.get()) {
    *out
      << "\n***"
      << "\n*** Testing Thyra::AmesosLinearOpWithSolveFactory (and Thyra::AmesosLinearOpWithSolve)"
      << "\n***\n"
      << "\nEchoing input options:"
      << "\n  matrixFile             = " << matrixFile
      << "\n  testTranspose          = " << testTranspose
      << "\n  numRandomVectors       = " << numRandomVectors
      << "\n  maxFwdError            = " << maxFwdError
      << "\n  maxError               = " << maxError
      << "\n  maxResid               = " << maxResid
      << "\n  showAllTests           = " << showAllTests
      << "\n  dumpAll                = " << dumpAll
      << std::endl;
    if(amesosLOWSFPL) {
      OSTab tab(out);
      *out
        << "amesosLOWSFPL:\n";
      amesosLOWSFPL->print(*OSTab(out).getOStream(),0,true);
    }
  }
  
  if(out.get()) *out << "\nA) Reading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";
  
  Epetra_SerialComm comm;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

  Teuchos::RefCountPtr<LinearOpBase<double> >
    A = Teuchos::rcp(new EpetraLinearOp(epetra_A));

  if(out.get() && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A,Teuchos::VERB_EXTREME);

  if(out.get()) *out << "\nB) Creating a AmesosLinearOpWithSolveFactory object opFactory ...\n";

  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
    opFactory = Teuchos::rcp(new AmesosLinearOpWithSolveFactory());

  opFactory->setParameterList(Teuchos::rcp(amesosLOWSFPL,false));

  if(out.get()) *out << "\nC) Creating a AmesosLinearOpWithSolve object nsA ...\n";

  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
    nsA = opFactory->createOp();

  opFactory->initializeOp( A, &*nsA );

  if(out.get()) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

  Thyra::seed_randomize<double>(0);

  LinearOpTester<double> linearOpTester;
  linearOpTester.check_adjoint(testTranspose);
  linearOpTester.num_random_vectors(numRandomVectors);
  linearOpTester.set_all_error_tol(maxFwdError);
  linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
  linearOpTester.show_all_tests(showAllTests);
  linearOpTester.dump_all(dumpAll);
  result = linearOpTester.check(*nsA,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  LinearOpWithSolveTester<double> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_default(true);
  linearOpWithSolveTester.forward_default_residual_error_tol(1.1*maxResid);
  linearOpWithSolveTester.forward_default_residual_warning_tol(2.0*maxResid);
  linearOpWithSolveTester.check_forward_residual(true);
  linearOpWithSolveTester.forward_residual_solve_tol(maxResid);
  linearOpWithSolveTester.forward_residual_slack_error_tol(1e-1*maxResid);
  linearOpWithSolveTester.forward_residual_slack_warning_tol(maxResid);
  linearOpWithSolveTester.check_forward_solution_error(true);
  linearOpWithSolveTester.forward_solution_error_solve_tol(maxError);
  linearOpWithSolveTester.forward_solution_error_slack_error_tol(1e-1*maxError);
  linearOpWithSolveTester.forward_solution_error_slack_warning_tol(maxError);
  if(testTranspose) {
    linearOpWithSolveTester.check_adjoint_default(true);
    linearOpWithSolveTester.adjoint_default_residual_error_tol(1.1*maxResid);
    linearOpWithSolveTester.adjoint_default_residual_warning_tol(2.0*maxResid);
    linearOpWithSolveTester.check_adjoint_residual(true);
    linearOpWithSolveTester.adjoint_residual_solve_tol(maxResid);
    linearOpWithSolveTester.adjoint_residual_slack_error_tol(1e-1*maxResid);
    linearOpWithSolveTester.adjoint_residual_slack_warning_tol(maxResid);
    linearOpWithSolveTester.check_adjoint_solution_error(true);
    linearOpWithSolveTester.adjoint_solution_error_solve_tol(maxError);
    linearOpWithSolveTester.adjoint_solution_error_slack_error_tol(1e-1*maxError);
    linearOpWithSolveTester.adjoint_solution_error_slack_warning_tol(maxError);
  }
  linearOpWithSolveTester.num_random_vectors(numRandomVectors);
  linearOpWithSolveTester.show_all_tests(showAllTests);
  linearOpWithSolveTester.dump_all(dumpAll);
  result = linearOpWithSolveTester.check(*nsA,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nF) Uninitialize the matrix object nsA, scale the epetra_A object by 2.5, and then refactor nsA with epetra_A ...\n";

  opFactory->uninitializeOp(&*nsA); // Optional call but a good idea if changing the operator
  epetra_A->Scale(2.5);
  opFactory->initializeOp(A,&*nsA);
  
  if(out.get()) *out << "\nG) Testing the LinearOpBase interface of nsA ...\n";

  Thyra::seed_randomize<double>(0);

  result = linearOpTester.check(*nsA,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nH) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  result = linearOpWithSolveTester.check(*nsA,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nI) Uninitialize the matrix object nsA, create a scaled (by 2.5) copy  epetra_A2 of epetra_A, and then refactor nsA with epetra_A2 ...\n";

  Teuchos::RefCountPtr<Epetra_CrsMatrix>
    epetra_A2 = Teuchos::rcp(new Epetra_CrsMatrix(*epetra_A));
  epetra_A2->Scale(2.5);
  Teuchos::RefCountPtr<LinearOpBase<double> >
    A2 = Teuchos::rcp(new EpetraLinearOp(epetra_A2));
  opFactory->initializeOp(A2,&*nsA);
  
  if(out.get()) *out << "\nJ) Testing the LinearOpBase interface of nsA ...\n";

  Thyra::seed_randomize<double>(0);

  result = linearOpTester.check(*nsA,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nK) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  result = linearOpWithSolveTester.check(*nsA,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nL) Create an implicitly scaled (by 2.5) and transposed matrix A3 = scale(2.5,transpose(A)) and initialize nsA2 ...\n";

  Teuchos::RefCountPtr<const LinearOpBase<double> >
    A3 = scale(2.5,transpose(A));
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
    nsA2 = createAndInitializeLinearOpWithSolve(*opFactory,A3);
  
  if(out.get()) *out << "\nM) Testing the LinearOpBase interface of nsA2 ...\n";

  Thyra::seed_randomize<double>(0);

  result = linearOpTester.check(*nsA2,OSTab(out).getOStream().get());
  if(!result) success = false;

  if(out.get()) *out << "\nN) Testing the LinearOpWithSolveBase interface of nsA2 ...\n";
    
  result = linearOpWithSolveTester.check(*nsA2,OSTab(out).getOStream().get());
  if(!result) success = false;
  
  if(out.get()) *out << "\nO) Testing that LinearOpBase interfaces of transpose(nsA) == nsA2 ...\n";

  result = linearOpTester.compare(
    *transpose(Teuchos::rcp_implicit_cast<const LinearOpBase<double> >(nsA)),*nsA2
    ,OSTab(out).getOStream().get()
    );
  if(!result) success = false;

#else // __sun
  
  if(out.get()) *out << "\nTest failed since is was not even compiled since __sun was defined!\n";
  success = false;

#endif // __sun

  return success;

}
