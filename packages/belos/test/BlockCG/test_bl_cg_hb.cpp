// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "createEpetraProblem.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

int main(int argc, char *argv[]) {
  //
#ifdef EPETRA_MPI	
  // Initialize MPI	
  MPI_Init(&argc,&argv); 	
  Belos::MPIFinalize mpiFinalize; // Will call finalize with *any* return
#endif
  //
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, proc_verbose = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;  // total number of right-hand sides to solve for
  int blocksize = 1;  // blocksize used by solver
  int maxiters = -1;  // maximum number of iterations for solver to use
  std::string filename("bcsstk14.hb");
  double tol = 1.0e-5;  // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size to be used by CG solver.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  //
  // Get the problem
  //
  int MyPID;
  RefCountPtr<Epetra_CrsMatrix> A;
  RefCountPtr<Epetra_MultiVector> B, X;
  int return_val =Belos::createEpetraProblem(filename,NULL,&A,&B,&X,&MyPID);
  if(return_val != 0) return return_val;
  proc_verbose = ( verbose && (MyPID==0) );
  //
  // Solve using Belos
  //
  typedef double                          ST;
  typedef Epetra_Operator                 OP;
  typedef Epetra_MultiVector              MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  //
  // *****Construct initial guess and random right-hand-sides *****
  //
  if (numrhs != 1) {
    X = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
    MVT::MvRandom( *X );
    B = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );
  }
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  //
  ParameterList belosList;
  belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  if (verbose) {
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + 
		   Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  else
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<double,MV,OP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
    return -1;
  }
  //
  // Create an iterative solver manager.
  //
  RefCountPtr< Belos::SolverManager<double,MV,OP> > newSolver
    = rcp( new Belos::BlockCGSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );
  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << blocksize << endl;
    cout << "Max number of CG iterations: " << maxiters << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = newSolver->solve();
  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(A->Map(), numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid ); 
  MVT::MvNorm( resid, &actual_resids );
  MVT::MvNorm( *B, &rhs_norm );
  if (proc_verbose) {
    cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      cout<<"Problem "<<i<<" : \t"<< actRes <<endl;
      if (actRes > tol) badRes = true;
    }
  }

  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose)
      cout << endl << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose)
    cout << endl << "End Result: TEST PASSED" << endl;
  return 0;
  //
} // end test_bl_cg_hb.cpp




