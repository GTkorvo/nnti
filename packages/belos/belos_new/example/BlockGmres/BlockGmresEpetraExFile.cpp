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
// This driver reads a problem from a file, which can be in Harwell-Boeing (*.hb),
// Matrix Market (*.mtx), or triplet format (*.triU, *.triS).  The right-hand side
// from the problem, if it exists, will be used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this example. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestOutputter.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

int main(int argc, char *argv[]) {
  //
  int MyPID = 0;
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
#endif
  //
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  using Teuchos::ParameterList;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;

  bool verbose = false, proc_verbose = false;
  int frequency = -1;        // frequency of status test output.
  int blocksize = 1;         // blocksize
  int numrestarts = 15;      // number of restarts allowed 
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  int maxsubspace = 25;      // maximum number of blocks the solver can use for the subspace
  std::string filename("orsirr1.hb");
  MT tol = 1.0e-5;           // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-restarts",&numrestarts,"Number of restarts allowed for GMRES solver.");
  cmdp.setOption("blocksize",&blocksize,"Block size used by GMRES.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  //
  // Get the problem
  //
  RefCountPtr<Epetra_Map> Map;
  RefCountPtr<Epetra_CrsMatrix> A;
  RefCountPtr<Epetra_Vector> B, X;
  EpetraExt::readEpetraLinearSystem(filename, Comm, &A, &Map, &X, &B);
  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters = -1)
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  //
  ParameterList My_PL;
  My_PL.set( "Num Blocks", maxiters );               // Maximum number of blocks in Krylov factorization
  My_PL.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  My_PL.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  My_PL.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  if (verbose)
    My_PL.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary );
  else
    My_PL.set( "Verbosity", Belos::Errors + Belos::Warnings );
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<double,MV,OP>
    My_LP( A, 
	   Teuchos::rcp_implicit_cast<Epetra_MultiVector>(X), 
	   Teuchos::rcp_implicit_cast<Epetra_MultiVector>(B) 
	   );
  int numrhs = B->NumVectors();     // number of right-hand sides
  My_LP.SetBlockSize( blocksize );
  //
  // *******************************************************************
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  Belos::OutputManager<double> My_OM();
 
  // Create an iterative solver manager.
  RefCountPtr< Belos::SolverManager<double,MV,OP> > newSolver
    = rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(rcp(&My_LP,false), My_PL) );

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << blocksize << endl;
    cout << "Number of restarts allowed: " << numrestarts << endl;
    cout << "Max number of Gmres iterations per restart cycle: " << maxiters << endl; 
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
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(*Map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid ); 
  MVT::MvNorm( resid, &actual_resids );
  MVT::MvNorm( *B, &rhs_norm );
  if (proc_verbose) {
    cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
    for ( int i=0; i<numrhs; i++) {
      cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<endl;
    }
  }

  if (ret!=Belos::Converged) {
    if (proc_verbose)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  //
} 
