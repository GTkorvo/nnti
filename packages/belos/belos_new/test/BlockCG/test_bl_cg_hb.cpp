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
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCG.hpp"
#include "createEpetraProblem.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_Time.hpp"

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
  Teuchos::Time timer("Belos CG");
  //
  // Get the problem
  //
  RefCountPtr<Epetra_CrsMatrix> A;
  RefCountPtr<Epetra_MultiVector> B, X;
  int MyPID;
  bool verbose;
  int return_val =Belos::createEpetraProblem(argc,argv,NULL,&A,&B,&X,&MyPID,&verbose);
  if(return_val != 0) return return_val;
  //
  // Solve using Belos
  //
  typedef Belos::Operator<double> OP;
  typedef Belos::MultiVec<double> MV;
  //
  // Construct a Belos::Operator instance through the Petra interface.
  //
  Belos::EpetraOp Amat( A );
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  int numrhs = 1;  // total number of right-hand sides to solve for
  int block = 1;  // blocksize used by solver
  int maxits = NumGlobalElements/block; // maximum number of iterations to run
  double tol = 1.0e-5;  // relative residual tolerance
  //
  // Construct the right-hand side and solution multivectors.
  //
  Belos::EpetraMultiVec rhs(*B);
  Belos::EpetraMultiVec soln(*X);
  //
  //  Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<double,MV,OP> My_LP( rcp(&Amat,false), rcp(&soln,false), rcp(&rhs,false) );
  My_LP.SetBlockSize( block );
  //
  // *******************************************************************
  // *************Start the block CG iteration*************************
  // *******************************************************************
  //
  Belos::StatusTestMaxIters<double,MV,OP> test1( maxits );
  Belos::StatusTestResNorm<double,MV,OP> test2( tol );
  Belos::StatusTestCombo<double,MV,OP> My_Test( Belos::StatusTestCombo<double,MV,OP>::OR, test1, test2 );
  //
  Belos::OutputManager<double> My_OM( MyPID );
  if (verbose)
    My_OM.SetVerbosity( 2 );	
  //
  Belos::BlockCG<double,MV,OP>
    MyBlockCG(rcp(&My_LP,false), rcp(&My_Test,false), rcp(&My_OM,false) );
  //
  // **********Print out information about problem*******************
  //
  if (verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << block << endl;
    cout << "Max number of CG iterations: " << maxits << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }
  //
  //
  if (verbose) {
    cout << endl << endl;
    cout << "Running Block CG -- please wait" << endl;
    cout << (numrhs+block-1)/block 
	 << " pass(es) through the solver required to solve for " << endl; 
    cout << numrhs << " right-hand side(s) -- using a block size of " << block
	 << endl << endl;
  }
  timer.start(true);
  MyBlockCG.Solve();
  timer.stop();
  
  if (verbose) {
    cout << "Solution time : "<< timer.totalElapsedTime()<<endl;
  }
  
  // Release all objects  

  if (My_Test.GetStatus() == Belos::Converged) {
    cout<< "***************** The test PASSED !!!********************"<<endl;
    return 0;
  }
  cout<< "********************The test FAILED!!! ********************"<<endl;
  return 1; 
  //
} // end test_bl_cg_hb.cpp




