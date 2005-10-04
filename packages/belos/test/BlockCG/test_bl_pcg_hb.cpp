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
// Multiple right-hand-sides are created randomly.
// The initial guesses are all set to zero. 
//
// As currently set up, this driver tests the case when the number of right-hand
// sides (numrhs = 15) is greater than the blocksize (block = 10) used by 
// the solver. Here, 2 passes through the solver are required to solve 
// for all right-hand sides. This information can be edited (see below - other
// information used by block solver - can be user specified) to solve for
// other sizes of systems. For example, one could set numrhs = 1 and block = 1,
// to solve a single right-hand side system in the traditional way, or, set
// numrhs = 1 and block > 1 to solve a single rhs-system with a block implementation. 
//
// 
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosPetraInterface.hpp"
#include "BelosBlockCG.hpp"
#include "Trilinos_Util.h"
#include "Ifpack_CrsIct.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Time.hpp"
//
//
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
	//
	int i;
	int n_nonzeros, N_update;
	int *bindx=0, *update=0, *col_inds=0;
	double *val=0, *row_vals=0;
	Teuchos::Time timer("Belos Preconditioned CG");
	
#ifdef EPETRA_MPI	
	// Initialize MPI	
	MPI_Init(&argc,&argv); 	
	Epetra_MpiComm Comm( MPI_COMM_WORLD );	
#else	
	Epetra_SerialComm Comm;	
#endif
	
	int MyPID = Comm.MyPID();
	bool verbose = 0;
	//
        if((argc < 2)&& MyPID==0) {
        cerr << "Usage: " << argv[0]
         << " [ -v ] [ HB_filename ]" << endl
         << "where:" << endl
         << "-v                 - run test in verbose mode" << endl
         << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
         << endl;
        return(1);
        }
        //
        // Find verbosity flag
        //
        int file_arg = 1;
        for(i = 1; i < argc; i++)
        {
          if(argv[i][0] == '-' && argv[i][1] == 'v') {
            verbose = (MyPID == 0);
            if(i==1) file_arg = 2;
          }
	}
	//
	//**********************************************************************
	//******************Set up the problem to be solved*********************
	//**********************************************************************
    	//
    	int NumGlobalElements;  // total # of rows in matrix
	//
	// *****Read in matrix from HB file******
	//
	Trilinos_Util_read_hb(argv[file_arg], MyPID, &NumGlobalElements, &n_nonzeros, &val, 
		                    &bindx);
	//
	// *****Distribute data among processors*****
	//
	Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalElements, &n_nonzeros, &N_update,
		                             &update, &val, &bindx);
	//
	//
    	// ********Other information used by block solver***********
	//*****************(can be user specified)******************
	//
	int numrhs = 15;  // total number of right-hand sides to solve for
    	int block = 10;  // blocksize used by solver
    	int maxits = NumGlobalElements/block - 1; // maximum number of iterations to run
    	double tol = 1.0e-6;  // relative residual tolerance
	//
	//*************************************************************
	//
	// *****Construct the matrix*****
	//
	int NumMyElements = N_update; // # local rows of matrix on processor
	//
    	// Create an integer vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
	// on this processor
	//
	int * NumNz = new int[NumMyElements];
	for (i=0; i<NumMyElements; i++) {
		NumNz[i] = bindx[i+1] - bindx[i] + 1;
	}
	//
	Epetra_Map Map(NumGlobalElements, NumMyElements, update, 0, Comm);
	//
	// Create a Epetra_Matrix
	//
	Epetra_CrsMatrix A(Copy, Map, NumNz);
	//
	// Add rows one-at-a-time
	//
	int NumEntries;
	for (i=0; i<NumMyElements; i++) {
		row_vals = val + bindx[i];
		col_inds = bindx + bindx[i];
		NumEntries = bindx[i+1] - bindx[i];
		assert(A.InsertGlobalValues(update[i], NumEntries, row_vals, col_inds)==0);
		assert(A.InsertGlobalValues(update[i], 1, val+i, update+i)==0);
	}
	//
	// Finish up
	//
	assert(A.FillComplete()==0);
	//
	// call the ctor that calls the petra ctor for a matrix
	//
	Belos::PetraMat<double> Amat( &A );
	//
	A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
	//
	//*****Select the Preconditioner*****
	//
	if (verbose) cout << endl << endl;
	if (verbose) cout << "Constructing ICT preconditioner" << endl;
	int Lfill = 0;
	// if (argc > 2) Lfill = atoi(argv[2]);
	if (verbose) cout << "Using Lfill = " << Lfill << endl;
	int Overlap = 0;
	// if (argc > 3) Overlap = atoi(argv[3]);
	if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
	double Athresh = 0.0;
	// if (argc > 4) Athresh = atof(argv[4]);
	if (verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
	double Rthresh = 1.0;
	// if (argc >5) Rthresh = atof(argv[5]);
	if (verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
	double dropTol = 1.0e-6;
	//
	Ifpack_CrsIct* ICT = 0;
	//
	if (Lfill > -1) {
		ICT = new Ifpack_CrsIct(A, dropTol, Lfill);
		ICT->SetAbsoluteThreshold(Athresh);
		ICT->SetRelativeThreshold(Rthresh);
		int initerr = ICT->InitValues(A);
		if (initerr != 0) cout << "InitValues error = " << initerr;
		assert(ICT->Factor() == 0);
	}
	//
	bool transA = false;
	double Cond_Est;
	ICT->Condest(transA, Cond_Est);
	if (verbose) {
		cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
		cout << endl;
	}
	Epetra_Operator& prec = dynamic_cast<Epetra_Operator&>(*ICT);
	//
	// call the ctor for the preconditioning object
	//
	Belos::PetraPrec<double> EpetraOpPrec( &prec );
	//
	//*****Construct initial guess and random right-hand-sides *****
	//
	Belos::PetraVec<double> soln(Map, numrhs);
	Belos::PetraVec<double> rhs(Map, numrhs);
	rhs.MvRandom();
	//
	//*****Create Linear Problem for Belos Solver
	//
	Belos::LinearProblemManager<double> My_LP(&Amat, &soln, &rhs);
	My_LP.SetLeftPrec( &EpetraOpPrec );
	My_LP.SetBlockSize( block );
	//
	//
	//*******************************************************************
	// *************Start the block CG iteration*************************
	//*******************************************************************
	//
        Belos::StatusTestMaxIters<double> test1( maxits );
        Belos::StatusTestResNorm<double> test2( tol );
        Belos::StatusTestCombo<double> My_Test( Belos::StatusTestCombo<double>::OR, test1, test2 );

	Belos::OutputManager<double> My_OM( MyPID );
	if (verbose)
	  My_OM.SetVerbosity( 2 );

	Belos::BlockCG<double> MyBlockCG(My_LP, My_Test, My_OM);
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

	if (verbose) {
	   cout << endl << endl;
	   cout << "Running Block CG -- please wait" << endl;
	   cout << (numrhs+block-1)/block 
		    << " pass(es) through the solver required to solve for " << endl; 
	   cout << numrhs << " right-hand side(s) -- using a block size of " << block
			<< endl << endl;
	}
	timer.start();
	MyBlockCG.Solve();	
	timer.stop();
        //
        // Compute actual residuals.
        //
        double* actual_resids = new double[numrhs];
        double* rhs_norm = new double[numrhs];
        Belos::PetraVec<double> resid( Map, numrhs );
        Amat.Apply( soln, resid );
        resid.MvAddMv( -1.0, resid, 1.0, rhs );
        resid.MvNorm( actual_resids );
        rhs.MvNorm( rhs_norm );
        if (verbose) {
          cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
          for (i=0; i<numrhs; i++) {
                cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<endl;
          }
	  cout<<endl;
        }
	if (verbose)
	  cout << "Solution time : "<< timer.totalElapsedTime()<<endl;
	
  // Release all objects  
  if (ICT) { delete ICT; ICT = 0; }
  delete [] NumNz;
  delete [] bindx;
  delete [] update;
  delete [] val; 
  delete [] actual_resids;
  delete [] rhs_norm;
	
#ifdef EPETRA_MPI

  MPI_Finalize() ;

#endif

  if (My_Test.GetStatus()!=Belos::Converged) {
	if (verbose)
      		cout << "End Result: TEST FAILED" << endl;	
	return -1;
  }
  //
  // Default return value
  //
  if (verbose)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
  //
} // end test_bl_pcg_hb.cpp

