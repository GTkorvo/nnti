#include "TSFCoreSolversGmresSolver.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
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

	bool success = true;
	bool verbose = false;

	try {

		int i, j;
		int n_nonzeros, N_update;
		int *bindx=0, *update=0, *col_inds=0;
		double *val=0, *row_vals=0;
		double *xguess=0, *b=0, *xexact=0, *xsolve=0;

#ifdef EPETRA_MPI	
		// Initialize MPI	
		MPI_Init(&argc,&argv); 	
		Epetra_MpiComm Comm( MPI_COMM_WORLD );	
#else	
		Epetra_SerialComm Comm;	
#endif
	
		int MyPID = Comm.MyPID();
		int NumProc = Comm.NumProc();
		// Set verbosity of output
		verbose = false;
		for( i=1; i<argc; i++ ) {
			if (argv[i][0]=='-' && argv[i][1]=='v' && MyPID==0 ) { verbose = true; };
		}
		if ( verbose ) { argc--; } // Decrement argument counter if one of the arguments is the verbosity flag.
		//
		if( argc < 4 ) {
			if ( verbose ) {
				cerr
					<< "Usage: " << argv[0] 
					<< " HB_filename max_iter tol [-v] " << endl
					<< "where:" << endl
					<< "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
					<< "max_iter           - maximum number of iterations allowed in GMRES solve" <<endl
					<< "tol                - relative residual tolerance for GMRES solve" << endl
					<< "[-v]               - verbosity flag for debugging" << endl
					<< endl;
			}
			return(1);
		}

		if ( verbose ) {
			std::cout
				<< std::endl << std::endl
				<< "***\n"
				<< "*** Testing simple GMRES solver on test problem \'" << argv[1] << "\'\n"
				<< "***\n\n";
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
		Trilinos_Util_read_hb(argv[1], MyPID, &NumGlobalElements, &n_nonzeros,
							  &val, &bindx, &xguess, &b, &xexact);
		//
		// *****Distribute data among processors*****
		//
		Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalElements, &n_nonzeros, &N_update,
										 &update, &val, &bindx, &xguess, &b, &xexact);
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
		Teuchos::RefCountPtr<Epetra_CrsMatrix> A =
			Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, NumNz) );
		//
		// Create some Epetra_Vectors
		//
		Teuchos::RefCountPtr<Epetra_Vector> xx = Teuchos::rcp( new Epetra_Vector(Copy, Map, xexact) );
		Teuchos::RefCountPtr<Epetra_Vector> bb = Teuchos::rcp( new Epetra_Vector(Copy, Map, b) );
		// Solution vector, initialize to zero.
		Teuchos::RefCountPtr<Epetra_Vector> x = Teuchos::rcp( new Epetra_Vector(Map) ); 
		x->PutScalar( 0.0 ); 
		// 
		// Add rows one-at-a-time
		//
		int NumEntries;
		for (i=0; i<NumMyElements; i++) {
			row_vals = val + bindx[i];
			col_inds = bindx + bindx[i];
			NumEntries = bindx[i+1] - bindx[i];
			assert(A->InsertGlobalValues(update[i], NumEntries, row_vals, col_inds)==0);
			assert(A->InsertGlobalValues(update[i], 1, val+i, update+i)==0);
		}
		//
		// Finish up
		//
		assert(A->TransformToLocal()==0);
		//
		A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
		//
		// Create the TSFCore Linear Operator and Vectors
		//
		TSFCore::EpetraLinearOp ELOp( A ); 
		TSFCore::EpetraVector RHS( bb );
		TSFCore::EpetraVector Soln( x );
		//
		// ********Other information used by the GMRES solver***********
		//
		int max_iter = atoi(argv[2]);  // maximum number of iterations to run
		double tol = atof(argv[3]);  // relative residual tolerance
		//
		//*******************************************************************
		// *************Start the GMRES iteration*************************
		//*******************************************************************
		//
		TSFCore::Solvers::GMRESSolver<double> MySolver( max_iter, tol );
		MySolver.solve( ELOp, RHS, &Soln, TSFCore::NOTRANS, max_iter, tol );
		//
		// Compute actual residual norm.
		//
		double bnorm = norm_2( RHS );
		ELOp.apply( TSFCore::NOTRANS, Soln, &RHS, 1.0, -1.0 );
		//
		double final_rel_err = norm_2( RHS )/bnorm;
		if (verbose) {
			cout << "******************* Results ************************"<<endl;
			cout << "Iteration "<< MySolver.currIteration()<<" of "<<max_iter<< endl;
			cout << "Final Computed GMRES Residual Norm : " << 
				MySolver.currEstRelResidualNorm() << endl;
			cout << "Actual Computed GMRES Residual Norm : " <<
				final_rel_err << endl;
			cout << "****************************************************"<<endl;
		}
      	
		// Release all objects  	
		delete [] NumNz;
		delete [] bindx;

		success = (final_rel_err <= tol);

	} // end try
	catch( const std::exception &excpt ) {
		std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
		success = false;
	}
	catch( ... ) {
		std::cerr << "*** Caught an unknow exception\n";
		success = false;
	}

	if (verbose) {
		if(success)
			std::cout << "\nCongratulations! the system was solved to the specified tolerance!\n";
		else
			std::cout << "\nOh no! the system was not solved to the specified tolerance!\n";
	}
	
	return success ? 0 : -1;

} // end TSFCoreSolversGmresTest.cpp
