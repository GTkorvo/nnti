#include "copyright.h"
//
//  CallAmesosFactory.ccp shows how to call Amesos_Factory() with 
//  different right hand sides for each solve.  It performs three solves:
//  Solve Ax = b
//  Solve Ax1 = x
//  Solve Ax2 = x1
//
#include "Amesos_config.h"
#include "Amesos_Factory.h"
#include "Amesos_Parameter_List.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "CreateTridi.h"

int main(int argc, char *argv[])
{


#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif



  int iam = Comm.MyPID() ; 

  const int NumPoints = 10;  // Must be between 2 and 100 (on large matrices,
                             // the problem is quite ill-conditioned) 

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Epetra_Map Map(NumPoints, 0, Comm);

  //  Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix A(Copy, Map, 0);

  //
  //  Populate A with a [-1,2,-1] tridiagonal matrix
  //  See CreateTridi.cpp in this directory 
  CreateTridi( A ) ; 
  
  Epetra_Vector x2(Map), x1(Map), x(Map), b(Map), residual(Map), temp(Map);

  //
  //  Solve Ax = b using Amesos_KLU via the Amesos_Factory interface
  //
  AMESOS::Parameter::List ParamList ;
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos_Factory Afactory;
  //
  //  Note that Abase is created with an empty Problem, none of A, x or b
  //  have been specified at this point.  
  Abase = Afactory.Create( AMESOS_MUMPS, Problem, ParamList ) ; 
  if ( Abase == 0 ) {
    cout << " AMESOS_KLU not implemented " << endl ; 
    exit(13);
  }

  //
  //  Factor A
  //
  Problem.SetOperator( &A );
  EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
  EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

  b.Random();
  //
  //  Solve Ax = b 
  //
  Problem.SetLHS( &x );
  Problem.SetRHS( &b );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //
  //  Solve Ax1 = x 
  //
  Problem.SetLHS( &x1 );
  Problem.SetRHS( &x );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //
  //  Solve Ax2 = x1
  //
  Problem.SetLHS( &x2 );
  Problem.SetRHS( &x1 );
  EPETRA_CHK_ERR( Abase->Solve(  ) ); 

  //
  //  Compute the residual: A^3 x2 - b
  //

  A.Multiply( false, x2, temp ) ; //  temp = A x2
  A.Multiply( false, temp, x2 ) ; //  x2 = A^2 x2
  A.Multiply( false, x2, temp ) ; //  temp = A^3 x2
  residual.Update( 1.0, temp, -1.0, b, 0.0 ) ;

  double norm_residual ;
  residual.Norm2( &norm_residual ) ; 

  if (iam == 0 ) {
    cout << " norm2(A^3 x-b) = " << norm_residual << endl ; 
    //
    //  This is an ill-conditioned problem
    //
    if ( norm_residual < (1e-15)*(1.0*NumPoints*NumPoints*NumPoints*
				  NumPoints*NumPoints*NumPoints) )
      cout << " Test Passed " << endl ;
    else
      cout << " TEST FAILED " << endl ;
  }

  delete Abase;

#ifdef EPETRA_MPI
  MPI_Finalize() ; 
#endif
  return 0;
}

