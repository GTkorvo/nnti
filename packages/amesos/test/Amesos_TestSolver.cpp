#include "Amesos_ConfigDefs.h"
#include "Amesos_Parameter_List.h"
#include <string>
#include "Trilinos_Util_ReadTriples2Epetra.h"
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef HAVE_AMESOS_KUNDERT
#include "KundertOO.h"
#endif
#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SLUS
#include "Epetra_SLU.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef TEST_AZTEC
#include "AztecOO.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "DscpackOO.h"
#include "Amesos_Dscpack.h"
#endif
#if 0 
#include "UmfpackOO.h"
#include "SpoolesserialOO.h"
#include "TimeMemory.h"
#include "SparseSolverResult.h"
#endif

#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h" 

#include <vector>
//
//  Amesos_TestSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers and computes the error 
//  and residual.  
//
//  It reads the matrix in on a single processor and can pass that
//  matrix to the solver or it can convert that matrix to a 
//  distributed matrix and pass the distributed matrix to the solver.
//
//  Amesos_TestSolver can test either A x = b or A^T x = b.
//  This can be a bit confusing because sparse direct solvers 
//  use compressed column storage - the transpose of Trilinos'
//  sparse row storage.
//
//  Matrices:
//    readA - Serial.  As read from the file.
//    transposeA - Serial.  The transpose of readA.
//    serialA - if (transpose) then transposeA else readA 
//    distributedA - readA distributed to all processes
//    passA - if ( distributed ) then distributedA else serialA
//
//
int Amesos_TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		       SparseSolverType SparseSolver,
		       bool transpose, 
		       int special, AMESOS_MatrixType matrix_type ) {

  int iam = Comm.MyPID() ;

  //  int whatever;
  //  if ( iam == 0 )  cin >> whatever ; 
  Comm.Barrier();


  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  string FileName = matrix_file ;
  int FN_Size = FileName.size() ; 
  string LastFiveBytes = FileName.substr( EPETRA_MAX(0,FN_Size-5), FN_Size );
  string LastFourBytes = FileName.substr( EPETRA_MAX(0,FN_Size-4), FN_Size );
  if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, false, Comm, readMap, readA, readx, 
						      readb, readxexact) );
  } else {
    if ( LastFiveBytes == ".triS" ) { 
      // Call routine to read in symmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( matrix_file, true, Comm, readMap, readA, readx, 
							readb, readxexact) );
    } else {
      if (  LastFourBytes == ".mtx" ) { 
	EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( matrix_file, Comm, readMap, 
							       readA, readx, readb, readxexact) );
      } else {
	// Call routine to read in HB problem
	Trilinos_Util_ReadHb2Epetra( matrix_file, Comm, readMap, readA, readx, 
						     readb, readxexact) ;
      }
    }
  }

  Epetra_CrsMatrix transposeA(Copy, *readMap, 0);
  Epetra_CrsMatrix *serialA ; 

  if ( transpose ) {
    assert( CrsMatrixTranspose( readA, &transposeA ) == 0 ); 
    serialA = &transposeA ; 
  } else {
    serialA = readA ; 
  }

  Epetra_RowMatrix * passA = 0; 
  Epetra_Vector * passx = 0; 
  Epetra_Vector * passb = 0;
  Epetra_Vector * passxexact = 0;
  Epetra_Vector * passresid = 0;
  Epetra_Vector * passtmp = 0;

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);

  Epetra_Export exporter(*readMap, map);
  Epetra_CrsMatrix A(Copy, map, 0);

  Epetra_Vector x(map);
  Epetra_Vector b(map);
  Epetra_Vector xexact(map);
  Epetra_Vector resid(map);
  Epetra_Vector readresid(*readMap);
  Epetra_Vector tmp(map);
  Epetra_Vector readtmp(*readMap);

  //  Epetra_Vector xcomp(map);      // X as computed by the solver
  bool distribute_matrix = ( matrix_type == AMESOS_Distributed ) ; 
  if ( distribute_matrix ) { 
    // Create Exporter to distribute read-in matrix and vectors
    //
    //  Initialize x, b and xexact to the values read in from the file
    //
    x.Export(*readx, exporter, Add);
    b.Export(*readb, exporter, Add);
    xexact.Export(*readxexact, exporter, Add);
    Comm.Barrier();
    
    A.Export(*serialA, exporter, Add);
    Comm.Barrier();

    assert(A.TransformToLocal()==0);    
    
    Comm.Barrier();

    passA = &A; 
    passx = &x; 
    passb = &b;
    passxexact = &xexact;
    passresid = &resid;
    passtmp = &tmp;
  } else { 
    passA = serialA; 
    passx = readx; 
    passb = readb;
    passxexact = readxexact;
    passresid = &readresid;
    passtmp = &readtmp;
  }


  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA, 
				 (Epetra_MultiVector *) passx, 
				 (Epetra_MultiVector *) passb );
  


  Epetra_Time TotalTime( Comm ) ; 
  //  switch( SparseSolver ) { 
  //  case UMFPACK:
  if ( false ) { 
#ifdef TEST_UMFPACK
  } else if ( SparseSolver == UMFPACK ) { 
    UmfpackOO umfpack( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    umfpack.SetTrans( transpose ) ; 
    umfpack.Solve() ; 
#endif
#ifdef HAVE_AMESOS_SLUS
  } else if ( SparseSolver == SuperLU ) { 
    Epetra_SLU superluserial( &Problem ) ;
    
    superluserial.Solve() ; 
#endif
#ifdef HAVE_AMESOS_SLUD
  } else if ( SparseSolver == SuperLUdist ) {

    for ( int i = 0; i < 1+special ; i++ ) { 
      SuperludistOO superludist( Problem ) ; 
      superludist.SetTrans( transpose );  
      EPETRA_CHK_ERR( superludist.Solve( true ) ); 
   }
#endif
#ifdef HAVE_AMESOS_SLUD2
  } else if ( SparseSolver == SuperLUdist2 ) {

    for ( int i = 0; i < 1+special ; i++ ) { 
      Superludist2_OO superludist2( Problem ) ; 
      superludist2.SetTrans( transpose );  
      EPETRA_CHK_ERR( superludist2.Solve( true ) ); 
   }
#endif
#ifdef HAVE_AMESOS_DSCPACK
  } else if ( SparseSolver == DSCPACKOLD ) {

    for ( int i = 0; i < 1+special ; i++ ) { 
      DscpackOO dscpack( Problem ) ; 
      EPETRA_CHK_ERR( dscpack.Solve( true ) ); 
   }
  } else if ( SparseSolver == DSCPACK ) {

    for ( int i = 0; i < 1+special ; i++ ) { 
      AMESOS::Parameter::List ParamList ;
      Amesos_Dscpack A_dscpack( Problem, ParamList ) ; 
      EPETRA_CHK_ERR( A_dscpack.Solve(  ) ); 
      //      cout << " After the call to A_dscpack.Solve() " << endl ; 
   }
#endif
#ifdef TEST_SPOOLES
  } else if ( SparseSolver == SPOOLES ) { 
    SpoolesOO spooles( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    spooles.SetTrans( transpose ) ; 
    spooles.Solve() ; 
#endif
#ifdef HAVE_AMESOS_KUNDERT
  } else if ( SparseSolver == KUNDERT ) { 
    KundertOO kundert( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    kundert.SetTrans( transpose ) ; 
    kundert.Solve() ; 
#endif
#ifdef TEST_SPOOLESSERIAL 
  } else if ( SparseSolver == SPOOLESSERIAL ) { 
    SpoolesserialOO spoolesserial( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    spoolesserial.Solve() ;
#endif
#ifndef TEST_AZTEC
  } else if ( SparseSolver == Aztec ) { 
    SparseDirectTimingVars::log_file 
      << "Aztec Solver won't compile for me on this platform" << endl ;
    string errormsg = "Aztec Solver won't compile for me on this platform" ;
    throw( errormsg ); 
#else
  } else if ( SparseSolver == Aztec ) { 
    AztecOO aztec( (Epetra_RowMatrix *) passA, 
		       (Epetra_MultiVector *) passx, 
		       (Epetra_MultiVector *) passb ) ; 
    
    aztec.Iterate(1000, 1.0e-10 * Anorm ) ; 
#endif
  } else { 
    SparseDirectTimingVars::log_file << "Solver not implemented yet" << endl ;
    cerr << "\n\n####################  Requested solver not available on this platform #####################\n" << endl ;
  }

  SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 
  SparseDirectTimingVars::SS_Result.Set_First_Time( 0.0 ); 
  SparseDirectTimingVars::SS_Result.Set_Middle_Time( 0.0 ); 
  SparseDirectTimingVars::SS_Result.Set_Last_Time( 0.0 ); 

  //
  //  Compute the error = norm(xcomp - xexact )
  //
  double error;
  passresid->Update(1.0, *passx, -1.0, *passxexact, 0.0);

  passresid->Norm2(&error);
  SparseDirectTimingVars::SS_Result.Set_Error(error) ;

  //  passxexact->Norm2(&error ) ; 
  //  passx->Norm2(&error ) ; 

  //
  //  Compute the residual = norm(Ax - b)
  //
  double residual ; 
  passA->Multiply( transpose, *passx, *passtmp);
  //  passresid->Norm2(&error);    BOGUS - I claim
  passresid->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
  passresid->Norm2(&residual);

  SparseDirectTimingVars::SS_Result.Set_Residual(residual) ;
    
  double bnorm; 
  passb->Norm2( &bnorm ) ; 
  SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm) ;

  double xnorm; 
  passx->Norm2( &xnorm ) ; 
  SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm) ;

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  
  Comm.Barrier();

  return 0;
}
