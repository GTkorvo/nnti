#include "Amesos_ConfigDefs.h"
#include "Amesos_Parameter_List.h"
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
#ifdef HAVE_AMESOS_SLUD
#include "SuperludistOO.h"
#endif
#ifdef HAVE_AMESOS_SLUD2
#include "Superludist2_OO.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "DscpackOO.h"
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef TEST_SPOOLES
#include "SpoolesOO.h"
#endif
#ifdef TEST_AZTEC
#include "AztecOO.h"
#endif

#include "Amesos_TestSolver.h"
#include "CrsMatrixTranspose.h"
#include "SparseDirectTimingVars.h"

#include <vector>
//
//  Amesos_TestMultiSolver.cpp reads in a matrix in Harwell-Boeing format, 
//  calls one of the sparse direct solvers, using blocked right hand sides
//  and computes the error and residual.  
//
//  TestSolver ignores the Harwell-Boeing right hand sides, creating
//  random right hand sides instead.  
//
//  Amesos_TestMultiSolver can test either A x = b or A^T x = b.
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
int Amesos_TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves, 
		      SparseSolverType SparseSolver, bool transpose,
		      int special, AMESOS_MatrixType matrix_type ) {


  int iam = Comm.MyPID() ;

  
  //  int hatever;
  //  if ( iam == 0 )  cin >> hatever ; 
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

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);

  // Create Exporter to distribute read-in matrix and vectors
  Epetra_Export exporter(*readMap, map);
  Epetra_CrsMatrix A(Copy, map, 0);

  Epetra_RowMatrix * passA = 0; 
  Epetra_MultiVector * passx = 0; 
  Epetra_MultiVector * passb = 0;
  Epetra_MultiVector * passxexact = 0;
  Epetra_MultiVector * passresid = 0;
  Epetra_MultiVector * passtmp = 0;

  Epetra_MultiVector x(map,numsolves);
  Epetra_MultiVector b(map,numsolves);
  Epetra_MultiVector xexact(map,numsolves);
  Epetra_MultiVector resid(map,numsolves);
  Epetra_MultiVector tmp(map,numsolves);

  Epetra_MultiVector serialx(*readMap,numsolves);
  Epetra_MultiVector serialb(*readMap,numsolves);
  Epetra_MultiVector serialxexact(*readMap,numsolves);
  Epetra_MultiVector serialresid(*readMap,numsolves);
  Epetra_MultiVector serialtmp(*readMap,numsolves);

  bool distribute_matrix = ( matrix_type == AMESOS_Distributed ) ; 
  if ( distribute_matrix ) { 
    //
    //  Initialize x, b and xexact to the values read in from the file
    //
    
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
    passx = &serialx; 
    passb = &serialb;
    passxexact = &serialxexact;
    passresid = &serialresid;
    passtmp = &serialtmp;
  }

  passxexact->SetSeed(1.31) ; 
  passxexact->Random();
  passx->SetSeed(11.231) ; 
  passx->Random();

  passb->PutScalar( 0.0 );
  passA->Multiply( transpose, *passxexact, *passb ) ; 

  double Anorm = passA->NormInf() ; 
  SparseDirectTimingVars::SS_Result.Set_Anorm(Anorm) ;

  Epetra_LinearProblem Problem(  (Epetra_RowMatrix *) passA, 
				 (Epetra_MultiVector *) passx, 
				 (Epetra_MultiVector *) passb );

  double max_resid = 0.0;
  for ( int i = 0 ; i < special+1 ; i++ ) { 
    
    Epetra_Time TotalTime( Comm ) ; 
    if ( false ) { 
#ifdef TEST_UMFPACK

      unused code

    } else if ( SparseSolver == UMFPACK ) { 
      UmfpackOO umfpack( (Epetra_RowMatrix *) passA, 
			 (Epetra_MultiVector *) passx, 
			 (Epetra_MultiVector *) passb ) ; 
    
      umfpack.SetTrans( transpose ) ; 
      umfpack.Solve() ; 
#endif
#ifdef TEST_SUPERLU
    } else if ( SparseSolver == SuperLU ) { 
      SuperluserialOO superluserial( (Epetra_RowMatrix *) passA, 
				     (Epetra_MultiVector *) passx, 
				     (Epetra_MultiVector *) passb ) ; 

      superluserial.SetPermc( SuperLU_permc ) ; 
      superluserial.SetTrans( transpose ) ; 
      superluserial.SetUseDGSSV( special == 0 ) ; 
      superluserial.Solve() ; 
#endif
#ifdef HAVE_AMESOS_SLUD
    } else if ( SparseSolver == SuperLUdist ) { 
      SuperludistOO superludist( Problem ) ; 
      superludist.SetTrans( transpose ) ; 
      EPETRA_CHK_ERR( superludist.Solve( true ) ) ;
#endif 
#ifdef HAVE_AMESOS_SLUD2
    } else if ( SparseSolver == SuperLUdist2 ) { 
      Superludist2_OO superludist2( Problem ) ; 
      superludist2.SetTrans( transpose ) ; 
      EPETRA_CHK_ERR( superludist2.Solve( true ) ) ;
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
#ifdef HAVE_AMESOS_DSCPACK
    } else if ( SparseSolver == DSCPACK ) { 
      AMESOS::Parameter::List ParamList ;
      Amesos_Dscpack dscpack( Problem, ParamList ) ; 
    
      EPETRA_CHK_ERR( dscpack.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_UMFPACK
    } else if ( SparseSolver == UMFPACK ) { 
      AMESOS::Parameter::List ParamList ;
      Amesos_Umfpack umfpack( Problem, ParamList ) ; 
      EPETRA_CHK_ERR( umfpack.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( umfpack.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_KLU
    } else if ( SparseSolver == KLU ) { 
      AMESOS::Parameter::List ParamList ;
      Amesos_Klu klu( Problem, ParamList ) ; 
      EPETRA_CHK_ERR( klu.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( klu.SymbolicFactorization(  ) ); 
      EPETRA_CHK_ERR( klu.NumericFactorization(  ) ); 
      EPETRA_CHK_ERR( klu.Solve( ) ); 
#endif
#ifdef HAVE_AMESOS_MUMPS
    } else if ( SparseSolver == MUMPS ) { 
      AMESOS::Parameter::List ParamList ;
      Amesos_Mumps mumps( Problem, ParamList ) ; 
      EPETRA_CHK_ERR( mumps.SetUseTranspose( transpose ) ); 
    
      EPETRA_CHK_ERR( mumps.Solve( ) ); 
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
      cerr << "\n\n####################  Requested solver not available (Or not tested with blocked RHS) on this platform #####################\n" << endl ;
    }

    SparseDirectTimingVars::SS_Result.Set_Total_Time( TotalTime.ElapsedTime() ); 
    SparseDirectTimingVars::SS_Result.Set_First_Time( 0.0 ); 
    SparseDirectTimingVars::SS_Result.Set_Middle_Time( 0.0 ); 
    SparseDirectTimingVars::SS_Result.Set_Last_Time( 0.0 ); 

    //
    //  Compute the error = norm(xcomp - xexact )
    //
    vector <double> error(numsolves) ; 
    double max_error = 0.0;
  
    passresid->Update(1.0, *passx, -1.0, *passxexact, 0.0);

    passresid->Norm2(&error[0]);
    for ( int i = 0 ; i< numsolves; i++ ) 
      if ( error[i] > max_error ) max_error = error[i] ; 
    SparseDirectTimingVars::SS_Result.Set_Error(max_error) ;

    //  passxexact->Norm2(&error[0] ) ; 
    //  passx->Norm2(&error ) ; 

    //
    //  Compute the residual = norm(Ax - b)
    //
    vector <double> residual(numsolves) ; 
  
    passtmp->PutScalar(0.0);
    passA->Multiply( transpose, *passx, *passtmp);
    passresid->Update(1.0, *passtmp, -1.0, *passb, 0.0); 
    passresid->Norm2(&residual[0]);

    for ( int i = 0 ; i< numsolves; i++ ) 
      if ( residual[i] > max_resid ) max_resid = residual[i] ; 


    SparseDirectTimingVars::SS_Result.Set_Residual(max_resid) ;
    
    vector <double> bnorm(numsolves); 
    passb->Norm2( &bnorm[0] ) ; 
    SparseDirectTimingVars::SS_Result.Set_Bnorm(bnorm[0]) ;

    vector <double> xnorm(numsolves); 
    passx->Norm2( &xnorm[0] ) ; 
    SparseDirectTimingVars::SS_Result.Set_Xnorm(xnorm[0]) ;

#if 0
    if ( iam == 0 )     cout << " Here is B: " << endl; 
    passb->Print( cout ) ;
    Comm.Barrier(); 
    if ( iam == 0 )    cout << endl <<" That was B " << endl ; 
#endif

    if ( false && iam == 0 ) { 

      cout << " Amesos_TestMutliSolver.cpp " << endl ; 
      for ( int i = 0 ; i< numsolves && i < 10 ; i++ ) {
	cout << "i=" << i 
	     << " error = " << error[i] 
	     << " xnorm = " << xnorm[i] 
	     << " residual = " << residual[i] 
	     << " bnorm = " << bnorm[i] 
	     << endl ; 
      
      }
    
      cout << endl << " max_resid = " << max_resid ; 
      cout << " max_error = " << max_error << endl ; 
      cout << " Get_residual() again = " << SparseDirectTimingVars::SS_Result.Get_Residual() << endl ;

    }
  }
  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  
  Comm.Barrier();

return 0 ;
}
