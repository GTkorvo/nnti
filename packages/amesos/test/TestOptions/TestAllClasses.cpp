#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "TestAllClasses.h"
#include "TestOtherClasses.h"
#include "TestSuperludist.h"
#include "TestScalapack.h"
#include "TestKlu.h"
 
int TestAllClasses( const vector<string> AmesosClasses,
		    int EpetraMatrixType,
		    const vector<bool> AmesosClassesInstalled,
		    Epetra_CrsMatrix *& Amat, 
		    const bool transpose, 
		    const bool verbose, 
		    const bool symmetric, 
		    const int Levels,
		    const double Rcond,
		    int Diagonal,
		    int ReindexRowMap,
		    int ReindexColMap,
		    int RangeMapType,
		    int DomainMapType,
		    bool distribute,
		    char *filename,
		    double &maxrelerror, 
		    double &maxrelresidual,
		    int &NumTests ) {

  assert( NumTests == 0 ) ; 

  bool RowMapEqualsColMap = ( ReindexColMap == 0 ) ; 

  string StringFilename = filename ; 
  bool bcsstk04 = ( StringFilename.find("bcsstk04") < StringFilename.find("xdz_notaname_garbage") );
  bool Khead = ( StringFilename.find("Khead") < StringFilename.find("xdz_notaname_garbage") );
  bool Superlu_rua = ( StringFilename.find("Superlu_rua") < StringFilename.find("xdz_notaname_garbage") );
  bool ImpcolB = ( StringFilename.find("ImpcolB") < StringFilename.find("xdz_notaname_garbage") );
  bool a662_bus_out = ( StringFilename.find("662_bus_out") < StringFilename.find("xdz_notaname_garbage") );
  bool MissingADiagonal = ( StringFilename.find("MissingADiagonal") < StringFilename.find("xdz_notaname_garbage") ) ||
    ( StringFilename.find("ImpcolB.rua") < StringFilename.find("xdz_notaname_garbage") );

  const Epetra_Map& row_map = Amat->RowMap() ; 

  const int NumAmesosClasses = AmesosClasses.size();
  int errors = 0 ;

  if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && EpetraMatrixType == 1 ) 
    return 0 ;   //  Can't reindex a RowMatrix because we don't know the indices up front 

  for (int i=0; i < NumAmesosClasses; i++ ) {
    if ( AmesosClassesInstalled[i] ) { 
      int Errors = 0 ; 
      int NumTheseTests = 0 ; 
      if ( Amat->Comm().MyPID() == 0 ) {
	bool ReIndex = ReindexRowMap || ReindexColMap ; 
	if ( ( verbose  &&  ( ! ReIndex ) ) { 
	
	  cout << __FILE__ << "::"  << __LINE__
	       << " Perhaps about to test " 
	       << AmesosClasses[i] << " "  
	       << " EpetraMatrixType = " <<  EpetraMatrixType 
	       << " transpose = " <<  transpose 
	       << " symmetric = " <<  symmetric 
	       << " Levels = " <<  Levels 
	       << " Diagonal = " <<  Diagonal 
	       << " ReindexRowMap = " <<  ReindexRowMap 
	       << " ReindexColMap = " <<  ReindexColMap 
	       << " DomainMapType = " <<  DomainMapType 
	       << " RangeMapType = " <<  RangeMapType 
	       << " distribute = " <<  distribute 
	       << " filename = " <<  filename 
	       << endl ;  
	}
      }
      if ( AmesosClasses[i] == "Amesos_Scalapack") { 
	bool RunScalapackTest = true;
	if ( ReindexRowMap || ReindexColMap ) RunScalapackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunScalapackTest = false ;   //  Bug #1403
	if ( RunScalapackTest && verbose) cout << " Testing SCALAPACK " << endl ; 
	if ( RunScalapackTest )	Errors = TestScalapack( Amat, 
							 EpetraMatrixType,
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 maxrelerror, 
							 maxrelresidual, 
							 NumTheseTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Umfpack" ) {
	bool RunUmfpackTest = true;
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) && row_map.DistributedGlobal() ) 
	  RunUmfpackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunUmfpackTest = false ;   //  Bug #1403

	if ( RunUmfpackTest && verbose) cout << " Testing UMFPACK " << endl ; 
	
	if ( RunUmfpackTest ) Errors = TestOtherClasses("Amesos_Umfpack",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 false,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTheseTests ) ;
#if 0

	This ought to work, but it crashes

      } else if ( AmesosClasses[i] == "Amesos_Klu" ) {
	bool RunKluTest = true;
	//  We only test reindexing on klu and paraklete
	if ( ( verbose  &&  ( ReIndex ) ) { 
	  
	  cout << __FILE__ << "::"  << __LINE__
	       << " Perhaps about to test " 
	       << AmesosClasses[i] << " "  
	       << " EpetraMatrixType = " <<  EpetraMatrixType 
	       << " transpose = " <<  transpose 
	       << " symmetric = " <<  symmetric 
	       << " Levels = " <<  Levels 
	       << " Diagonal = " <<  Diagonal 
	       << " ReindexRowMap = " <<  ReindexRowMap 
	       << " ReindexColMap = " <<  ReindexColMap 
	       << " DomainMapType = " <<  DomainMapType 
	       << " RangeMapType = " <<  RangeMapType 
	       << " distribute = " <<  distribute 
	       << " filename = " <<  filename 
	       << endl ;  
	}
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) && row_map.DistributedGlobal() ) 
	  RunKluTest = false ;   //  Bug #969
	if ( (   ReindexColMap != 0  ) )  //  Bug #969
	  RunKluTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunKluTest = false ;   //  Bug #1403




	if ( RunKluTest && verbose) cout << " Testing KLU " << endl ; 
	
	if ( RunKluTest ) Errors = TestOtherClasses("Amesos_Klu",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 false,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTheseTests ) ;
#endif
      } else if ( AmesosClasses[i] == "Amesos_Lapack" ) {
	bool RunLapackTest = true;
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) && row_map.DistributedGlobal() ) 
	  RunLapackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunLapackTest = false ;   //  Bug #1403

	if ( RunLapackTest && verbose) cout << " Testing LAPACK " << endl ; 
	
	if ( RunLapackTest ) Errors = TestOtherClasses("Amesos_Lapack",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 false,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTheseTests ) ;
      } else if ( AmesosClasses[i] == "Amesos_Taucs" ) {
	bool RunTaucsTest = true;
	if ( ( ReindexRowMap != 0  || ReindexColMap != 0 ) ) 
	  RunTaucsTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunTaucsTest = false ;   //  Bug #1403
	//	if ( MissingADiagonal ) RunTaucsTest = false ; // Bug #1449
	//	if ( transpose ) RunTaucsTest = false ; // Bug #1579
	if ( a662_bus_out)  RunTaucsTest = false ; // Bug #1449
	if ( ! symmetric ) RunTaucsTest = false ; 
	if ( Khead ) RunTaucsTest = false ;   // Bug #1449

	if ( RunTaucsTest && verbose) cout << " Testing TAUCS " << endl ; 
	
	
	if ( RunTaucsTest ) Errors = TestOtherClasses("Amesos_Taucs",
						      EpetraMatrixType,
						      Amat, 
						      transpose, 
						      verbose, 
						      Levels, 
						      Rcond, 
						      RowMapEqualsColMap,
						      false, 
						      maxrelerror, 
						      maxrelresidual, 
						      NumTheseTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Pardiso" ) {
	bool RunPardisoTest = true;
	if ( ReindexRowMap != 0  || ReindexColMap != 0 ) // Bug #969 
	  RunPardisoTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) 
	  RunPardisoTest = false ;   //  Bug #1403
	if ( bcsstk04 ) RunPardisoTest = false ;   // Bug #1916
	if ( a662_bus_out )  RunPardisoTest = false ; // Bug #1916
	if ( transpose ) RunPardisoTest = false ;   // Bug #1992
	if ( MissingADiagonal ) RunPardisoTest = false ; // Bug #1916 
	if ( Khead ) RunPardisoTest = false ; // Bug #1916 
	if ( EpetraMatrixType == 1 )  RunPardisoTest = false ; // Bug #1994 
	if ( distribute )  RunPardisoTest = false ; // Bug #1995
	if ( RunPardisoTest && verbose) cout << " Testing PARDISO " << endl ; 
	if ( Amat->Comm().NumProc() > 1 ) RunPardisoTest = false ; 

	if ( RunPardisoTest ) Errors = TestOtherClasses("Amesos_Pardiso",
							 EpetraMatrixType,
							 Amat, 
							 transpose, 
							 verbose, 
							 Levels, 
							 Rcond, 
							 RowMapEqualsColMap,
							 false,
							 maxrelerror, 
							 maxrelresidual, 
							 NumTheseTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Mumps" ) {
	bool RunMumpsTest = true;
	if ( ( ReindexRowMap || ReindexColMap ) ) 
	  RunMumpsTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunMumpsTest = false ;   //  Bug #1403
	if ( MissingADiagonal ) RunMumpsTest = false ; // Bug #1435
	if ( distribute )  RunMumpsTest = false ; // Bug #
	if (  RunMumpsTest && verbose) cout << " Testing MUMPS " << endl ; 

	if ( RunMumpsTest ) Errors = TestOtherClasses("Amesos_Mumps",
							 EpetraMatrixType,
						       Amat, 
						       transpose, 
						       verbose, 
						       Levels, 
						       Rcond, 
						       RowMapEqualsColMap,
						      false,
						       maxrelerror, 
						       maxrelresidual, 
						       NumTheseTests ) ;

      } else if ( AmesosClasses[i] == "Amesos_Klu" ) {
	bool RunKluTest = true;
	if ( (   ReindexColMap != 0  ) )  //  Bug #969
	  RunKluTest = false ;   //  Bug #969

	//	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunKluTest = false ;   //  Bug #1403

	Teuchos::ParameterList ParamList;
	if ( ReindexRowMap != 0 )  ParamList.set( "Reindex", true );
	if ( ( RangeMapType != 0 || DomainMapType != 0 || distribute ) ) 
	  ParamList.set( "DontTrustMe", true );
#ifndef HAVE_AMESOS_EPETRAEXT
	if ( ( ReindexRowMap || ReindexColMap ) ) 
	  RunKluTest = false ;  
#endif
	if ( ImpcolB ) RunKluTest = false ;   // See bug #1928 

	if ( RunKluTest && verbose) cout << " Testing KLU " << endl ; 

	if ( RunKluTest ) Errors = TestKlu( Amat, 
					    EpetraMatrixType,
					    transpose, 
					    verbose, 
					    Levels,
					    Rcond, 
					    ParamList, 
					    RowMapEqualsColMap, 
					    false,
					    maxrelerror, 
					    maxrelresidual, 
					    NumTheseTests ) ;
  
	if ( Amat->Comm().MyPID() == 0 && Errors ) 
	  cout << " FAILURE in "  
	       << __FILE__ << "::"  << __LINE__
	       << " Amesos_Klu" 
	       << " EpetraMatrixType = " <<  EpetraMatrixType 
	       << " transpose = " <<  transpose 
	       << " symmetric = " <<  symmetric 
	       << " Levels = " <<  Levels 
	       << " Diagonal = " <<  Diagonal 
	       << " ReindexRowMap = " <<  ReindexRowMap 
	       << " ReindexColMap = " <<  ReindexColMap 
	       << " DomainMapType = " <<  DomainMapType 
	       << " RangeMapType = " <<  RangeMapType 
	       << " distribute = " <<  distribute 
	       << " filename = " <<  filename 
	       << " NumTheseTests = " <<  NumTheseTests 
	       << " Errors = " <<  Errors << endl ;  

      } else if ( AmesosClasses[i] == "Amesos_Superlu" ) {
	bool RunSuperluTest = true;
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && Amat->Comm().NumProc() > 1  )  //  Bug #969
	  RunSuperluTest = false ;   //  Bug #969
	if ( MissingADiagonal ) RunSuperluTest = false ; // Bug #1404
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunSuperluTest = false ;   //  Bug #1403
	if ( bcsstk04 && transpose ) RunSuperluTest = false ;  // Bug #1927 
	if ( a662_bus_out && transpose ) RunSuperluTest = false ;  // Bug #1927 
	if ( Khead ) RunSuperluTest= false ;  // Bug #1927 

	if ( RunSuperluTest ) {
	  if ( verbose) cout << " Testing SUPERLU " << endl ; 
	  Errors = TestOtherClasses("Amesos_Superlu",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     false,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTheseTests ) ;
	}
	if ( Amat->Comm().MyPID() == 0 && Errors ) 
	  cout << " FAILURE in " 
	       << __FILE__ << "::"  << __LINE__
	       << " Amesos_Superlu" 
	       << " EpetraMatrixType = " <<  EpetraMatrixType 
	       << " transpose = " <<  transpose 
	       << " symmetric = " <<  symmetric 
	       << " Levels = " <<  Levels 
	       << " Diagonal = " <<  Diagonal 
	       << " ReindexRowMap = " <<  ReindexRowMap 
	       << " ReindexColMap = " <<  ReindexColMap 
	       << " DomainMapType = " <<  DomainMapType 
	       << " RangeMapType = " <<  RangeMapType 
	       << " distribute = " <<  distribute 
	       << " filename = " <<  filename 
	       << " NumTheseTests = " <<  NumTheseTests 
	       << " Errors = " <<  Errors << endl ;  
      } else if ( AmesosClasses[i] == "Amesos_Pastix" ) {
	bool RunPastixTest = true;
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) && Amat->Comm().NumProc() > 1  )  //  Bug #969
	  RunPastixTest = false ;   //  Bug #969
	if ( MissingADiagonal ) RunPastixTest = false ; // Bug #1404
	//	RowMapEqualsColMap = false ; // Bug #1405 - this turns off the AddToDiag test 
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunPastixTest = false ;   //  Bug #1403

	if ( RunPastixTest ) {
	  if ( verbose) cout << " Testing Pastix " << endl ; 
	  Errors = TestOtherClasses("Amesos_Pastix",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     false,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTheseTests ) ;
	}

      } else if ( AmesosClasses[i] == "Amesos_Paraklete" ) {

	//  We only test reindexing on klu and paraklete
	if ( ( verbose  &&  ( ReIndex ) ) { 
	  
	  cout << __FILE__ << "::"  << __LINE__
	       << " Perhaps about to test " 
	       << AmesosClasses[i] << " "  
	       << " EpetraMatrixType = " <<  EpetraMatrixType 
	       << " transpose = " <<  transpose 
	       << " symmetric = " <<  symmetric 
	       << " Levels = " <<  Levels 
	       << " Diagonal = " <<  Diagonal 
	       << " ReindexRowMap = " <<  ReindexRowMap 
	       << " ReindexColMap = " <<  ReindexColMap 
	       << " DomainMapType = " <<  DomainMapType 
	       << " RangeMapType = " <<  RangeMapType 
	       << " distribute = " <<  distribute 
	       << " filename = " <<  filename 
	       << endl ;  
	}
	bool RunParakleteTest = true;
	if ( (   ReindexColMap != 0  ) )  //  Bug #969
	  RunParakleteTest = false ;   //  Bug #969

	//	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunParakleteTest = false ;   //  Bug #1403
	if ( RunParakleteTest && verbose) cout << " Testing PARAKLETE " << endl ; 

	Teuchos::ParameterList ParamList;
	if ( ReindexRowMap != 0 )  ParamList.set( "Reindex", true );
	if ( ( RangeMapType != 0 || DomainMapType != 0 || distribute ) ) 
	  ParamList.set( "DontTrustMe", true );
	if ( ! transpose ) RunParakleteTest = false ; // Bug #1953 
#ifndef HAVE_AMESOS_EPETRAEXT
	if ( ( ReindexRowMap || ReindexColMap ) ) 
	  RunParakleteTest = false ;  
#endif
	//	if ( ImpcolB ) RunParakleteTest = false ;   // See bug #1928 

	if ( RunParakleteTest ) {
	  if ( verbose) cout << " Testing Paraklete " << endl ; 
	  Errors = TestOtherClasses("Amesos_Paraklete",
							 EpetraMatrixType,
				     Amat, 
				     transpose, 
				     verbose, 
				     Levels, 
				     Rcond, 
				     RowMapEqualsColMap,
				     false,
				     maxrelerror, 
				     maxrelresidual, 
				     NumTheseTests ) ;
	}
      } else if ( AmesosClasses[i] == "Amesos_Dscpack" ) {
	//
	//  A quick sanity check - make sure symmetric is the same on all processes
	//
	const int sym_int = symmetric?0:1 ; 
	int sym_int_out = sym_int; 
	Amat->Comm().Broadcast( &sym_int_out, 1, 0 ) ; 
	assert( sym_int == sym_int_out ) ; 

	bool RunDscpackTest = true;
	if ( ! symmetric ) RunDscpackTest = false ; 
	if ( (  ReindexRowMap != 0 ||  ReindexColMap != 0  ) )  //  Bug #969
	  RunDscpackTest = false ;   //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunDscpackTest = false ;   //  Bug #1403
	if ( Khead ) RunDscpackTest = false ;   // Bug #1234


	if ( RunDscpackTest ) { 
	  if ( verbose) cout << " Testing DSCPACK " << endl ; 
    
	  Errors = TestOtherClasses("Amesos_Dscpack",
				    EpetraMatrixType,
				    Amat, 
				    transpose, 
				    verbose, 
				    Levels, 
				    Rcond, 
				    RowMapEqualsColMap,
				    false,
				    maxrelerror, 
				    maxrelresidual, 
				    NumTheseTests ) ;
	} 
      } else if ( AmesosClasses[i] == "Amesos_Superludist" ) {
	bool RunSuperludistTest = true;
	if ( transpose ) { 
	  RunSuperludistTest = false ;    // Bug #822
	}
	if ( ReindexRowMap || ReindexColMap ) RunSuperludistTest = false ;    //  Bug #969
	if ( ( RangeMapType != 0 || DomainMapType != 0 ) ) RunSuperludistTest = false ;   //  Bug #1403
	//	if ( MissingADiagonal ) RunSuperludistTest = false ; // Bug #1404 NOT
	if ( Khead ) RunSuperludistTest= false ;  // Bug #368
	if ( RunSuperludistTest ) { 
	  if ( verbose) cout << " Testing Superludist " << endl ; 
  
	  Errors = TestSuperludist(Amat, 
				   EpetraMatrixType,
				    transpose, 
				    verbose, 
				    Levels, 
				    Rcond, 
				    maxrelerror, 
				    maxrelresidual, 
				    NumTheseTests ) ;
	}
      }
      if ( Amat->Comm().MyPID() == 0 ) {
	if ( Errors || ( verbose && NumTheseTests > 0 ) ) { 
	  if ( Errors ) { 
	    cout << " FAILURE in " ; 
	  } else { 
	    cout << " NO FAILURE in " ; 
	  }
	
	  cout << __FILE__ << "::"  << __LINE__
	       << " " << AmesosClasses[i] << " "  
	       << " EpetraMatrixType = " <<  EpetraMatrixType 
	       << " transpose = " <<  transpose 
	       << " symmetric = " <<  symmetric 
	       << " Levels = " <<  Levels 
	       << " Diagonal = " <<  Diagonal 
	       << " ReindexRowMap = " <<  ReindexRowMap 
	       << " ReindexColMap = " <<  ReindexColMap 
	       << " DomainMapType = " <<  DomainMapType 
	       << " RangeMapType = " <<  RangeMapType 
	       << " distribute = " <<  distribute 
	       << " filename = " <<  filename 
	       << " NumTheseTests = " <<  NumTheseTests 
	       << " Errors = " <<  Errors << endl ;  
	}
      }
      errors += Errors ;
      NumTests += NumTheseTests ;
    }
  }
  
  if ( verbose) cout << " TestAllClasses errors = " << errors << endl ; 

  return errors;
}

