//
//  SSTest 
//
//  usage: 
//     Solver InputMatrix MatrixType Special Numsolves Transpose MaxError MaxResid 
//     Where solver is:  SuperLU, SuperLUdist, UMFPACK, KUNDERT, SPOOLES, SPOOLESERIAL or AZTEC 
//     special is, at present, only used in SuperLU, where 0 means dgssv
//     and 1 means dgssvx 
//  examples:
//     SSTest.exe SPOOLES ImpcolA.rua 0 1 1 0 1e-10 1e-10 
//     source SmallTest.csh
//
//  output:  
//    SST.log (append) 
//    SST.summary (append) 
//
//  exits with 0 if test completed (does not imply that the test passed)
//  exits with -1 if command line options or file permissions are wrong 
//
  int SPOOLESmsglvl ;
  int SPOOLES_front_matrix ;
  int SPOOLES_pivoting ;
  double SPOOLES_tau ; 
  double SPOOLES_droptol ; 
  int SPOOLES_lookahead ; 
  int SuperLU_permc ; 


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#if (defined(SOLARIS) || defined(TFLOP) ) 
struct utsname { 
char sysname [ 257 ] ; 
char nodename [ 257 ] ; 
char release [ 257 ] ; 
char version [ 257 ] ; 
char machine [ 257 ] ; 
} ; 
extern "C" {
extern int uname ( struct utsname * ) ; 
}
#else
#include <sys/utsname.h>
#endif
#include <string.h>
#include <assert.h>
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "SparseDirectTimingVars.h"
#include "TestSolver.h"

#ifdef SOLARIS
#include <unistd.h>
#endif

#if 0
extern "C" {
#include "BridgeMPI.h"
}
#endif

//  #include "TSF.h"
//  using std::exception ;

main(int argc, char **argv)
{

  //
  //  The following are the values returned from the tester
  //
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumMpiProcs = Comm.NumProc(); 

  const int MAX_TOLERANCE_RATIO = 1000 ;
  int exit_value = 0 ; 
  const int MAXNAMELENGTH = 800;
  utsname uname_buf; 
  char timebuffer[MAXNAMELENGTH];
  string Sprogram ;
  if ( argc >1 ) Sprogram = argv[1] ;
  const int NUM_PARAMS = 9 ; 
  const int NUM_SUPERLU_PARAMS = NUM_PARAMS + 1 ; 
  const int NUM_SPOOLES_PARAMS = NUM_PARAMS + 6 ; 
  bool argc_ok = ( argc == NUM_PARAMS ) ; 
  if ( argc == NUM_SPOOLES_PARAMS && Sprogram == "SPOOLES" ) argc_ok = true ; 
  if ( argc == NUM_SPOOLES_PARAMS && Sprogram == "SPOOLESSERIAL" ) 
    argc_ok = true ; 
  
  if ( argc == NUM_SUPERLU_PARAMS && Sprogram == "SuperLU" ) 
    argc_ok = true ; 

//
//  The usage print should be a subroutine and printed everywhere 
//  that we find a problem with the command line arguments
//
  if ( ! argc_ok ) {
    if ( MyPID == 0 ) {
      cerr << " argc = " << argc << " Sprogam= " << Sprogram << 
	" SPOOLES? " << (int) (Sprogram=="SPOOLES") << endl ; 
      cerr << "Usage: " << argv[0] <<" SolverName InputMatrix special numsolves transpose maxerror maxresid" << endl ; 
      cerr << "    Solvername = UMFPACK, SuperLUdist, AZTEC. SPOOLES, SPOOLESSERIAL, KUNDER or SuperLU " << endl;
      cerr << "    InputMatrix must be a file in Harwell Boeing format"<< endl;
      cerr << "    special is, at presetn, only used in SuperLU " << endl ; 
      cerr << "      special == 0 means use DGSSV, special == 1 means use DGSSVX " << endl ; 
      cerr << "    numsolves = number of right hand sidess (<0 means MRHS, >1 means BRHS) " << endl ; 
      cerr << "    transpose = 1 means test A^T x = b instead of Ax = b" << endl ; 
      cerr << "    maxerror = maximum allowed error  < 0 == no check " << endl ; 
      cerr << "    maxresid = maximum allowed residual < 0 == no check" << endl ; 
      cerr << "    if maxerror == -2 and maxresid == -2, failure (hang or abort) is expected" << endl ; 
      cerr << "    if maxerror == 1e30 and maxresid == 1e30, the solver is expected to finish but prodcue incorrect results" << endl ; 
      
    }
#ifdef EPETRA_MPI
    MPI_Finalize();
#endif
    exit( -1 ) ; 
  }
 
  if ( MyPID == 0 ) {
    int uname_stat = uname( &uname_buf ) ; 
    
    time_t now = time( NULL ) ; 
    tm *localnow = localtime( &now ) ; 
    (void) strftime( timebuffer, MAXNAMELENGTH, "20%g%b%d@%H:%M", localnow ) ;
  }

  //
  //  Open SST.log and SST.summary
  //
  char *ShortOutputFileName = (char *) "SST.summary" ;
  char *LongOutputFileName = (char *) "SST.log" ;

  bool summary = MyPID == 0 ;
  bool verbose = MyPID == 0 ; 
  bool log = MyPID == 0 ; 
#ifdef SOLARIS
  //  log = false ;                                     // On Solaris mpich, the second ofstream.open fails 
#endif

  FILE *matrix_fd;
  ofstream summary_file;

  if ( verbose ) { 
    matrix_fd = fopen( argv[2], "r" ) ; 
    if ( matrix_fd == NULL ) {
      cerr << "Unable to open " << argv[2] << " for reading" << endl ; 
      exit_value = - 1; 
    } else {
      fclose( matrix_fd ) ; 
    }

    if ( log ) { 
      SparseDirectTimingVars::log_file.open( LongOutputFileName, ios::app ) ; 
      if ( SparseDirectTimingVars::log_file.fail() ) {
	cerr << "Unable to open " << LongOutputFileName << " for writing" << endl ; 
	exit_value = - 1; 
      }
    }

    if ( summary ) {
      summary_file.open( ShortOutputFileName, ios::app ) ; 
      if ( summary_file.fail() ) {
	cerr << "Unable to open " << ShortOutputFileName << " for writing" << endl ; 
	exit_value = - 1; 
      }
    }
  }

  //
  //  Check command line parameters
  //
  SparseSolverType SparseSolver ; 
  
  int MatType = atoi( argv[3] ) ; 
  int special = atoi( argv[4] ) ; 
  int numsolves = atoi( argv[5] ) ; 
  int transpose =  atoi( argv[6] ) ; 
  double maxerror = atof( argv[7] ) ;
  double maxresid = atof( argv[8] ) ;

  if ( Sprogram == "UMFPACK" ) 
    SparseSolver = UMFPACK ; 
  else if  ( Sprogram == "AZTEC" ) 
    SparseSolver = Aztec ; 
  else if  ( Sprogram == "KUNDERT" ) 
    SparseSolver = KUNDERT ; 
  else if  ( Sprogram == "SuperLUdist" ) 
    SparseSolver = SuperLUdist ; 
  else if  ( Sprogram == "SPOOLES" ) 
    SparseSolver = SPOOLES ; 
  else if  ( Sprogram == "SPOOLESSERIAL" ) 
    SparseSolver = SPOOLESSERIAL ; 
  else if  ( Sprogram == "SuperLU" ) 
    SparseSolver = SuperLU ; 
  else {
    if (verbose) cerr << "Unknown program: " << Sprogram << endl ; 
    exit_value = -1 ; 
  }

  //
  //  Note - the following are untested and they are also 
  //  unused.  We have not yet passed them to the SpoolesOO
  //  and SpoolesserialOO
  //
  SPOOLESmsglvl = 0 ;
  SPOOLES_front_matrix = 0;
  SPOOLES_pivoting = 0 ;
  SPOOLES_tau = 100.0 ; 
  SPOOLES_droptol = 0.001 ; 
  SPOOLES_lookahead = 0 ; 
  if ( argc == NUM_SPOOLES_PARAMS && 
       ( Sprogram == "SPOOLES" ||  Sprogram == "SPOOLESSERIAL" ) )
    { 
      SPOOLESmsglvl = atoi( argv[8] );
      SPOOLES_front_matrix = atoi( argv[9] );
      SPOOLES_pivoting = atoi( argv[10] );
      SPOOLES_tau = atof( argv[11] );
      SPOOLES_droptol = atof( argv[12] );
      SPOOLES_lookahead = atoi( argv[13] );
    }
  
  SuperLU_permc = 1 ;    // Set the default to MMD on A'*A
  if ( argc == NUM_SUPERLU_PARAMS && 
       ( Sprogram == "SuperLU" || Sprogram == "SuperLUdist" ) ) { 
      SuperLU_permc = atoi( argv[8] );
      assert( SuperLU_permc >= 0 && SuperLU_permc <= 3 ) ; 
  }




  const int MaxNumSolves = 3200 ; 
  if ( MatType < 0 || MatType > 1  ) { 
    if ( verbose ) 
      cerr << " MatType must be 0 or 1, is: " 
	<< MatType << endl ; 
    exit_value = -1 ; 
  }
  if ( special < 0 || special > 1000  ) { 
    if ( verbose ) 
      cerr << " Special must be 0 or 1, is: " 
	<< special << endl ; 
    exit_value = -1 ; 
  }
  if ( numsolves< -MaxNumSolves || numsolves > MaxNumSolves ) { 
    if ( verbose ) 
      cerr << "The number of solves must be between 0 and " << MaxNumSolves 
	<< " is: "
	  << numsolves << endl ; 
    exit_value = -1 ; 
  }
  if ( MatType == 1 && numsolves != 1 ) {
    if ( verbose ) {
      cerr << "MatType == 1 (distributed) is only supported for numsolves == 1.   " <<
	" Here numsolves = " << numsolves << endl ; 
    }
    exit_value = -1 ;     
  }
  if ( transpose< 0 ||  transpose > 1) { 
    if ( verbose ) 
      cerr << "transpose must be 0 (no trans) or 1" 
	<< ", it is: "
	  << transpose << endl ; 
    exit_value = -1 ; 
  }
  if ( transpose< 0 ||  transpose > 1) { 
    if ( verbose ) 
      cerr << "transpose must be 0 (no trans) or 1" 
	<< ", it is: "
	  << transpose << endl ; 
    exit_value = -1 ; 
  }
#if 1
  if ( transpose != 0 && SparseSolver == SPOOLESSERIAL ) { 
    if ( verbose ) 
      cerr << "Our use of SPOOLESSERIAL does not support transpose yet" << endl ;
    exit_value = -1 ; 
  }
#endif
  if ( transpose != 0 && SparseSolver == KUNDERT ) { 
    if ( verbose ) 
      cerr << "Our use of KUNDERT does not support transpose yet" << endl ;
    exit_value = -1 ; 
  }
  if ( transpose != 0 && SparseSolver == Aztec ) { 
    if ( verbose ) 
      cerr << "Our use of AZTEC does not support transpose yet" << endl ;
    exit_value = -1 ; 
  }
  if ( NumMpiProcs != 1 && MatType != 1 && SparseSolver == Aztec ) { 
    if ( verbose ) 
      cerr << "AZTEC accepts only distributed matrices on multiple processes" << endl ;
    exit_value = -1 ; 
  }
  if ( numsolves != 1 &&  SparseSolver != SuperLU && SparseSolver != SuperLUdist ) {
    if ( verbose ) 
      cerr << "Only SuperLU and SuperLUdist support MRHS and BRHS" << endl ;
    exit_value = -1 ; 
  }
    

#ifdef TFLOP
	char *hostname = "TFLOP";
	char *releasenum = "unknown";
#else
	char *hostname = uname_buf.nodename ; 
	char *releasenum = uname_buf.release;
#endif
	


 Comm.Broadcast( &exit_value, 1, 0 ) ; 

  if ( exit_value == 0 ) { 

    AMESOS_MatrixType MatrixType = AMESOS_Serial ; 
    if ( MatType == 1 ) MatrixType = AMESOS_Distributed ; 

    if ( log ) { 
      //
      //  Log time stamp and machine information 
      //

      SparseDirectTimingVars::log_file << endl << "TIMESTAMP:" << hostname << " " 
				       << argv[1] << " " << timebuffer 
				       << " BEGIN RUN" << endl ; 
      SparseDirectTimingVars::log_file << uname_buf.sysname << 
	hostname << releasenum << uname_buf.version << 
	  uname_buf.machine << endl ;
    }
    if (summary ) { 
      summary_file << endl << setw(12) << hostname
		   << setw(12) <<  argv[1] 
		   << " " << setw(-1) << timebuffer
		   << setw(15) << argv[2] << setw(6)
		   << MatType << " " 
		   << special << " " 
		   << NumMpiProcs <<  setw(6) 
		   << numsolves << setw(3) << transpose << setprecision(12) ;
      if ( maxresid == -2 && maxerror == -2 ) summary_file << "Failure OK" ; 
      flush( summary_file ) ; 
    }
    if (MyPID == 0 ) { 
      cerr << endl << setw(12) << hostname
		   << setw(12) <<  argv[1] 
		   << " " << setw(-1) << timebuffer
		   << setw(15) << argv[2] << setw(6)
		   << MatType << " " 
		   << special << " " 
		   << NumMpiProcs <<  setw(6) 
		   << numsolves << setw(3) << transpose << setprecision(12) ;
      if ( maxresid == -2 && maxerror == -2 ) cerr << "Failure OK" ; 
      flush( cerr ) ; 
    }
    //
    //  Perform the test
    //    
    SparseDirectTimingVars::log_file << SparseDirectTimingVars::SS_Result << endl ; 

    try { 

	if ( numsolves < 0 ) { 
	  TestMrhsSolver( Comm, argv[2], - numsolves, SparseSolver, (transpose==1), special ) ; 
	} else { 
	  if ( numsolves > 1 ) { 
	    TestMultiSolver( Comm, argv[2], numsolves, SparseSolver, (transpose==1), special ) ; 
	  } else {
	    TestSolver( Comm, argv[2], SparseSolver, (transpose==1), special, MatrixType ) ; 
	  }
	} 
      //
      //  Log time and memory estimates
      //    
      if ( log ) {
	SparseDirectTimingVars::log_file << SparseDirectTimingVars::SS_Result << endl ; 
	
	//
	//  Print a single line to the summary file (and a copy of same to 
	//  the log file (details_fd) then print a final line to the log 
	//  file.  
	//
	SparseDirectTimingVars::log_file << endl << "TIMESTAMP:" << hostname 
					 << argv[1] << timebuffer 
					 << " END RUN" << endl ; 
	
	SparseDirectTimingVars::log_file 
	  << setw(12) << hostname << setw(9) <<  argv[1] 
	  << " " << setw(-1) << timebuffer
	  << setw(15) << argv[2] << setw(6) <<  NumMpiProcs <<  setw(6) 
	  << special << " " 
	  << numsolves << setw(6)
	  << transpose << setprecision(12) ;
	SparseDirectTimingVars::SS_Result.PrintSummary(SparseDirectTimingVars::log_file) ;
	SparseDirectTimingVars::log_file << "SS_Result = " 
					 << SparseDirectTimingVars::SS_Result 
					 << endl ; 

      }
      if (summary ) { 
	SparseDirectTimingVars::SS_Result.PrintSummary(summary_file) ;
	SparseDirectTimingVars::SS_Result.PrintSummary(cerr) ;
	bool ErrorOK = maxerror <= -1 ||  
	  SparseDirectTimingVars::SS_Result.Get_Error() < maxerror ;
	bool ResidualOK = maxresid <= -1 ||  
	  SparseDirectTimingVars::SS_Result.Get_Residual() < maxresid ;
	if ( ErrorOK && ResidualOK ) summary_file << " OK" ; 
	if ( ErrorOK && ResidualOK ) cerr << " OK" ; 
	if ( ! ErrorOK ) {
	  summary_file << " Error too large is: " << 
	    SparseDirectTimingVars::SS_Result.Get_Error() <<
	    " should be < " << maxerror  ; 
	  cerr << " Error too large is: " << 
	    SparseDirectTimingVars::SS_Result.Get_Error() <<
	    " should be < " << maxerror  ; 
	}
	//
	//  Here we check to see if the answer is better than we expect.
	//
	//  1e30 means that the code promises to complete but makes no 
	//  promise about the answer
	//  If maxerror is 1e30, we set it to 10, meaning that if the actual 
	//  answer is good to 10/ MAX_TOLERANCE_RATIO (presently 10/1000 = .01)
	//  we print a TOLERANCE is too large message.
	//
	if (maxerror == 1e30 ) maxerror = 10 ; 
	if (SparseDirectTimingVars::SS_Result.Get_Error() < maxerror / MAX_TOLERANCE_RATIO && 
	    maxerror > 1.1e-15 ) {
	  summary_file << " Error TOLERANCE is too large: " << 
	    SparseDirectTimingVars::SS_Result.Get_Error() <<
	    " is allowed to be " << maxerror  ; 
	  cerr << " Error TOLERANCE is too large: " << 
	    SparseDirectTimingVars::SS_Result.Get_Error() <<
	    " is allowed to be " << maxerror  ; 
	}
	if ( ! ResidualOK ) {
	  summary_file << " Residual too large is:" <<
	    SparseDirectTimingVars::SS_Result.Get_Residual() <<
	    " should be < " << maxresid  ; 
	  cerr << " Residual too large is:" <<
	    SparseDirectTimingVars::SS_Result.Get_Residual() <<
	    " should be < " << maxresid  ; 
	}

	if (maxresid == 1e30 ) maxresid = 10 ; 
	if (SparseDirectTimingVars::SS_Result.Get_Residual() < maxresid / MAX_TOLERANCE_RATIO && 
	    maxresid > 1.1e-15 ) {
	  summary_file << " Residual TOLERANCE is too large: " << 
	    SparseDirectTimingVars::SS_Result.Get_Residual() <<
	    " is allowed to be " << maxresid  ; 
	  cerr << " Residual TOLERANCE is too large: " << 
	    SparseDirectTimingVars::SS_Result.Get_Residual() <<
	    " is allowed to be " << maxresid  ; 
	}

      flush( summary_file ) ; 
      cerr << endl ; // Atlantis won't print anything without this.
      flush( cerr ) ; 
      }
    }
    catch(string errormsg)
      {
	if ( summary ) { summary_file << errormsg ; } 
	if ( log ) SparseDirectTimingVars::log_file << errormsg ; 
	if ( verbose || summary ) cerr << errormsg << endl;
      }

  }

  if ( summary )   summary_file.close() ; 
  if ( log )   SparseDirectTimingVars::log_file.close() ; 

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  exit( exit_value ) ; 
}

  

