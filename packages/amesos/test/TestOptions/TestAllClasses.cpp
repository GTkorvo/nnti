#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "TestAllClasses.h"
#include "TestOtherClasses.h"
#include "TestSuperludist.h"
 
int TestAllClasses(Epetra_CrsMatrix *& Amat, 
		   bool transpose, 
		   bool verbose, 
		   bool symmetric, 
		   int Levels,
		   const double Rcond,
		   double &maxrelerror, 
		   double &maxrelresidual,
		   int &NumTests) {

  int errors = 0 ;

  if ( transpose ) { 
    if ( verbose ) cout << "Superludist not test with transpose enabled " << endl ; 
  } else {
    if ( verbose) cout << " Testing Superludist " << endl ; 
  
    errors += TestSuperludist(Amat, 
			      transpose, 
			      verbose, 
			      Levels, 
			      Rcond, 
			      maxrelerror, 
			      maxrelresidual, 
			      NumTests ) ;
  }
  
  if ( verbose) cout << " Testing UMFPACK " << endl ; 

  errors += TestOtherClasses("Amesos_Umfpack",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;

  //
  //  A quick sanity check - make you symmetric is the same on all processes
  //
  const int sym_int = symmetric?0:1 ; 
  int sym_int_out = sym_int; 
  Amat->Comm().Broadcast( &sym_int_out, 1, 0 ) ; 
  assert( sym_int == sym_int_out ) ; 

  if ( symmetric ) { 
    if ( verbose) cout << " Testing DSCPACK " << endl ; 
    
    errors += TestOtherClasses("Amesos_Dscpack",
			       Amat, 
			       transpose, 
			       verbose, 
			       Levels, 
			       Rcond, 
			       maxrelerror, 
			       maxrelresidual, 
			       NumTests ) ;
  } else {
    if ( verbose ) cout << " DSCPACK not tested on unsymmetric matrices " 
			<< endl ; 
  }
    
  if ( verbose) cout << " Testing SCALAPACK " << endl ; 

  errors += TestOtherClasses("Amesos_Scalapack",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  

  if ( verbose) cout << " Testing KLU " << endl ; 

  errors += TestOtherClasses("Amesos_Klu",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;
  
  if ( verbose) cout << " Testing MUMPS " << endl ; 

  errors += TestOtherClasses("Amesos_Mumps",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 
			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;

  if ( verbose) cout << " Testing SUPERLU " << endl ; 

  errors += TestOtherClasses("Amesos_Superlu",
			     Amat, 
			     transpose, 
			     verbose, 
			     Levels, 
			     Rcond, 

			     maxrelerror, 
			     maxrelresidual, 
			     NumTests ) ;

  if ( verbose) cout << " TestAllClasses errors = " << errors << endl ; 

  return errors;
}

