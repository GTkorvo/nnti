#include "Epetra_CrsMatrix.h"
 
int TestOtherClasses( const char* AmesosClass,
		     int EpetraMatrixType,
		      Epetra_CrsMatrix *& Amat, 
		      const bool transpose, 
		      const bool verbose, 
		      const int Levels,
		      const double Rcond,
		      bool RowMapEqualsColMap, 
		      double &maxrelerror, 
		      double &maxrelresidual,
		      int &NumTests ) ;

