
#include "RBGen_MatrixMarketFileIOHandler.h"

#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "EpetraExt_mmio.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"	

#include "Teuchos_Utils.hpp"

#include <stdio.h>
#include <stdlib.h>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


namespace RBGen {
  
  MatrixMarketFileIOHandler::MatrixMarketFileIOHandler()
    : num_nodes(0), isInit(false)
  {
  }

  Teuchos::RefCountPtr<Epetra_MultiVector> MatrixMarketFileIOHandler::Read( const std::vector<std::string>& filenames )
  {

    Teuchos::RefCountPtr<Epetra_MultiVector> newMV;

    if (isInit) {

#ifdef EPETRA_MPI
      Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
      Epetra_SerialComm comm;
#endif

      int i, rows = 0;
      int num_vecs = 0;
      int num_files = filenames.size();
      std::vector<int> cols(num_files,0);
      FILE * handle = 0;
      for (i=0; i<num_files; ++i) {
	int info = 0, rows_i = 0;

	// Open the data file
	handle = fopen(filenames[i].c_str(), "r");
	if (handle == 0) {
	  // TO DO:  THROW EXCEPTION!
	}

	// Get the array dimensions
	info = EpetraExt::mm_read_mtx_array_size( handle, &rows_i, &cols[i] );
	if (info != 0) {
	  // TO DO:  THROW EXCEPTION!
	}
	if (i==0) {
	  rows = rows_i;  // Get the number of rows from the first file
	}
	else {
	  // Check to make sure the number of rows is the same.
	  if (rows_i != rows) {
	    // TO DO:  THROW EXCEPTION!
	  }
	}	 
	// Add the number of columns up.
	num_vecs += cols[i];
	
	// Close the data file
	fclose( handle );
      }

      // Create the map and full multivector.
      Epetra_Map Map( rows, 0, comm );
      newMV = Teuchos::rcp( new Epetra_MultiVector( Map, num_vecs ) );

      // Create a pointer to a multivector 
      int col_ptr = 0;
      Epetra_MultiVector* fileMV = 0;

      for ( i=0; i<num_files; i++ ) {
	//
	//  Read in Epetra_MultiVector from file.
	//
	int info = EpetraExt::MatrixMarketFileToMultiVector( filenames[i].c_str(), Map, fileMV );
	if (info != 0) {
	  // TO DO:  THROW EXCEPTION!
	}
	//
	//  Get a view of the multivector columns.
	//
	Epetra_MultiVector subMV( View, *newMV, col_ptr, cols[i] );
	// 
	//  Put the multivector read in from the file into this subview.
	//
	subMV.Update( 1.0, *fileMV, 0.0 );
	//
	//  Clean up the multivector
	//
	if (fileMV) { delete fileMV; fileMV=0; }
      }

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
    // Return.
    return newMV;
  }
  
  void MatrixMarketFileIOHandler::Write( Teuchos::RefCountPtr<const Epetra_MultiVector> MV, const std::string& filename )
  {
    if (isInit) {

      EpetraExt::MultiVectorToMatrixMarketFile( filename.c_str(), *MV );

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
  }
  
} // namespace RBGen


