
#include "RBGen_MatrixMarketFileIOHandler.h"
#include "RBGen_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "EpetraExt_mmio.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"	

#include "Teuchos_Utils.hpp"

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

  void MatrixMarketFileIOHandler::Initialize( const Teuchos::RefCountPtr<Teuchos::ParameterList>& params )
  {
    // Get the "File I/O" sublist.
    Teuchos::ParameterList& fileio_params = params->sublist( "File IO" );

    // Get the input path.
    in_path = "";
    if ( fileio_params.isParameter( "Data Input Path" ) ) {
      in_path = Teuchos::getParameter<std::string>( fileio_params, "Data Input Path" );
    }

    // Get the output path.
    out_path = "";
    if ( fileio_params.isParameter( "Data Output Path" ) ) {
      out_path = Teuchos::getParameter<std::string>( fileio_params, "Data Output Path" );
    }

    // This file i/o handler is now initialized.
    isInit = true;
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
        std::string temp_filename = in_path + filenames[i];
	handle = fopen(temp_filename.c_str(), "r");
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
	std::string curr_filename = in_path + filenames[i];
        int info = EpetraExt::MatrixMarketFileToMultiVector( curr_filename.c_str(), Map, fileMV );
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

      std::string temp_filename = out_path + filename;
      EpetraExt::MultiVectorToMatrixMarketFile( temp_filename.c_str(), *MV );

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
  }
  
} // namespace RBGen


