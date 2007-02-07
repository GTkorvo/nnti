
#include "RBGen_EpetraCrsMatrixFileIOHandler.h"
#include "RBGen_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "EpetraExt_RowMatrixOut.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


namespace RBGen {
  
  EpetraCrsMatrixFileIOHandler::EpetraCrsMatrixFileIOHandler()
    : isInit(false)
  {
  }

  void EpetraCrsMatrixFileIOHandler::Initialize( const Teuchos::RefCountPtr<Teuchos::ParameterList>& params )
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

  Teuchos::RefCountPtr<Epetra_CrsMatrix> EpetraCrsMatrixFileIOHandler::Read( const std::vector<std::string>& filenames )
  {

    Teuchos::RefCountPtr<Epetra_CrsMatrix> newMTX;

    if (isInit) {

#ifdef EPETRA_MPI
      Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
      Epetra_SerialComm comm;
#endif

      if (filenames.size() > 1) {
	  // TO DO:  THROW EXCEPTION!
      }
	
      // Open the data file
      std::string temp_filename = in_path + filenames[0];

      // Create a null pointer to the Epetra_Map
      Teuchos::RefCountPtr<Epetra_Map> Map;

      // Read in the matrix from file
      EpetraExt::readEpetraLinearSystem( temp_filename, comm, &newMTX, &Map );

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
    // Return.
    return newMTX;
  }
  
  void EpetraCrsMatrixFileIOHandler::Write( Teuchos::RefCountPtr<const Epetra_CrsMatrix> MTX, const std::string& filename )
  {
    if (isInit) {

      std::string temp_filename = out_path + filename;
      EpetraExt::RowMatrixToMatrixMarketFile( temp_filename.c_str(), *MTX );

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
  }
  
} // namespace RBGen


