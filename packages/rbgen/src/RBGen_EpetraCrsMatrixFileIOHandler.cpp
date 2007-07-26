
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

  void EpetraCrsMatrixFileIOHandler::Initialize( const Teuchos::RCP<Teuchos::ParameterList>& params )
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

  Teuchos::RCP<Epetra_Operator> EpetraCrsMatrixFileIOHandler::Read( const std::vector<std::string>& filenames )
  {

    Teuchos::RCP<Epetra_CrsMatrix> newMTX;

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
      Teuchos::RCP<Epetra_Map> Map;

      // Read in the matrix from file
      EpetraExt::readEpetraLinearSystem( temp_filename, comm, &newMTX, &Map );

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
    // Return.
    return newMTX;
  }
  
  void EpetraCrsMatrixFileIOHandler::Write( const Teuchos::RCP<const Epetra_Operator>& MTX, const std::string& filename )
  {
    if (isInit) {

      std::string temp_filename = out_path + filename;
      EpetraExt::RowMatrixToMatrixMarketFile( temp_filename.c_str(), *(Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(MTX)) );

    }
    else {
      // TO DO:  THROW EXCEPTION!
    }      
  }
  
} // namespace RBGen


