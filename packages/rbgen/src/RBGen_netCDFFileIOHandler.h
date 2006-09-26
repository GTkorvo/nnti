

#ifndef NETCDF_FILE_IO_HANDLER_H
#define NETCDF_FILE_IO_HANDLER_H


#include "RBGen_FileIOHandler.hpp"
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

// Forward declaration of Epetra_MultiVector class
class Epetra_MultiVector;

namespace RBGen {
 
  class netCDFFileIOHandler : public virtual FileIOHandler< Epetra_MultiVector > 
  {  
  public:
    //@{ @name Constructor/Destructor.
    
    //! Default constructor.
    netCDFFileIOHandler();
    
    //! Destructor.
    virtual ~netCDFFileIOHandler();
    
    //@}
    
    //@{ @name Initialization/Reset Methods
    
    void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params ) { params_ = params; }
    
    void Reset();
    
    //@}
    
    //@{ @name File Reading Methods
    
    //! Method for reading multiple files and putting them into an Epetra_MultiVector.
    Teuchos::RefCountPtr< Epetra_MultiVector > Read( const std::vector<std::string>& filenames );

    //@}

    //@{ @name Writing Methods

    //! Method for writing one Epetra_MultiVector into a file.
    void Write( Teuchos::RefCountPtr< const Epetra_MultiVector > MV, const std::string& filename );

    //@}

    //@{ @name Handler Status Methods

    //! Return initialized status of the handler
    bool isInitialized() const { return true; };

    //@}

  private:

    bool isInitialized_;
    int num_nodes, num_nod_var, len_string;
    char **var_name;
    Teuchos::RefCountPtr< Teuchos::ParameterList > params_;

    // Method for handling error from netCDF.
    void handle_error( int status ) {
      if (status != NC_NOERR) {
	fprintf(stderr,"%s\n", nc_strerror(status));
	exit(-1);
      }
    }
  };
  
} // namespace RBGen

#endif // NETCDF_FILE_IO_HANDLER_H

