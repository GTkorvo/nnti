#include "RBGen_EpetraMVFileIOFactory.h"

namespace RBGen {

  EpetraMVFileIOFactory::EpetraMVFileIOFactory() 
  { 
     // Insert the acceptable input file types for this factory
     file_formats.push_back("Burkardt");
#ifdef HAVE_RBGEN_NETCDF
     file_formats.push_back("NetCDF");
#endif
#ifdef HAVE_RBGEN_EPETRAEXT
     file_formats.push_back("Matrix Market");
#endif
  }

  Teuchos::RefCountPtr< FileIOHandler< Epetra_MultiVector > >
  EpetraMVFileIOFactory::create( const Teuchos::ParameterList& params )
  {
    // See if the "File I/O" sublist exists
    if ( !params.isSublist( "File IO" ) ) {
      //  TO DO:  THROW EXCEPTION!!!!
    }

    // Get the "File I/O" sublist.
    const Teuchos::ParameterList& fileio_params = params.sublist( "File IO" );

    // Get the file format type
    std::string file_format = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(fileio_params),
                                                                  "Type" );

    Teuchos::RefCountPtr< FileIOHandler< Epetra_MultiVector > > RBFileIO;

    // File input format based on Burkardt's input files
    if ( file_format == "Burkardt" ) {
      RBFileIO = Teuchos::rcp( new BurkardtFileIOHandler() );
    } else
    // File input format for NetCDF files
#ifdef HAVE_RBGEN_NETCDF
    if ( file_format == "NetCDF" ) {
      RBFileIO = Teuchos::rcp( new NetCDFFileIOHandler() );
    } else
#endif
    // File input format for Matrix Market files
#ifdef HAVE_RBGEN_EPETRAEXT
    if ( file_format == "Matrix Market" ) {
      RBFileIO = Teuchos::rcp( new MatrixMarketFileIOHandler() );
    } else 
#endif
    {
    // Throw an exception because the format type is not recognized by this factory
    }
    //
    // Return the method created
    //
    return RBFileIO;
  }
  
} // end of RBGen namespace

