
#ifndef RBGEN_EPETRAMV_FILEIO_FACTORY_HPP
#define RBGEN_EPETRAMV_FILEIO_FACTORY_HPP

#include "RBGen_FileIOFactory.hpp"
#include "RBGen_BurkardtFileIOHandler.h"
#include "RBGen_MatrixMarketFileIOHandler.h"
#include "RBGen_netCDFFileIOHandler.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include <string>

namespace RBGen {

  class EpetraMVFileIOFactory : public virtual FileIOFactory<Epetra_MultiVector> {
 
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    EpetraMVFileIOFactory();

    //! Destructor.
    virtual ~EpetraMVFileIOFactory() {};
    //@}

    //@{ @name Factory methods

    Teuchos::RefCountPtr< FileIOHandler< Epetra_MultiVector > > create( const Teuchos::ParameterList& params );

    //@}

  private:

    // Available file formats
    std::vector<std::string> file_formats;

  };



  template< class DataSetType > 
  Teuchos::RefCountPtr< FileIOHandler< DataSetType > >
  FileIOFactory< DataSetType >::create( const Teuchos::ParameterList& params )
  {
    if ( !params.isParameter( "File IO Type" ) ) {
      //  TO DO:  THROW EXCEPTION!!!! 
    }
    
    // Get the file format type
    std::string file_format = Teuchos::getParameter<std::string>( const_cast<Teuchos::ParameterList&>(params), 
								  "File IO Type" );

    Teuchos::RefCountPtr< FileIOHandler< DataSetType > > RBFileIO;
    
    // File input format based on Burkardt's input files
    if ( file_format == "Burkardt" )
      RBFileIO = Teuchos::rcp( new BurkardtFileIOHandler() );
    
    // File input format for netCDF files
    if ( file_format == "netCDF" )
      RBFileIO = Teuchos::rcp( new netCDFFileIOHandler() );
    
    // File input format for Matrix Market files
    if ( file_format == "Matrix Market" )
      RBFileIO = Teuchos::rcp( new MatrixMarketFileIOHandler() );
    //
    // Return the method created
    //
    return RBFileIO;
  }  

} // end of RBGen namespace

#endif
