
#ifndef RBGEN_EPETRAMV_FILEIO_FACTORY_HPP
#define RBGEN_EPETRAMV_FILEIO_FACTORY_HPP

#include "RBGen_FileIOFactory.hpp"
#include "RBGen_BurkardtFileIOHandler.h"
#include "RBGen_ConfigDefs.h"

#ifdef HAVE_RBGEN_EPETRAEXT
#include "RBGen_MatrixMarketFileIOHandler.h"
#endif

#ifdef HAVE_NETCDF
#include "RBGen_netCDFFileIOHandler.h"
#endif

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

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

} // end of RBGen namespace

#endif
