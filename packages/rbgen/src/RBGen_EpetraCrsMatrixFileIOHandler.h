#ifndef EPETRA_CRSMATRIX_FILE_IO_HANDLER_H
#define EPETRA_CRSMATRIX_FILE_IO_HANDLER_H


#include "RBGen_FileIOHandler.hpp"

// Forward declaration of Epetra_MultiVector class
class Epetra_CrsMatrix;

namespace RBGen {
  
  //! FileIOHandler for reading EpetraCrsMatrix data from a file using EpetraExt.
  class EpetraCrsMatrixFileIOHandler : public virtual FileIOHandler< Epetra_CrsMatrix > {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    EpetraCrsMatrixFileIOHandler();

    //! Destructor.
    virtual ~EpetraCrsMatrixFileIOHandler() {};

    //@}

    //! @name Initialization/Reset Methods
    //@{

    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params );

    void Reset() { isInit = false; };

    //@}

    //! @name File Reading Methods
    //@{

    //! Method for reading a file and constructing an Epetra_CrsMatrix.
    Teuchos::RCP<Epetra_CrsMatrix> Read( const std::vector<std::string>& filenames );

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one Epetra_CrsMatrix into a file using the same type as was.
    void Write( Teuchos::RCP<const Epetra_CrsMatrix> MTX, const std::string& filename );

    //@}
    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    bool isInitialized() const { return isInit; };

    //@}

  private:

    // Whether or not we know the file format.
    bool isInit;

    // File input / output paths
    std::string in_path, out_path;

    // ParameterList that this file handler was initialized with.
    Teuchos::RCP< Teuchos::ParameterList > params_;

  };
  
} // namespace RBGen

#endif // EPETRA_CRSMATRIX_FILE_IO_HANDLER_H

