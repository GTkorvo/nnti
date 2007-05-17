#ifndef RBGEN_FILE_IO_HANDLER_HPP
#define RBGEN_FILE_IO_HANDLER_HPP

#include "RBGen_ConfigDefs.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace RBGen {

  //! Abstract base class for reading datasets from files.
  /*!
   */
  template<class DataSetType>  
  class FileIOHandler {  
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOHandler() {};

    //! Destructor.
    virtual ~FileIOHandler() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize file reader using 
    virtual void Initialize( const Teuchos::RefCountPtr<Teuchos::ParameterList>& params ) = 0;
    
    void Reset() {};
    //@}
    
    //! @name File Reading Methods
    //@{
    
    //! Method for reading multiple files and putting them into an data set.
    virtual Teuchos::RefCountPtr< DataSetType > Read( const std::vector<std::string>& filenames ) = 0;

    //@}

    //! @name Writing Methods
    //@{

    //! Method for writing one data set into a file.
    virtual void Write( Teuchos::RefCountPtr< const DataSetType > MV, const std::string& filename ) = 0;

    //@}

    //! @name Handler Status Methods
    //@{

    //! Return initialized status of the handler
    virtual bool isInitialized() const = 0;

    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_FILE_IO_HANDLER_HPP
