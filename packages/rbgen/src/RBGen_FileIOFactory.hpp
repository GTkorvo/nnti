#ifndef RBGEN_FILEIO_FACTORY_HPP
#define RBGEN_FILEIO_FACTORY_HPP

// Forward declarations for Teuchos.
namespace Teuchos {
  class ParameterList;

  template <class T>
  class RefCountPtr;
}

namespace RBGen {

  // Forward declarations for FileIOHandler.
  template< class DataSetType >
  class FileIOHandler;

  //! Abstract factory for instantiating FileIOHandler concrete classes.
  /*!
   */
  template< class DataSetType > 
  class FileIOFactory {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    FileIOFactory() {};

    //! Destructor.
    virtual ~FileIOFactory() {};
    //@}

    //! @name Factory methods
    //@{

    virtual Teuchos::RefCountPtr< FileIOHandler< DataSetType > > create( const Teuchos::ParameterList& params ) = 0;

    //@}

  };

} // end of RBGen namespace

#endif // RBGEN_FILEIO_FACTORY_HPP
