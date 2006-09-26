

#ifndef RBGEN_PREPROCESSOR_FACTORY_HPP
#define RBGEN_PREPROCESSOR_FACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace RBGen {

  template< class DataSetType >
  class Preprocessor;

  template< class DataSetType >
  class PreprocessorFactory {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    PreprocessorFactory() {};

    //! Destructor.
    virtual ~PreprocessorFactory() {};
    //@}

    //@{ @name Factory methods

    virtual Teuchos::RefCountPtr<Preprocessor<DataSetType> > create( const Teuchos::ParameterList& params ) = 0;
    
    //@}

  };
  
} // end of RBGen namespace

#endif // RBGEN_PREPROCESSOR_FACTORY_HPP
