

#ifndef RBGEN_METHOD_FACTORY_HPP
#define RBGEN_METHOD_FACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "RBGen_Method.hpp"

namespace RBGen {
 
  template< class DataSetType, class OperatorType > 
  class MethodFactory {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    MethodFactory() {};

    //! Destructor.
    virtual ~MethodFactory() {};
    //@}

    //@{ @name Factory methods

    virtual Teuchos::RefCountPtr<Method< DataSetType, OperatorType > > create( const Teuchos::ParameterList& params ) = 0;
    
    //@}

  };
  
} // end of RBGen namespace

#endif // RBGEN_METHOD_FACTORY_HPP
