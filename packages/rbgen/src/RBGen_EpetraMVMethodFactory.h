

#ifndef RBGEN_EPETRAMV_METHOD_FACTORY_H
#define RBGEN_EPETRAMV_METHOD_FACTORY_H

#include "Teuchos_ParameterList.hpp"
#include "RBGen_MethodFactory.hpp"
#include "RBGen_ConfigDefs.h"

// Forward declaration of Epetra_Multivector
class Epetra_MultiVector;

namespace RBGen {
 
  class EpetraMVMethodFactory : public virtual MethodFactory<Epetra_MultiVector> {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    EpetraMVMethodFactory() {};

    //! Destructor.
    virtual ~EpetraMVMethodFactory() {};
    //@}

    //@{ @name Factory methods

    Teuchos::RefCountPtr<Method< Epetra_MultiVector > > create( const Teuchos::ParameterList& params );
    
    //@}

  };
  
} // end of RBGen namespace

#endif // RBGEN_EPETRAMV_METHOD_FACTORY_H
