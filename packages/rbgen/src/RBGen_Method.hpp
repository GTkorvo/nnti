

#ifndef RBGEN_METHOD_HPP
#define RBGEN_METHOD_HPP

#include "RBGen_FileIOHandler.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

namespace RBGen {
  
  template<class DataSetType, class OperatorType>
  class Method {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    Method() {};


    //! Destructor.
    virtual ~Method() {};
    //@}

    //@{ @name Computation Methods

    virtual void computeBasis() = 0;

    virtual void updateBasis( const Teuchos::RefCountPtr< DataSetType >& update_ss ) = 0;
    
    //@}

    //@{ @name Get Methods
    
    //! Get the basis computed by the reduced basis method.
    virtual Teuchos::RefCountPtr<const DataSetType> getBasis() const = 0;

    //! Returns the computational time taken to compute the reduced basis.
    virtual double getCompTime() const = 0;

    //@}
    
    //@{ @name Set Methods

    //! Initialize the method with the given parameter list and snapshot set.
    virtual void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                             const Teuchos::RefCountPtr< DataSetType >& ss,
			     const Teuchos::RefCountPtr< RBGen::FileIOHandler< OperatorType > >& fileio = Teuchos::null ) = 0;
    
    //! Reset the snapshot set used to compute the reduced basis.
    virtual void Reset( const Teuchos::RefCountPtr< DataSetType >& new_ss ) = 0;

    //@}

    //@{ @name Status Methods

    virtual bool isInitialized() = 0;

    //@}
  };
  
} // end of RBGen namespace

#endif // RBGEN_METHOD_HPP
