#ifndef MUELU_NOFACTORY_DECL_HPP
#define MUELU_NOFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <Teuchos_Assert.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryBase.hpp"

namespace MueLu {

  class Level;
  
  /*!
    @class NoFactory class.
    @brief NoFactory that is used for data stored in level class for that no generating factory is available/necessary.

    Use Singleton pattern.
  */
  class NoFactory : public FactoryBase {

    //! Constructor.
    NoFactory() ;

  public:

    //! Destructor.
    virtual ~NoFactory() ;

    //! Implementation of FactoryBase interface
    //@{
    
    //!
    void CallBuild(Level & requestedLevel) const ;

    //!
    void CallDeclareInput(Level & requestedLevel) const ;

    //@}
    
    //! Static Get() functions
    //@{
    
    //! 
    static const RCP<const NoFactory> getRCP() ;

    //! 
    static const NoFactory* get() ;

    //@}

  private:
    static RCP<const NoFactory> noFactory_; // static NoFactory instance for user defined "factories"

  }; // class NoFactory

} // namespace MueLu

#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_NOFACTORY_DECL_HPP
