#ifndef MUELU_NEEDS_DECL_HPP
#define MUELU_NEEDS_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <Teuchos_ParameterEntry.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_TwoKeyMap.hpp"

namespace MueLu {
  
  /*!
    @class Needs
    @brief Class that allows cross-factory communication of data needs.

    Maintains a list of 'Needs' for a given Level. For example, a restriction factory that
    transposes the tentative prolongator 'Needs' the prolongator factory to save this.

    The data is stored using a variable name and a pointer to the generating factory. The
    pointer to the generating factory is only used as "key". For the Needs class it doesn't
    matter if the given factory pointer is valid or not. So the NULL pointer can also be used.

    A reference counter keeps track of the storage and automatically frees the memory if
    the data is not needed any more. In the standard mode, the data first has to be requested
    by calling the Request function. Then the data can be set by calling Set.
    With Get the user can fetch data when needed. Release decrements the reference counter
    for the given variable.
  */
  class Needs : public BaseClass {

  private:

    UTILS::TwoKeyMap<std::string, const FactoryBase*, int>                     countTable_; //<! Stores number of outstanding requests for a need.
    UTILS::TwoKeyMap<std::string, const FactoryBase*, Teuchos::ParameterEntry> dataTable_;  //<! Stores data associated with a need.
    UTILS::TwoKeyMap<std::string, const FactoryBase*, bool>                    keepTable_;  //<! keep status of a variable

    //TODO: Key1 = const std::string?

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    Needs() ;

    virtual ~Needs() ;

    //@}

    //! @name Set
    //! @brief functions for setting data in data storage
    //@{

    //      void Set(const Key1 & key1, const Key2 & key2, const Value & entry) 

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void Set(const std::string & ename, const T & entry, const FactoryBase* factory) ; // Set

    //@}

    //! @name Request/Release data

    //@{

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string & ename, const FactoryBase* factory) ; //Request

    //! Decrement the storage counter.
    void Release(const std::string & ename, const FactoryBase* factory) ; //Release

    //@}

    //! @name Get data
    //@{

    //! @brief Get data without decrementing associated storage counter (i.e., read-only access)
    // Usage: Level->Get< RCP<Operator> >("A", factoryPtr)
    template <class T>
    const T & Get(const std::string & ename, const FactoryBase* factory) const ;

    //! @brief Get data without decrementing associated storage counter (i.e., read-only access)
    // Usage: Level->Get< RCP<Operator> >("A", factoryPtr)
    template <class T>
    T & Get(const std::string & ename, const FactoryBase* factory) ;

    //@}

    //! @name Permanent storage
    //@{

    //! @brief Keep variable 'ename' generated by 'factory'. This variable is not handled by the internal reference counter system
    void Keep(const std::string & ename, const FactoryBase* factory, bool keep = true) ;

    //! @brief returns true, if variable 'ename' generated by 'factory' is permanently stored
    //! returns false, if variable 'ename' is not requested or requested, but not for being permanently stored.
    bool IsKept(const std::string & ename, const FactoryBase* factory) const ;
    //@}

    //! @name Utilities.
    //@{

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string & ename, const FactoryBase* factory) const ;

    //! Test whether some data has been requested.
    // Note1: IsRequested() can be used inside of a factory to know if it worth it to build some data.
    // (if a factory is called, at least one ename of the factory is requested but maybe not all of them)
    // Note2: this tells nothing about whether the data actually exists.
    bool IsRequested(const std::string & ename, const FactoryBase* factory) const ;

    //! Test whether a factory is generating factory of a requested variable in Needs
    bool IsRequestedFactory(const FactoryBase* factory) ;

    //! returns how often factory is generating factory of requested variables in Needs
    //! used to decide whether Level.Release(factory) for releasing factory dependencies can
    //! be called safely
    int CountRequestedFactory(const FactoryBase* factory) ;

    //! Test whether a factory is among the generating factories of data that is already available
    bool IsAvailableFactory(const FactoryBase* factory) ;

    //! @brief Return the number of outstanding requests for a need.
    //!  Throws a <tt>Exceptions::RuntimeError</tt> exception if the need either hasn't been requested or hasn't been saved.
    int NumRequests(const std::string & ename, const FactoryBase* factory) const ;

    //! @name Helper functions
    //@{

    //! Returns a vector of strings containing all key names of requested variables
    std::vector<std::string> RequestedKeys() const ;

    std::vector<const FactoryBase*> RequestedFactories(const std::string & ename) const ;

    std::string GetType(const std::string & ename, const FactoryBase* factory) const ;

    //@}

    //! @name I/O Functions
    //@{

    //! Printing method
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;

    //@}

  private:

    //! Copy constructor
    Needs(const Needs & source) ;

  }; //class Needs

} //namespace MueLu

#define MUELU_NEEDS_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_NEEDS_DECL_HPP
