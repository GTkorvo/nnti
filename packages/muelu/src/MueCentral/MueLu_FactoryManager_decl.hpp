#ifndef MUELU_FACTORYMANAGER_DECL_HPP
#define MUELU_FACTORYMANAGER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_UCAggregationFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"

namespace MueLu {

  /*!
    @class FactoryManager
    @brief This class specifies the default factory that should generate some data on a Level if the data does not exist and
    the generating factory has not been specified.
    
    Consider the following example.

    @code
      RCP<SingleLevelFactory> Afact;
      Level currentLevel;
      RCP<Operator> thisLevelA;
      thisLevelA = currentLevel.Get<Operator>("A",Afact.get());
    @endcode

    @todo If Afact is null (actually, Teuchos::null), then the FactoryManager associated with currentLevel will determine whether a default factory has
    been specified for creating A.  If "yes", then that factory will be called, A will be stored in currentLevel, and an RCP will be returned by
    the Get call.  If "no", then the FactoryManager will <b>throw an exception indicating that it does not know how to generate A</b>.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class FactoryManager : public FactoryManagerBase {
#undef MUELU_FACTORYMANAGER_SHORT
#include "MueLu_UseShortNames.hpp"
      
  public:

    //@{

    /*! @brief Constructor.

            @param[in] PFact Factory to generate the prolongation operator.
            @param[in] RFact Factory to generate the restriction operator.
            @param[in] AcFact Factory to generate the coarse grid operator A.
    */
    FactoryManager(const RCP<const FactoryBase> PFact = Teuchos::null, const RCP<const FactoryBase> RFact = Teuchos::null, const RCP<const FactoryBase> AcFact = Teuchos::null);
    
    //! Constructor used by HierarchyFactory (temporary, will be removed)
    FactoryManager(const std::map<std::string, RCP<const FactoryBase> >& factoryTable) {
      factoryTable_ = factoryTable;
      SetIgnoreUserData(false); // set IgnorUserData flag to false (default behaviour)
    }

    //! Destructor.
    virtual ~FactoryManager();

    //@}

    //@{ Get/Set functions.

    /*! @brief Set Factory

        Register the factory that should generate data if said factory is not specified in the request.

        @param[in] name of variable
        @param[in] factory that generates the data
    */
    void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory);

    /*! @brief Get factory associated with a particular data name.

       @param[in] varName name of variable.

    */
    const RCP<const FactoryBase> & GetFactory(const std::string & varName) const;

    //!
    const RCP<const FactoryBase> & GetDefaultFactory(const std::string & varName) const;

    //@}

    void Clean() const;

  private:
    
    //! @name Helper functions
    //@{

    /*! Add a factory to the default factory list and return it. This helper function is used by GetDefaultFactory()

     @todo TODO factory->setObjectLabel("Default " + varName + "Factory");
    */

    const RCP<const FactoryBase> & SetAndReturnDefaultFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const;

    //! Test if factoryTable_[varName] exists
    static bool IsAvailable(const std::string & varName, const std::map<std::string, RCP<const FactoryBase> > & factoryTable);
    //@}

    /*! @brief User-defined factories.

      Note: we distinguish 'user defined factory' and 'default factory' to allow the deallocation of default factories separately.
    */
    std::map<std::string, RCP<const FactoryBase> > factoryTable_;

    /*! @brief Table that holds default factories.

      -# We distinguish 'user defined factory' and 'default factory' to allow the deallocation of default factories separately.
      -# <tt>defaultFactoryTable_</tt> is mutable because default factories are only added to the list when they are requested
      to avoid allocation of unused factories.
    */
    mutable 
    std::map<std::string, RCP<const FactoryBase> > defaultFactoryTable_;
    
  }; // class

} // namespace MueLu

#define MUELU_FACTORYMANAGER_SHORT
#endif // MUELU_FACTORYMANAGER_DECL_HPP
