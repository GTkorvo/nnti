#ifndef MUELU_HIERARCHYHELPERS_DECL_HPP
#define MUELU_HIERARCHYHELPERS_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

//! An exception safe way to call the method 'Level::SetFactoryManager()'
class SetFactoryManager {

public:

  //@{

  //!
  SetFactoryManager(Level & level, const RCP<const FactoryManagerBase> & factoryManager)
  : level_(level)
  ;

  //! Destructor.
  virtual ~SetFactoryManager() ;

  //@}

private:
  Level & level_;
  //TODO: const RCP<const FactoryManagerBase> prevFactoryManager_;, save & restore previous factoryManager instead of reset to Teuchos::null.
};

// On the first level, 'A' is not generated by a factory (but is user-defined)
// On other levels, we can use either NoFactory or the real 'A" factory.
class InternalFactoryManager : public FactoryManagerBase {

public:

  //!
  InternalFactoryManager(RCP<const FactoryManagerBase> & parentFactoryManager)
  : factoryManager_(parentFactoryManager), noFact_(NoFactory::getRCP())
  ;

  //! Destructor.
  virtual ~InternalFactoryManager() ;

  //! GetFactory
  const RCP<const FactoryBase> & GetFactory(const std::string & varName) const ;

  //! Clean
  void Clean() const ;

private:
  RCP<const FactoryManagerBase> factoryManager_;
  RCP<const FactoryBase> noFact_; //TODO: remove
};


template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
class TopRAPFactory : public TwoLevelFactoryBase {
#include "MueLu_UseShortNames.hpp"

public:

  TopRAPFactory(RCP<const FactoryManagerBase> parentFactoryManager)
  : factoryManager_(rcp( new InternalFactoryManager(parentFactoryManager))), PFact_(parentFactoryManager->GetFactory("P")), RFact_(parentFactoryManager->GetFactory("R")), AcFact_(parentFactoryManager->GetFactory("A"))
  ;

  virtual ~TopRAPFactory() ;

  void DeclareInput(Level & fineLevel, Level & coarseLevel) const ;

  void Build(Level & fineLevel, Level & coarseLevel) const ;

private:
  RCP<const FactoryManagerBase> factoryManager_;
  RCP<const FactoryBase> PFact_;
  RCP<const FactoryBase> RFact_;
  RCP<const FactoryBase> AcFact_;
};

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
class TopSmootherFactory : public SingleLevelFactoryBase { //TODO: inherit from SmootherFactoryBase ?
#include "MueLu_UseShortNames.hpp"

public:

  TopSmootherFactory(RCP<const FactoryManagerBase> parentFactoryManager, const std::string & varName)
  : factoryManager_(rcp( new InternalFactoryManager(parentFactoryManager))), smootherFact_(parentFactoryManager->GetFactory(varName))
  ;

  virtual ~TopSmootherFactory() ;

  void DeclareInput(Level & level) const ;

  void Build(Level & level) const ;

private:
  RCP<const FactoryManagerBase> factoryManager_;
  RCP<const FactoryBase> smootherFact_;
};

} // namespace MueLu

// TODO: remove 'RCP' for TopRAPFactory::factoryManager_ and TopSmootherFactory::factoryManager_?

#define MUELU_HIERARCHY_HELPERS_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_HIERARCHYHELPERS_DECL_HPP
