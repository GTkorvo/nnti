#ifndef MUELU_UCAGGREGATIONFACTORY_DECL_HPP
#define MUELU_UCAGGREGATIONFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <Teuchos_Utils.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_LocalAggregationAlgorithm.hpp"
#include "MueLu_LeftoverAggregationAlgorithm.hpp"

namespace MueLu {

  /*!
    @class UCAggregationFactory class.
    @brief Factory for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.
  */

  /* Factory input:
     - a graph ("Graph") generated by GraphFact_

     Factory output:
     - aggregates ("Aggegates")

     Factory options:
     - 
     - 
     - 
  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class UCAggregationFactory : public SingleLevelFactoryBase {
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    UCAggregationFactory(RCP<FactoryBase> graphFact = Teuchos::null)
      : graphFact_(graphFact)
    ;

    //! Destructor.
    virtual ~UCAggregationFactory() ;
    //@}

    //! @name Set/get methods.
    //@{

    // Options algo1
    void SetOrdering(Ordering ordering) ;
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) ;
    Ordering GetOrdering() const ;
    int GetMaxNeighAlreadySelected() const ;

    // Options algo2
    void SetPhase3AggCreation(double phase3AggCreation) ;
    double GetPhase3AggCreation() const ;

    // Options shared algo1 and algo2
    void SetMinNodesPerAggregate(int minNodesPerAggregate) ;
    int GetMinNodesPerAggregate() const ;
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Build aggregates. */
    void Build(Level &currentLevel) const
    ;

    //@}

  private:

    //! Graph Factory
    RCP<FactoryBase> graphFact_;
 
    //! Algorithms
    LocalAggregationAlgorithm algo1_;
    LeftoverAggregationAlgorithm algo2_;

  }; // class UCAggregationFactory

} //namespace MueLu

//TODO: can be more generic:
// - allow to choose algo 
// - base class for algorithm and options forward to algorithm as parameter list

#define MUELU_UCAGGREGATIONFACTORY_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_UCAGGREGATIONFACTORY_DECL_HPP
