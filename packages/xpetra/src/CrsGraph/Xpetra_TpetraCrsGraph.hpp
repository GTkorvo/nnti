#ifndef XPETRA_TPETRACRSGRAPH_HPP
#define XPETRA_TPETRACRSGRAPH_HPP

/* this file is automatically generated - do not edit (see script/tpetra.py) */

#include "Xpetra_TpetraConfigDefs.hpp"

#include "Xpetra_CrsGraph.hpp"

#include "Tpetra_CrsGraph.hpp"

#include "Xpetra_Utils.hpp"

#include "Xpetra_TpetraMap.hpp"

#include "Xpetra_TpetraImport.hpp"

#include "Xpetra_TpetraExport.hpp"

using Teuchos::ParameterList;

namespace Xpetra {

  // TODO: move that elsewhere
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP< const CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > toXpetra(RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > );

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP< const Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > toTpetra(const RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph);
  //

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class S, class LO, class GO, class N, class SpMatOps>
  class CrsMatrix;
#endif
  
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class TpetraCrsGraph
    : public CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

    // The following typedef is used by the XPETRA_DYNAMIC_CAST() macro.
    typedef TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> TpetraCrsGraphClass;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor specifying fixed number of entries for each row.
    TpetraCrsGraph(const RCP< const map_type > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const RCP< ParameterList > &params=null)
      : graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps >(toTpetra(rowMap), maxNumEntriesPerRow, toTpetra(pftype), params))) {  }

    //! Constructor specifying (possibly different) number of entries in each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const RCP< ParameterList > &params=null)
      : graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps >(toTpetra(rowMap), NumEntriesPerRowToAlloc, toTpetra(pftype), params))) {  }

    //! Constructor specifying column Map and fixed number of entries for each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const RCP< ParameterList > &params=null)
      : graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps >(toTpetra(rowMap), toTpetra(colMap), maxNumEntriesPerRow, toTpetra(pftype), params))) {  }

    //! Constructor specifying column Map and number of entries in each row.
    TpetraCrsGraph(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const RCP< ParameterList > &params=null)
      : graph_(Teuchos::rcp(new Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps >(toTpetra(rowMap), toTpetra(colMap), NumEntriesPerRowToAlloc, toTpetra(pftype), params))) {  }

    //! Destructor.
    virtual ~TpetraCrsGraph() {  }

    //@}

    //! @name Insertion/Removal Methods
    //@{

    //! Insert graph indices, using global IDs.
    void insertGlobalIndices(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &indices) { XPETRA_MONITOR("TpetraCrsGraph::insertGlobalIndices"); graph_->insertGlobalIndices(globalRow, indices); }

    //! Insert graph indices, using local IDs.
    void insertLocalIndices(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &indices) { XPETRA_MONITOR("TpetraCrsGraph::insertLocalIndices"); graph_->insertLocalIndices(localRow, indices); }

    //! Remove all graph indices from the specified local row.
    void removeLocalIndices(LocalOrdinal localRow) { XPETRA_MONITOR("TpetraCrsGraph::removeLocalIndices"); graph_->removeLocalIndices(localRow); }

    //@}

    //! @name Transformational Methods
    //@{

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params=null) { XPETRA_MONITOR("TpetraCrsGraph::fillComplete"); graph_->fillComplete(toTpetra(domainMap), toTpetra(rangeMap), params); }

    //! Signal that data entry is complete.
    void fillComplete(const RCP< ParameterList > &params=null) { XPETRA_MONITOR("TpetraCrsGraph::fillComplete"); graph_->fillComplete(params); }

    //@}

    //! @name Methods implementing RowGraph.
    //@{

    //! Returns the communicator.
    const RCP< const Comm< int > >  getComm() const { XPETRA_MONITOR("TpetraCrsGraph::getComm"); return graph_->getComm(); }

    //! Returns the Map that describes the row distribution in this graph.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const { XPETRA_MONITOR("TpetraCrsGraph::getRowMap"); return toXpetra(graph_->getRowMap()); }

    //! Returns the Map that describes the column distribution in this graph.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const { XPETRA_MONITOR("TpetraCrsGraph::getColMap"); return toXpetra(graph_->getColMap()); }

    //! Returns the Map associated with the domain of this graph.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const { XPETRA_MONITOR("TpetraCrsGraph::getDomainMap"); return toXpetra(graph_->getDomainMap()); }

    //! Returns the Map associated with the domain of this graph.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const { XPETRA_MONITOR("TpetraCrsGraph::getRangeMap"); return toXpetra(graph_->getRangeMap()); }

    //! Returns the importer associated with this graph.
    RCP< const Import< LocalOrdinal, GlobalOrdinal, Node > > getImporter() const { XPETRA_MONITOR("TpetraCrsGraph::getImporter"); return toXpetra(graph_->getImporter()); }

    //! Returns the number of global rows in the graph.
    global_size_t getGlobalNumRows() const { XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumRows"); return graph_->getGlobalNumRows(); }

    //! Returns the number of global columns in the graph.
    global_size_t getGlobalNumCols() const { XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumCols"); return graph_->getGlobalNumCols(); }

    //! Returns the number of graph rows owned on the calling node.
    size_t getNodeNumRows() const { XPETRA_MONITOR("TpetraCrsGraph::getNodeNumRows"); return graph_->getNodeNumRows(); }

    //! Returns the number of columns connected to the locally owned rows of this graph.
    size_t getNodeNumCols() const { XPETRA_MONITOR("TpetraCrsGraph::getNodeNumCols"); return graph_->getNodeNumCols(); }

    //! Returns the index base for global indices for this graph.
    GlobalOrdinal getIndexBase() const { XPETRA_MONITOR("TpetraCrsGraph::getIndexBase"); return graph_->getIndexBase(); }

    //! Returns the global number of entries in the graph.
    global_size_t getGlobalNumEntries() const { XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumEntries"); return graph_->getGlobalNumEntries(); }

    //! Returns the local number of entries in the graph.
    size_t getNodeNumEntries() const { XPETRA_MONITOR("TpetraCrsGraph::getNodeNumEntries"); return graph_->getNodeNumEntries(); }

    //! Returns the current number of entries on this node in the specified global row.
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { XPETRA_MONITOR("TpetraCrsGraph::getNumEntriesInGlobalRow"); return graph_->getNumEntriesInGlobalRow(globalRow); }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const { XPETRA_MONITOR("TpetraCrsGraph::getNumEntriesInLocalRow"); return graph_->getNumEntriesInLocalRow(localRow); }

    //! Returns the current number of allocated entries for this node in the specified global row .
    size_t getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const { XPETRA_MONITOR("TpetraCrsGraph::getNumAllocatedEntriesInGlobalRow"); return graph_->getNumAllocatedEntriesInGlobalRow(globalRow); }

    //! Returns the current number of allocated entries on this node in the specified local row.
    size_t getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const { XPETRA_MONITOR("TpetraCrsGraph::getNumAllocatedEntriesInLocalRow"); return graph_->getNumAllocatedEntriesInLocalRow(localRow); }

    //! Returns the number of global diagonal entries, based on global row/column index comparisons.
    global_size_t getGlobalNumDiags() const { XPETRA_MONITOR("TpetraCrsGraph::getGlobalNumDiags"); return graph_->getGlobalNumDiags(); }

    //! Returns the number of local diagonal entries, based on global row/column index comparisons.
    size_t getNodeNumDiags() const { XPETRA_MONITOR("TpetraCrsGraph::getNodeNumDiags"); return graph_->getNodeNumDiags(); }

    //! Returns the maximum number of entries across all rows/columns on all nodes.
    size_t getGlobalMaxNumRowEntries() const { XPETRA_MONITOR("TpetraCrsGraph::getGlobalMaxNumRowEntries"); return graph_->getGlobalMaxNumRowEntries(); }

    //! Returns the maximum number of entries across all rows/columns on this node.
    size_t getNodeMaxNumRowEntries() const { XPETRA_MONITOR("TpetraCrsGraph::getNodeMaxNumRowEntries"); return graph_->getNodeMaxNumRowEntries(); }

    //! Indicates whether the graph has a well-defined column map.
    bool hasColMap() const { XPETRA_MONITOR("TpetraCrsGraph::hasColMap"); return graph_->hasColMap(); }

    //! Indicates whether the graph is lower triangular.
    bool isLowerTriangular() const { XPETRA_MONITOR("TpetraCrsGraph::isLowerTriangular"); return graph_->isLowerTriangular(); }

    //! Indicates whether the graph is upper triangular.
    bool isUpperTriangular() const { XPETRA_MONITOR("TpetraCrsGraph::isUpperTriangular"); return graph_->isUpperTriangular(); }

    //! If graph indices are in the local range, this function returns true. Otherwise, this function returns false. */.
    bool isLocallyIndexed() const { XPETRA_MONITOR("TpetraCrsGraph::isLocallyIndexed"); return graph_->isLocallyIndexed(); }

    //! If graph indices are in the global range, this function returns true. Otherwise, this function returns false. */.
    bool isGloballyIndexed() const { XPETRA_MONITOR("TpetraCrsGraph::isGloballyIndexed"); return graph_->isGloballyIndexed(); }

    //! Returns true if fillComplete() has been called and the graph is in compute mode.
    bool isFillComplete() const { XPETRA_MONITOR("TpetraCrsGraph::isFillComplete"); return graph_->isFillComplete(); }

    //! Returns true if storage has been optimized.
    bool isStorageOptimized() const { XPETRA_MONITOR("TpetraCrsGraph::isStorageOptimized"); return graph_->isStorageOptimized(); }

    //! Extract a const, non-persisting view of global indices in a specified row of the graph.
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &Indices) const { XPETRA_MONITOR("TpetraCrsGraph::getGlobalRowView"); graph_->getGlobalRowView(GlobalRow, Indices); }

    //! Extract a const, non-persisting view of local indices in a specified row of the graph.
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices) const { XPETRA_MONITOR("TpetraCrsGraph::getLocalRowView"); graph_->getLocalRowView(LocalRow, indices); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { XPETRA_MONITOR("TpetraCrsGraph::description"); return graph_->description(); }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { XPETRA_MONITOR("TpetraCrsGraph::describe"); graph_->describe(out, verbLevel); }

    //@}

    //! @name Advanced methods, at increased risk of deprecation.
    //@{

    //! Get an ArrayRCP of the row-offsets.
    ArrayRCP< const size_t > getNodeRowPtrs() const { XPETRA_MONITOR("TpetraCrsGraph::getNodeRowPtrs"); return graph_->getNodeRowPtrs(); }

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const { XPETRA_MONITOR("TpetraCrsGraph::getMap"); return rcp( new TpetraMap< LocalOrdinal, GlobalOrdinal, Node >(graph_->getMap()) ); }

    //! Import.
    void doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source, 
                  const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM) { 
      XPETRA_MONITOR("TpetraCrsGraph::doImport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, source, tSource, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > v = tSource.getTpetra_CrsGraph();
      //graph_->doImport(toTpetraCrsGraph(source), *tImporter.getTpetra_Import(), toTpetra(CM));

      graph_->doImport(*v, toTpetra(importer), toTpetra(CM));
    }

    //! Export.
    void doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsGraph::doExport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, dest, tDest, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > v = tDest.getTpetra_CrsGraph();
      graph_->doExport(*v, toTpetra(importer), toTpetra(CM)); 

    }

    //! Import (using an Exporter).
    void doImport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsGraph::doImport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, source, tSource, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > v = tSource.getTpetra_CrsGraph();

      graph_->doImport(*v, toTpetra(exporter), toTpetra(CM));

    }

    //! Export (using an Importer).
    void doExport(const DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM) {
      XPETRA_MONITOR("TpetraCrsGraph::doExport");

      XPETRA_DYNAMIC_CAST(const TpetraCrsGraphClass, dest, tDest, "Xpetra::TpetraCrsGraph::doImport only accept Xpetra::TpetraCrsGraph as input arguments.");//TODO: remove and use toTpetra()
      RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > v = tDest.getTpetra_CrsGraph();

      graph_->doExport(*v, toTpetra(exporter), toTpetra(CM)); 

    }

    // @}

    //! @name Xpetra specific
    //@{

    //! TpetraCrsGraph constructor to wrap a Tpetra::CrsGraph object
    TpetraCrsGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > &graph) : graph_(graph) { }

    //! Get the underlying Tpetra graph
    RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > getTpetra_CrsGraph() const { return graph_; }
    
    //@}
    
  private:
    
    RCP< Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > graph_;

  }; // TpetraCrsGraph class

  // TODO: move that elsewhere
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP< const CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > toXpetra(RCP< const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > graph) { //TODO: return TpetraCrsGraph instead of CrsGraph
    // typedef TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsGraphClass;
    // XPETRA_RCP_DYNAMIC_CAST(const TpetraCrsGraphClass, graph, tGraph, "toTpetra");

    RCP<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tGraph = Teuchos::rcp_const_cast<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(graph); //TODO: can I avoid the const_cast ?
    return rcp ( new Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(tGraph) );
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP< const Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > toTpetra(const RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph) {
    typedef TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsGraphClass;
    XPETRA_RCP_DYNAMIC_CAST(const TpetraCrsGraphClass, graph, tpetraCrsGraph, "toTpetra");
    return tpetraCrsGraph->getTpetra_CrsGraph();
  }
  //

} // Xpetra namespace

#define XPETRA_TPETRACRSGRAPH_SHORT
#endif // XPETRA_TPETRACRSGRAPH_HPP
