// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_TPETRAMAP_HPP
#define XPETRA_TPETRAMAP_HPP

/* this file is automatically generated - do not edit (see script/tpetra.py) */

#include "Xpetra_TpetraConfigDefs.hpp"

#include <Tpetra_Map.hpp>

#include "Xpetra_Map.hpp"
#include "Xpetra_Utils.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  // TODO: move that elsewhere
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & toTpetra(const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &);

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > toTpetra(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &);

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > toXpetra(const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > &);
  //

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class TpetraMap
    : public virtual Map<LocalOrdinal,GlobalOrdinal,Node> {

  public:

    //! @name Constructors and destructor
    //@{

    //! Constructor with Tpetra-defined contiguous uniform distribution.
    TpetraMap(global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP< Node > &node=KokkosClassic::Details::getNode<Node>())
      : map_(Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalElements, indexBase, comm, toTpetra(lg), node))) {  }

    //! Constructor with a user-defined contiguous distribution.
    TpetraMap(global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=KokkosClassic::Details::getNode<Node>())
      : map_(Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalElements, numLocalElements, indexBase, comm, node))) {  }

    //! Constructor with user-defined arbitrary (possibly noncontiguous) distribution.
    TpetraMap(global_size_t numGlobalElements, const Teuchos::ArrayView< const GlobalOrdinal > &elementList, GlobalOrdinal indexBase, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node=KokkosClassic::Details::getNode<Node>())
      : map_(Teuchos::rcp(new Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node >(numGlobalElements, elementList, indexBase, comm, node))) {  }

    //! Destructor.
    ~TpetraMap() {  }

    //@}

    //! @name Attributes
    //@{

    //! The number of elements in this Map.
    global_size_t getGlobalNumElements() const { XPETRA_MONITOR("TpetraMap::getGlobalNumElements"); return map_->getGlobalNumElements(); }

    //! The number of elements belonging to the calling node.
    size_t getNodeNumElements() const { XPETRA_MONITOR("TpetraMap::getNodeNumElements"); return map_->getNodeNumElements(); }

    //! The index base for this Map.
    GlobalOrdinal getIndexBase() const { XPETRA_MONITOR("TpetraMap::getIndexBase"); return map_->getIndexBase(); }

    //! The minimum local index.
    LocalOrdinal getMinLocalIndex() const { XPETRA_MONITOR("TpetraMap::getMinLocalIndex"); return map_->getMinLocalIndex(); }

    //! The maximum local index on the calling process.
    LocalOrdinal getMaxLocalIndex() const { XPETRA_MONITOR("TpetraMap::getMaxLocalIndex"); return map_->getMaxLocalIndex(); }

    //! The minimum global index owned by the calling process.
    GlobalOrdinal getMinGlobalIndex() const { XPETRA_MONITOR("TpetraMap::getMinGlobalIndex"); return map_->getMinGlobalIndex(); }

    //! The maximum global index owned by the calling process.
    GlobalOrdinal getMaxGlobalIndex() const { XPETRA_MONITOR("TpetraMap::getMaxGlobalIndex"); return map_->getMaxGlobalIndex(); }

    //! The minimum global index over all processes in the communicator.
    GlobalOrdinal getMinAllGlobalIndex() const { XPETRA_MONITOR("TpetraMap::getMinAllGlobalIndex"); return map_->getMinAllGlobalIndex(); }

    //! The maximum global index over all processes in the communicator.
    GlobalOrdinal getMaxAllGlobalIndex() const { XPETRA_MONITOR("TpetraMap::getMaxAllGlobalIndex"); return map_->getMaxAllGlobalIndex(); }

    //! The local index corresponding to the given global index.
    LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const { XPETRA_MONITOR("TpetraMap::getLocalElement"); return map_->getLocalElement(globalIndex); }

    //! The global index corresponding to the given local index.
    GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const { XPETRA_MONITOR("TpetraMap::getGlobalElement"); return map_->getGlobalElement(localIndex); }

    //! Return the process IDs and corresponding local IDs for the given global IDs.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList, const Teuchos::ArrayView< LocalOrdinal > &LIDList) const { XPETRA_MONITOR("TpetraMap::getRemoteIndexList"); return toXpetra(map_->getRemoteIndexList(GIDList, nodeIDList, LIDList)); }

    //! Return the process IDs for the given global IDs.
    LookupStatus getRemoteIndexList(const Teuchos::ArrayView< const GlobalOrdinal > &GIDList, const Teuchos::ArrayView< int > &nodeIDList) const { XPETRA_MONITOR("TpetraMap::getRemoteIndexList"); return toXpetra(map_->getRemoteIndexList(GIDList, nodeIDList)); }

    //! Return a view of the global indices owned by this node.
    Teuchos::ArrayView< const GlobalOrdinal > getNodeElementList() const { XPETRA_MONITOR("TpetraMap::getNodeElementList"); return map_->getNodeElementList(); }

    //@}

    //! @name Boolean tests
    //@{

    //! True if the local index is valid for this Map on this node, else false.
    bool isNodeLocalElement(LocalOrdinal localIndex) const { XPETRA_MONITOR("TpetraMap::isNodeLocalElement"); return map_->isNodeLocalElement(localIndex); }

    //! True if the global index is found in this Map on this node, else false.
    bool isNodeGlobalElement(GlobalOrdinal globalIndex) const { XPETRA_MONITOR("TpetraMap::isNodeGlobalElement"); return map_->isNodeGlobalElement(globalIndex); }

    //! True if this Map is distributed contiguously, else false.
    bool isContiguous() const { XPETRA_MONITOR("TpetraMap::isContiguous"); return map_->isContiguous(); }

    //! Whether this Map is globally distributed or locally replicated.
    bool isDistributed() const { XPETRA_MONITOR("TpetraMap::isDistributed"); return map_->isDistributed(); }

    //! True if and only if map is compatible with this Map.
    bool isCompatible(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { XPETRA_MONITOR("TpetraMap::isCompatible"); return map_->isCompatible(toTpetra(map)); }

    //! True if and only if map is identical to this Map.
    bool isSameAs(const Map< LocalOrdinal, GlobalOrdinal, Node > &map) const { XPETRA_MONITOR("TpetraMap::isSameAs"); return map_->isSameAs(toTpetra(map)); }

    //@}

    //! @name
    //@{

    //! Get this Map's Comm object.
    const Teuchos::RCP< const Teuchos::Comm< int > >  getComm() const { XPETRA_MONITOR("TpetraMap::getComm"); return map_->getComm(); }

    //! Get this Map's Node object.
    const Teuchos::RCP< Node >  getNode() const { XPETRA_MONITOR("TpetraMap::getNode"); return map_->getNode(); }

    //@}

    //! @name
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const { XPETRA_MONITOR("TpetraMap::description"); return map_->description(); }

    //! Print this object with the given verbosity level to the given FancyOStream.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { XPETRA_MONITOR("TpetraMap::describe"); map_->describe(out, verbLevel); }

    RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > removeEmptyProcesses () const {
      return toXpetra(map_->removeEmptyProcesses());
    }
    RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > replaceCommWithSubset (const Teuchos::RCP<const Teuchos::Comm<int> >& newComm) const {
      return toXpetra(map_->replaceCommWithSubset(newComm));
    }

    //@}

    //! @name Xpetra specific
    //@{

    //! TpetraMap constructor to wrap a Tpetra::Map object
    TpetraMap(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node > > &map)
      : map_(map) { }

    //! Get the library used by this object (Tpetra or Epetra?)
    UnderlyingLib lib() const { return Xpetra::UseTpetra; }

    //! Get the underlying Tpetra map
    const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > & getTpetra_Map() const { return map_; }

    //@}

  protected:

    RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > map_;

  }; // TpetraMap class

  // TODO: move that elsewhere
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & toTpetra(const Map<LocalOrdinal,GlobalOrdinal,Node> &map) {
    // TODO: throw exception
    const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> & tpetraMap = dynamic_cast<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &>(map);
    return *tpetraMap.getTpetra_Map();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > toTpetra(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map) {
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    if (map != Teuchos::null) {
      XPETRA_RCP_DYNAMIC_CAST(const TpetraMapClass, map, tpetraMap, "toTpetra");
      return tpetraMap->getTpetra_Map();
    }
    return Teuchos::null;
  }

  // In some cases (for instance, in MueLu adapter to Tpetra operator), we need to return a reference. This is only possible if
  // we assume that the map argument is nonzero
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >& toTpetraNonZero(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map) {
    TEUCHOS_TEST_FOR_EXCEPTION(map.is_null(), std::invalid_argument, "map must be nonzero");
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    XPETRA_RCP_DYNAMIC_CAST(const TpetraMapClass, map, tpetraMap, "toTpetra");
    return tpetraMap->getTpetra_Map();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > toXpetra(const RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > &map) {
    if (map != Teuchos::null)
      return rcp( new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(map));
    return Teuchos::null;
  }
  //

  // TODO: removed (but currently used in unit test)
  namespace useTpetra {

    //! Non-member function to create a locally replicated Map with a specified node.
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >
    createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
      XPETRA_MONITOR("useTpetra::createLocalMapWithNode");

      return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
    }

    //! Non-member function to create a (potentially) non-uniform, contiguous Map with the default node.
    template <class LocalOrdinal, class GlobalOrdinal>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> >
    createContigMap(global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
      XPETRA_MONITOR("useTpetra::createContigMap");

      return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createContigMap<LocalOrdinal,GlobalOrdinal>(numElements, localNumElements, comm)));
    }

    //! Non-member function to create a (potentially) non-uniform, contiguous Map with a user-specified node.
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> >
    createContigMapWithNode(global_size_t numElements, size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
      XPETRA_MONITOR("useTpetra::createContigMap");
      return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
    }
  } // useTpetra namespace

} // Xpetra namespace

// TODO: remove?
//!  Returns true if \c map is identical to this map. Implemented in TpetraMap::isSameAs().
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator== (const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map1, const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map2) {
  XPETRA_MONITOR("TpetraMap==TpetraMap");
  return map1.isSameAs(map2);
}

//! Returns true if \c map is not identical to this map. Implemented in TpetraMap::isSameAs().
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!= (const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map1, const Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node> &map2) {
  XPETRA_MONITOR("TpetraMap!=TpetraMap");
  return !map1.isSameAs(map2);
}

#define XPETRA_TPETRAMAP_SHORT
#endif // XPETRA_TPETRAMAP_HPP
