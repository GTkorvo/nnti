// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_DISTOBJECT_HPP
#define TPETRA_DISTOBJECT_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Distributor.hpp"

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

namespace Tpetra {

  /// \class DistObject
  /// \brief A base class for distributed objects that support import and export operations.
  ///
  /// DistObject is a base class for all Tpetra distributed global
  /// objects.  It provides the basic mechanisms and interface
  /// specifications for importing and exporting operations using \c
  /// Tpetra::Import and \c Tpetra::Export objects.
  ///
  /// <b> Distributed Global vs. Replicated Local.</b>
  ///
  /// <ul>
  /// <li> Distributed Global objects - In most instances, a
  /// distributed object will be partitioned across multiple memory
  /// images associated with multiple processors.  In this case, there
  /// is a unique copy of each element and elements are spread across
  /// all images specified by the Teuchos::Comm communicator object.
  /// <li> Replicated Local Objects - Some algorithms use objects that
  /// are too small to be distributed across all processors.  The
  /// upper Hessenberg matrix in a GMRES iterative solve is a good
  /// example.  computation.  In other cases, such as with block
  /// iterative methods, block dot product functions produce small
  /// dense matrices that are required by all images.  Replicated
  /// local objects handle these situations.
  /// </ul>
  template <class Packet, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class DistObject : virtual public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! constructor
    explicit DistObject(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);

    //! copy constructor
    DistObject(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source);

    //! destructor
    virtual ~DistObject();

    //@}

    //! @name Import/Export Methods
    //@{ 

    //! Import
    void 
    doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source, 
	      const Import<LocalOrdinal,GlobalOrdinal,Node> &importer, 
	      CombineMode CM);

    //! Export
    void doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &dest, 
		   const Export<LocalOrdinal,GlobalOrdinal,Node> &exporter, 
		   CombineMode CM);

    //! Import (using an Exporter)
    void doImport(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source,
                  const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter, CombineMode CM);

    //! Export (using an Importer)
    void doExport(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &dest,
                  const Import<LocalOrdinal,GlobalOrdinal,Node>& importer, CombineMode CM);

    //@}

    //! @name Attribute Accessor Methods
    //@{ 

    /// \brief True if this is a globally distributed object, else false.
    ///
    /// For a definition of "globally distributed" (and its opposite,
    /// "locally replicated"), see the documentation of \c Map.
    inline bool isDistributed() const;

    //! The Map with which this DistObject was constructed.
    inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& 
    getMap() const { return map_; }

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method.

    void print(std::ostream &os) const;

    //@} 

  protected:

    /// \enum ReverseOption
    /// \brief Whether the data transfer should be performed in forward or reverse mode.
    enum ReverseOption {
      DoForward, //*!< Perform the transfer in forward mode.
      DoReverse  //*!< Perform the transfer in reverse mode.
    };

    //! Perform transfer (redistribution) of data across memory images.
    virtual void doTransfer(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> &source,
                            CombineMode CM,
                            size_t numSameIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &remoteLIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                            Distributor &distor,
                            ReverseOption revOp);

    // The following four methods must be implemented by the derived class

    /// \brief Compare the source and target (\e this) objects for compatibility.
    ///
    /// \return True if they are compatible, else false.
    virtual bool 
    checkSizes (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>& source) = 0;

    //! Perform copies and permutations that are local to this image.
    /*!
      \param source [in]
             On entry, the DistObject that we are importing from.
      \param numSameIDs [in]
             On entry, the number of elements that are the same on the source and dest objects.
         (i.e. The element is owned by the same image in both source and dest, 
         and no permutation occurs.)
      \param numPermuteIDs [in]
             On entry, the number of elements that are locally permuted between source and dest objects.
      \param permuteToLIDs [in]
             On entry, contains a list of the elements that are permuted. (Listed by their LID in the
         destination DistObject.)
      \param permuteFromLIDs [in]
             On entry, contains a list of the elements that are permuted. (Listed by their LID in the
         source DistObject.)
    */
    virtual void copyAndPermute(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & source,
                                size_t numSameIDs,
                                const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                                const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs) = 0;

    //! Perform any packing or preparation required for communication.
    /*!
      \param source [in]
             On entry, the DistObject that we are importing from.
      \param exportLIDs [in]
             On entry, a list of the entries we will be sending to other images.
             (Listed by their LID in the source DistObject.)
      \param exports [out]
             On exit, buffer for data we will be sending out.
      \param numPacketsPerLID [out]
             On exit, numPacketsPerLID[i] contains the number of packets to be exported for
             exportLIDs[i].
      \param constantNumPackets [out]
             On exit, 0 if numPacketsPerLID has variable contents (different size for each LID).
             If nonzero, then it is expected that num-packets-per-LID is constant, and
             constantNumPackets holds that value.
      \param distor [in]
             On entry, contains the Distributor object we are using.         
    */
    virtual void packAndPrepare(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & source,
                                const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                                Teuchos::Array<Packet> &exports,
                                const Teuchos::ArrayView<size_t> & numPacketsPerLID,
                                size_t& constantNumPackets,
                                Distributor &distor) = 0;

    //! Perform any unpacking and combining after communication.
    /*!
      \param importLIDs [in]
             On entry, a list of the entries we received from other images.
             (Listed by their LID in the target DistObject.)
      \param imports [in]
             Buffer containing data we received.
      \param numPacketsPerLID [in]
             numPacketsPerLID[i] contains the number of packets imported for
             importLIDs[i].
      \param constantNumPackets [in]
             If nonzero, then numPacketsPerLID is constant (same value in all entries)
             and constantNumPackets is that value.
      \param distor [in]
             The Distributor object we are using.
      \param CM [in]
             The Tpetra::CombineMode to use when combining the imported entries with existing entries.
    */
    virtual void 
    unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
		      const Teuchos::ArrayView<const Packet> &imports,
		      const Teuchos::ArrayView<size_t> &numPacketsPerLID,
		      size_t constantNumPackets,
		      Distributor &distor,
		      CombineMode CM) = 0;

    virtual void createViews() const {}
    virtual void createViewsNonConst(Kokkos::ReadWriteOption rwo) {}
    virtual void releaseViews() const {}

    //! The Map over which this object is distributed.
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > map_;

  private:
    // buffers into which packed data is imported
    Teuchos::Array<Packet> imports_;
    Teuchos::Array<size_t> numImportPacketsPerLID_;
    // buffers from which packed data is exported
    Teuchos::Array<Packet> exports_;
    Teuchos::Array<size_t> numExportPacketsPerLID_;

  }; // class DistObject

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::DistObject(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & map)
  : map_(map)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::DistObject(const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & source)
  : map_(source.map_)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::~DistObject() 
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A, 
	    const Import<LocalOrdinal,GlobalOrdinal,Node> & importer, 
	    CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(   *getMap() != *importer.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *importer.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    size_t numSameIDs = importer.getNumSameIDs();
    const Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = importer.getExportLIDs();
    const Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = importer.getRemoteLIDs();
    const Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = importer.getPermuteToLIDs();
    const Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
                     importer.getDistributor(), DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A, 
	    const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter, 
	    CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(   *getMap() != *exporter.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *exporter.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    size_t numSameIDs = exporter.getNumSameIDs();
    Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = exporter.getExportLIDs();
    Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = exporter.getRemoteLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = exporter.getPermuteToLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               exporter.getDistributor(), DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doImport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
	    const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter, 
	    CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(  * getMap() != *exporter.getSourceMap(), std::runtime_error, "Target Maps don't match.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *exporter.getTargetMap(), std::runtime_error, "Source Maps don't match.");
    size_t numSameIDs = exporter.getNumSameIDs();
    Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = exporter.getRemoteLIDs();
    Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = exporter.getExportLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = exporter.getPermuteFromLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               exporter.getDistributor(), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::
  doExport (const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & A,
	    const Import<LocalOrdinal,GlobalOrdinal,Node> & importer, 
	    CombineMode CM) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION( *  getMap() != *importer.getSourceMap(), std::runtime_error, "Target Maps don't match.");
    TEUCHOS_TEST_FOR_EXCEPTION( *A.getMap() != *importer.getTargetMap(), std::runtime_error, "Source Maps don't match.");
    size_t numSameIDs = importer.getNumSameIDs();
    Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = importer.getRemoteLIDs();
    Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = importer.getExportLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = importer.getPermuteFromLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               importer.getDistributor(), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::isDistributed() const {
    return map_->isDistributed();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::doTransfer(
      const DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node> & source,
      CombineMode CM,
      size_t numSameIDs, 
      const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs, 
      const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs,
      const Teuchos::ArrayView<const LocalOrdinal> &remoteLIDs,    
      const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
      Distributor &distor, ReverseOption revOp) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION( checkSizes(source) == false, std::runtime_error, 
        "Tpetra::DistObject::doTransfer(): checkSizes() indicates that DistOjbects are not size-compatible.");
    Kokkos::ReadWriteOption rwo = Kokkos::ReadWrite;
    if (CM == INSERT || CM == REPLACE) {
      if (numSameIDs + permuteToLIDs.size() + remoteLIDs.size() == this->getMap()->getNodeNumElements()) {
        // we are overwriting, and everything will be overwritten
        rwo = Kokkos::WriteOnly;
      }
    }
    source.createViews();
    this->createViewsNonConst(rwo); 
    if (numSameIDs + permuteToLIDs.size()) {
      copyAndPermute(source,numSameIDs,permuteToLIDs,permuteFromLIDs);
    }
    size_t constantNumPackets = 0;
    numExportPacketsPerLID_.resize(exportLIDs.size());
    numImportPacketsPerLID_.resize(remoteLIDs.size());
    packAndPrepare(source,exportLIDs,exports_,numExportPacketsPerLID_(),constantNumPackets,distor);
    source.releaseViews();
    if (constantNumPackets != 0) {
      size_t rbufLen = remoteLIDs.size()*constantNumPackets;
      imports_.resize(rbufLen);
    }
    if ((isDistributed() && revOp == DoReverse) || (source.isDistributed() && revOp == DoForward)) 
    {
      // call one of the doPostsAndWaits functions
      if (revOp == DoReverse) {
        if (constantNumPackets == 0) { //variable num-packets-per-LID:
          distor.doReversePostsAndWaits(numExportPacketsPerLID_().getConst(), 1,
                                        numImportPacketsPerLID_());
          size_t totalImportPackets = 0;
          for(Array_size_type i=0; i<numImportPacketsPerLID_.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID_[i];
          }
          imports_.resize(totalImportPackets);
          distor.doReversePostsAndWaits(exports_().getConst(),numExportPacketsPerLID_(),
                                        imports_(), numImportPacketsPerLID_());
        }
        else {
          distor.doReversePostsAndWaits(exports_().getConst(),constantNumPackets,imports_());
        }
      }
      else {
        if (constantNumPackets == 0) { //variable num-packets-per-LID:
          distor.doPostsAndWaits(numExportPacketsPerLID_().getConst(), 1,
                                 numImportPacketsPerLID_());
          size_t totalImportPackets = 0;
          for(Array_size_type i=0; i<numImportPacketsPerLID_.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID_[i];
          }
          imports_.resize(totalImportPackets);
          distor.doPostsAndWaits(exports_().getConst(),numExportPacketsPerLID_(),
                                 imports_(), numImportPacketsPerLID_());
        }
        else {
          distor.doPostsAndWaits(exports_().getConst(),constantNumPackets,imports_());
        }
      }
      unpackAndCombine(remoteLIDs,imports_(),numImportPacketsPerLID_(), constantNumPackets, distor,CM);
    }
    this->releaseViews();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void DistObject<Packet,LocalOrdinal,GlobalOrdinal,Node>::print(std::ostream &os) const
  {
    using std::endl;
    os << "Tpetra::DistObject" << endl
       << " export buffer size: " << exports_.size() << endl
       << " import buffer size: " << imports_.size() << endl
       << "Map:" << endl
       << map_;
  }

} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_HPP */
