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

#ifndef TPETRA_DIRECTORY_HPP
#define TPETRA_DIRECTORY_HPP

#include <Teuchos_as.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_DirectoryImpl.hpp>
#include <Tpetra_DirectoryImpl_def.hpp>

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_Directory_decl.hpp"
#endif

namespace Tpetra {

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::Directory (const Map<LO, GO, NT>& map) :
    impl_ (NULL)
  {
    // Create an implementation object of the appropriate type,
    // depending on whether the Map is distributed or replicated, and
    // contiguous or noncontiguous.
    const Details::Directory<LO, GO, NT>* dir = NULL;
    if (map.isDistributed ()) {
      if (map.isUniform ()) {
        dir = new Details::ContiguousUniformDirectory<LO, GO, NT> (map);
      }
      else if (map.isContiguous ()) {
        dir = new Details::DistributedContiguousDirectory<LO, GO, NT> (map);
      }
      else {
        dir = new Details::DistributedNoncontiguousDirectory<LO, GO, NT> (map);
      }
    }
    else {
      dir = new Details::ReplicatedDirectory<LO, GO, NT> (map);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(dir == NULL, std::logic_error, "Tpetra::"
      "Directory constructor failed to create Directory implementation.  "
      "Please report this bug to the Tpetra developers.");
    impl_ = dir;
  }

  template<class LO, class GO, class NT>
  TEUCHOS_DEPRECATED
  Directory<LO, GO, NT>::Directory (const Teuchos::RCP<const Map<LO, GO, NT> >& mapPtr) :
    impl_ (NULL)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      mapPtr.is_null (), std::invalid_argument,
      "Tpetra::Directory(const RCP<const Map>& mapPtr): "
      "Input Map pointer is null.");
    const Map<LO, GO, NT>& map = *mapPtr;

    // Create an implementation object of the appropriate type,
    // depending on whether the Map is distributed or replicated, and
    // contiguous or noncontiguous.
    const Details::Directory<LO, GO, NT>* dir = NULL;
    if (map.isDistributed ()) {
      if (map.isUniform ()) {
        dir = new Details::ContiguousUniformDirectory<LO, GO, NT> (map);
      }
      else if (map.isContiguous ()) {
        dir = new Details::DistributedContiguousDirectory<LO, GO, NT> (map);
      }
      else {
        dir = new Details::DistributedNoncontiguousDirectory<LO, GO, NT> (map);
      }
    }
    else {
      dir = new Details::ReplicatedDirectory<LO, GO, NT> (map);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(dir == NULL, std::logic_error, "Tpetra::"
      "Directory constructor failed to create Directory implementation.  "
      "Please report this bug to the Tpetra developers.");
    impl_ = dir;
  }

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::
  Directory (const Map<LO, GO, NT>& map,
             const Tpetra::Details::TieBreak<LO,GO>& tieBreak) :
    impl_ (NULL)
  {
    // Create an implementation object of the appropriate type,
    // depending on whether the Map is distributed or replicated, and
    // contiguous or noncontiguous.
    //
    // mfh 06 Apr 2014: When a distributed noncontiguous Directory
    // takes a TieBreak, all the entries (local indices and process
    // ranks) owned by the Directory on the calling process pass
    // through the TieBreak object.  This may have side effects, such
    // as the TieBreak object remembering whether there were any
    // duplicates on the calling process.  We want to extend use of a
    // TieBreak object to other kinds of Directories.  For a
    // distributed contiguous Directory, the calling process owns all
    // of the (PID,LID) pairs in the input Map.  For a locally
    // replicated contiguous Directory, Process 0 owns all of the
    // (PID,LID) pairs in the input Map.
    //
    // It may seem silly to pass in a TieBreak when there are no ties
    // to break.  However, the TieBreak object gets to see all
    // (PID,LID) pairs that the Directory owns on the calling process,
    // and interface of TieBreak allows side effects.  Users may wish
    // to exploit them regardless of the kind of Map they pass in.
    const Details::Directory<LO, GO, NT>* dir = NULL;
    bool usedTieBreak = false;
    if (map.isDistributed ()) {
      if (map.isUniform ()) {
        dir = new Details::ContiguousUniformDirectory<LO, GO, NT> (map);
      }
      else if (map.isContiguous ()) {
        dir = new Details::DistributedContiguousDirectory<LO, GO, NT> (map);
      }
      else {
        dir = new Details::DistributedNoncontiguousDirectory<LO, GO, NT> (map, tieBreak);
        usedTieBreak = true;
      }
    }
    else {
      dir = new Details::ReplicatedDirectory<LO, GO, NT> (map);

      if (tieBreak.mayHaveSideEffects () && map.getNodeNumElements () != 0) {
        // We need the second clause in the above test because Map's
        // interface provides an inclusive range of local indices.
        const int myRank = map.getComm ()->getRank ();
        // In a replicated Directory, Process 0 owns all the
        // Directory's entries.  This is an arbitrary assignment; any
        // one process would do.
        if (myRank == 0) {
          std::vector<std::pair<int, LO> > pidLidList (1);
          const LO minLocInd = map.getMinLocalIndex ();
          const LO maxLocInd = map.getMaxLocalIndex ();
          for (LO locInd = minLocInd; locInd <= maxLocInd; ++locInd) {
            pidLidList[0] = std::make_pair (myRank, locInd);
            const GO globInd = map.getGlobalElement (locInd);
            // We don't care about the return value; we just want to
            // invoke the side effects.
            (void) tieBreak.selectedIndex (globInd, pidLidList);
          }
        }
      }
      usedTieBreak = true;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(dir == NULL, std::logic_error, "Tpetra::"
      "Directory constructor failed to create Directory implementation.  "
      "Please report this bug to the Tpetra developers.");

    if (! usedTieBreak && tieBreak.mayHaveSideEffects () &&
        map.getNodeNumElements () != 0) {
      // We need the third clause in the above test because Map's
      // interface provides an inclusive range of local indices.
      std::vector<std::pair<int, LO> > pidLidList (1);
      const LO minLocInd = map.getMinLocalIndex ();
      const LO maxLocInd = map.getMaxLocalIndex ();
      const int myRank = map.getComm ()->getRank ();
      for (LO locInd = minLocInd; locInd <= maxLocInd; ++locInd) {
        pidLidList[0] = std::make_pair (myRank, locInd);
        const GO globInd = map.getGlobalElement (locInd);
        // We don't care about the return value; we just want to
        // invoke the side effects.
        (void) tieBreak.selectedIndex (globInd, pidLidList);
      }
    }

    impl_ = dir;
  }

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::Directory () : impl_ (NULL) {}

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::~Directory () {
    if (impl_ != NULL) {
      delete impl_;
      impl_ = NULL;
    }
  }

  template<class LO, class GO, class NT>
  LookupStatus
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Map<LO, GO, NT>& map,
                       const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs) const
  {
    const bool computeLIDs = false;
    return impl_->getEntries (map, globalIDs, nodeIDs, Teuchos::null, computeLIDs);
  }


  template<class LO, class GO, class NT>
  LookupStatus TEUCHOS_DEPRECATED
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "This method is DEPRECATED.  Please call the "
      "three-argument version of getDirectoryEntries that takes a Map, an "
      "ArrayView<const GO>, and an ArrayView<int>.");
  }

  template<class LO, class GO, class NT>
  LookupStatus
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Map<LO, GO, NT>& map,
                       const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs,
                       const Teuchos::ArrayView<LO>& localIDs) const
  {
    const bool computeLIDs = true;
    return impl_->getEntries (map, globalIDs, nodeIDs, localIDs, computeLIDs);
  }

  template<class LO, class GO, class NT>
  LookupStatus TEUCHOS_DEPRECATED
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs,
                       const Teuchos::ArrayView<LO>& localIDs) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "This method is DEPRECATED.  Please call the "
      "four-argument version of getDirectoryEntries that takes a Map, an "
      "ArrayView<const GO>, an ArrayView<int>, and an ArrayView<LO>.");
  }

  template<class LO, class GO, class NT>
  bool Directory<LO, GO, NT>::isOneToOne (const Teuchos::Comm<int>& comm) const {
    return impl_->isOneToOne (comm);
  }

  template<class LO, class GO, class NT>
  std::string
  Directory<LO, GO, NT>::description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "Directory"
       << "<" << TypeNameTraits<LO>::name ()
       << ", " << TypeNameTraits<GO>::name ()
       << ", " << TypeNameTraits<NT>::name () << ">";
    return os.str ();
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_DIRECTORY_INSTANT(LO,GO,NODE) \
  \
  template class Directory< LO , GO , NODE >; \

#endif // TPETRA_DIRECTORY_HPP
