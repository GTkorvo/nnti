// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#ifndef MUELU_REPARTITIONFACTORY_DEF_HPP
#define MUELU_REPARTITIONFACTORY_DEF_HPP

#include <iostream>

#include "MueLu_RepartitionFactory_decl.hpp" // TMP JG NOTE: before other includes, otherwise I cannot test the fwd declaration in _def

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Hashtable.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <MueLu_UCAggregationCommHelper.hpp>

#include "MueLu_Utilities.hpp" // TMP JG NOTE: only for maxAll, so no _fwd in _decl

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RepartitionFactory(
                LO minRowsPerProcessor, double nnzMaxMinRatio, GO startLevel, LO useDiffusiveHeuristic, GO minNnzPerProcessor) :
    minRowsPerProcessor_(minRowsPerProcessor),
    nnzMaxMinRatio_(nnzMaxMinRatio),
    startLevel_(startLevel),
    useDiffusiveHeuristic_(useDiffusiveHeuristic),
    minNnzPerProcessor_(minNnzPerProcessor)
  { }

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~RepartitionFactory() {}

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "Partition");
  }

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    using Teuchos::Array;
    using Teuchos::ArrayRCP;

    //TODO: We only need a CrsGraph. This class does not have to be templated on Scalar types.
    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

    { // scoping

      // ======================================================================================================
      // Determine whether partitioning is needed.
      // ======================================================================================================
      // TODO: Result of the test can be stored in the level to avoid doing this several time.

      // Repartitionning iff !(Test1 || Test2) && (Test3 || Test4)
      // This is implemented with a short-circuit evaluation

      // Test1: skip partitioning if level is too big
      std::ostringstream msg1; msg1 << std::endl << "    current level = " << currentLevel.GetLevelID() << ", first level where repartitioning can happen is " << startLevel_ << ".";
      if (currentLevel.GetLevelID() < startLevel_) {
        throw(MueLu::Exceptions::HaltRepartitioning("No repartitioning necessary:" + msg1.str()));
      }

      // Test 2: check whether A is spread over more than one process.
      // Note: using type 'int' as it is enough (comm->getSize() is an int). Test can also be done with 'bool' but numActiveProcesses is printed later
      // TODO: this global communication can be avoided if we store the information with the matrix (it is known when matrix is created)
      int numActiveProcesses = 0;
      sumAll(comm, (int)((A->getNodeNumRows() > 0) ? 1 : 0), numActiveProcesses);
      std::ostringstream msg2; msg2 << std::endl << "    # processes with rows = " << numActiveProcesses;
      if (numActiveProcesses == 1) {
        throw(MueLu::Exceptions::HaltRepartitioning("No repartitioning necessary:" + msg2.str()));
      }

      bool doRepartition = false;
      std::ostringstream msg3, msg4, *msgDoRepartition=NULL; // msgDoRepartition == msg3 or msg4 depending on Test3 and Test4 evalution

      // Test3: check whether any node has too few rows
      // Note: (!Test2) ensures that repartitionning is not done when only 1 proc and globalNumRow < minRowsPerProcessor_
      if (minRowsPerProcessor_ > 0) { // skip test if criteria not used
        LO     minNumRows;
        size_t numMyRows = A->getNodeNumRows();
        LO LOMAX = Teuchos::OrdinalTraits<LO>::max(); // processors without rows do not participate to minAll()
        minAll(comm, (LO)((numMyRows > 0) ? numMyRows : LOMAX), minNumRows);
        TEUCHOS_TEST_FOR_EXCEPTION(minNumRows >= LOMAX, Exceptions::RuntimeError, "internal error");
        msg3 << std::endl << "    min # rows per proc = "   << minNumRows << ", min allowable = " << minRowsPerProcessor_;
        if (minNumRows < minRowsPerProcessor_) {
          doRepartition = true;
          msgDoRepartition = &msg3;
        }
      }

      // Test4: check whether the number of nonzeros per process is imbalanced
      if (doRepartition == false) { // skip test if we already now that we will do repartitioning
        size_t numMyNnz  = A->getNodeNumEntries();
        GO minNnz, maxNnz;
        maxAll(comm, (GO)numMyNnz, maxNnz);
        minAll(comm, (GO)((numMyNnz > 0) ? numMyNnz : maxNnz), minNnz); // min nnz over all proc (disallow any processors with 0 nnz)
        double imbalance = ((double) maxNnz) / minNnz;
        msg4 << std::endl << "    nonzero imbalance = " << imbalance  << ", max allowable = " << nnzMaxMinRatio_;
        if (imbalance > nnzMaxMinRatio_) {
          doRepartition = true;
          msgDoRepartition = &msg4;
        }
      }

      if (!doRepartition) {
        throw(MueLu::Exceptions::HaltRepartitioning("No repartitioning necessary:" + msg1.str() + msg2.str() + msg3.str() + msg4.str()));
      }

      // print only conditions that triggered the repartitioning
      GetOStream(Statistics0, 0) << "Repartitioning necessary:" << msg1.str() << msg2.str() << msgDoRepartition->str() << std::endl;

    } // scoping

    // FIXME Quick way to figure out how many partitions there should be. (Same as what's done in ML.)
    // FIXME Should take into account nnz, maybe?  Perhaps only when user is using min #nnz per row threshold.
    GO numPartitions;
    if (currentLevel.IsAvailable("number of partitions")) {
      numPartitions = currentLevel.Get<GO>("number of partitions");
    } else {
      if ((GO)A->getGlobalNumRows() < minRowsPerProcessor_) numPartitions = 1;
      else                                                  numPartitions = A->getGlobalNumRows() / minRowsPerProcessor_;
      if (numPartitions > comm->getSize())
        numPartitions = comm->getSize();
      currentLevel.Set("number of partitions", numPartitions);
    }
    GetOStream(Statistics0, 0) << "Number of partitions to use = " << numPartitions << std::endl;

    // ======================================================================================================
    // Determine the global size of each partition.
    // ======================================================================================================
    // Length of vector "decomposition" is local number of DOFs.  Its entries are partition numbers each DOF belongs to.

    RCP<GOVector> decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");

    // Use a hashtable to record how many local rows belong to each partition.
    RCP<Teuchos::Hashtable<GO, GO> > hashTable;
    hashTable = rcp(new Teuchos::Hashtable<GO, GO>(numPartitions + numPartitions/2));
    ArrayRCP<const GO> decompEntries;
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    bool flag=false;
    for (int i=0; i<decompEntries.size(); ++i) {
      if (decompEntries[i] >= numPartitions) flag = true;
      if (hashTable->containsKey(decompEntries[i])) {
        GO count = hashTable->get(decompEntries[i]);
        ++count;
        hashTable->put(decompEntries[i], count);
      } else {
        hashTable->put(decompEntries[i], 1);
      }
    }
    int problemPid;
    int mypid = comm->getRank();
    maxAll(comm, (flag ? mypid : -1), problemPid);
    std::ostringstream buf; buf << problemPid;
    TEUCHOS_TEST_FOR_EXCEPTION(problemPid>-1, Exceptions::RuntimeError, "pid " + buf.str() + " encountered a partition number is that out-of-range");
    decompEntries = Teuchos::null;

    Teuchos::Array<GO> allPartitionsIContributeTo;
    Teuchos::Array<GO> allLocalPartSize;
    hashTable->arrayify(allPartitionsIContributeTo, allLocalPartSize);

    GO indexBase = decomposition->getMap()->getIndexBase();

    // Source map is overlapping.  GIDs owned by this pid are the partition numbers found above.
    RCP<Map> sourceMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                           allPartitionsIContributeTo(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

    // Store # of local DOFs in each partition in a vector based on above map.
    RCP<GOVector> localPartSizeVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(sourceMap, false);
    ArrayRCP<GO> data;
    if (localPartSizeVec->getLocalLength() > 0)
      data = localPartSizeVec->getDataNonConst(0);
    for (int i=0; i<hashTable->size(); ++i)
      data[i] = allLocalPartSize[i];
    data = Teuchos::null;

    // Target map is nonoverlapping.  Pid k has GID N if and only if k owns partition N.
    GO myPartitionNumber;
    Array<int> partitionOwners;
    DeterminePartitionPlacement(currentLevel, myPartitionNumber, partitionOwners);

/*
    // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
    std::cout << "pid " << mypid << " here" << std::endl;
    sleep(1); comm->barrier();
    for (int i=0; i<comm->getSize(); ++i) {
      if (mypid == i) {
        std::cout << "pid " << mypid << " owns partition " << myPartitionNumber << std::endl;
        for (int j=0; j<partitionOwners.size(); ++j)
          std::cout << "     partition " << j << " owned by pid " << partitionOwners[j] << std::endl;
      }
      comm->barrier();
    }
    sleep(1); comm->barrier();
    // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
*/

    GO numDofsThatStayWithMe=0;
    Teuchos::Array<GO> partitionsIContributeTo;
    Teuchos::Array<GO> localPartSize;
    for (int i=0; i<allPartitionsIContributeTo.size(); ++i)
    {
      if (allPartitionsIContributeTo[i] != myPartitionNumber) {
        partitionsIContributeTo.push_back(allPartitionsIContributeTo[i]);
        localPartSize.push_back(allLocalPartSize[i]);
      }
      else
        numDofsThatStayWithMe = allLocalPartSize[i];
    }

    // Note: "numPartitionsISendTo" does not include this PID
    GO numPartitionsISendTo = hashTable->size();
    // I'm a partition owner, so don't count my own.
    if (myPartitionNumber >= 0 && numPartitionsISendTo>0 && numDofsThatStayWithMe>0) numPartitionsISendTo--;
    assert(numPartitionsISendTo == partitionsIContributeTo.size());

    Array<int> partitionOwnersISendTo;
    for (int i=0; i<partitionsIContributeTo.size(); ++i) {
      partitionOwnersISendTo.push_back(partitionOwners[partitionsIContributeTo[i]]);
    }

    Array<GO> localMapElement;
    if (myPartitionNumber >= 0)
      localMapElement.push_back(myPartitionNumber);
    RCP<Map> targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           numPartitions,
                                           localMapElement(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

    RCP<const Export> exporter = ExportFactory::Build( sourceMap, targetMap);

    // If this pid owns a partition, globalPartSizeVec has one local entry that is the global size of said partition.
    // If this pid doesn't own a partition, globalPartSizeVec is locally empty.
    RCP<GOVector> globalPartSizeVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(targetMap);
    globalPartSizeVec->doExport(*localPartSizeVec, *exporter, Xpetra::ADD);
    int myPartitionSize = 0;
    ArrayRCP<const GO> constData;
    if (globalPartSizeVec->getLocalLength() > 0) {
      constData = globalPartSizeVec->getData(0);
      myPartitionSize = constData[0];
    }
    constData = Teuchos::null;

    // ======================================================================================================
    // Calculate how many PIDs (other than myself) contribute to my partition.
    // ======================================================================================================
    RCP<GOVector> howManyPidsSendToThisPartitionVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(targetMap);
    RCP<GOVector> partitionsISendTo = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(sourceMap, false);
    if (partitionsISendTo->getLocalLength() > 0)
      data = partitionsISendTo->getDataNonConst(0);
    for (int i=0; i<data.size(); ++i) {
      // don't count myself as someone I send to... (sourceMap is based on allPartitionsIContributeTo)
      if (sourceMap->getGlobalElement(i) != myPartitionNumber)
        data[i] = 1;
      else
        data[i] = 0;
    }
    data = Teuchos::null;

    // Note: "howManyPidsSendToThisPartition" does not include this PID
    int howManyPidsSendToThisPartition = 0;
    howManyPidsSendToThisPartitionVec->doExport(*partitionsISendTo, *exporter, Xpetra::ADD);
    if (howManyPidsSendToThisPartitionVec->getLocalLength() > 0) {
      constData = howManyPidsSendToThisPartitionVec->getDataNonConst(0);
      howManyPidsSendToThisPartition = constData[0];
    }
    constData = Teuchos::null;

    // ======================================================================================================
    // Calculate which PIDs contribute to my partition, and how much each contributes,
    // with ireceive/send/wait cycle.
    // FIXME Jan.12.2012 Teuchos::Comm methods ireceive and wait don't work (bugs 5483 and 5484), so for
    // now use the raw MPI methods.
    // ======================================================================================================

    Array<GO> pidsIReceiveFrom(howManyPidsSendToThisPartition);
    Array<GO> numDofsIReceiveFromOnePid(howManyPidsSendToThisPartition);
    for (int j=0; j<numDofsIReceiveFromOnePid.size(); ++j)
      numDofsIReceiveFromOnePid[j] = -99;

    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    if (tmpic == Teuchos::null)
      throw(Exceptions::RuntimeError("Cannot cast base Teuchos::Comm to Teuchos::MpiComm object."));
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();

    // First post non-blocking receives.
    Array<MPI_Request> requests(howManyPidsSendToThisPartition);
    for (int i=0; i<howManyPidsSendToThisPartition; ++i) {
      MPI_Irecv((void*)&(numDofsIReceiveFromOnePid[i]), 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, *rawMpiComm, &requests[i]);
    }

    // Next post sends.
    for (int i=0; i< partitionsIContributeTo.size(); ++i) {
      comm->send(sizeof(GO), (char*)&localPartSize[i], partitionOwnersISendTo[i]);
    }

    // Finally do waits.
    Array<MPI_Status> status(howManyPidsSendToThisPartition);
    for (int i=0; i<howManyPidsSendToThisPartition; ++i)
      MPI_Wait(&requests[i], &status[i]);

    for (int i=0; i<pidsIReceiveFrom.size(); ++i)
      pidsIReceiveFrom[i] = status[i].MPI_SOURCE;

    comm->barrier();

    // =================================================================================================
    // Calculate partial offsets for permutation row map, via MPI_Scan based on global partition sizes.
    // Communicate those offsets back to respective PIDS of unpermuted matrix using ireceive/send/wait cycle.
    // =================================================================================================
    // Partition numbers are used as tags to ensure messages go to correct PIDs.
    // Use MPI_DOUBLE type to avoid any overflow problems.
    // partitionSizeOffset is the first GID in this partition.
    // gidOffsets is an array of GID offsets.  gidOffsets[i] is the first GID for the dofs received from PID partitionOwnersISendTo[i].
    // Note: Quantities "numPartitionsISendTo" and "howManyPidsSendToThisPartition" do not include me
    double partitionSizeOffset;
    double ttt = myPartitionSize;
    MPI_Scan(&ttt, &partitionSizeOffset, 1, MPI_DOUBLE, MPI_SUM, *rawMpiComm);
    partitionSizeOffset -= myPartitionSize;

    Array<double> gidOffsets;
    if (howManyPidsSendToThisPartition > 0) {
      gidOffsets.resize(howManyPidsSendToThisPartition);
      gidOffsets[0] = partitionSizeOffset;
    }
    for (int i=1; i<howManyPidsSendToThisPartition; ++i) {
      gidOffsets[i] = gidOffsets[i-1] + numDofsIReceiveFromOnePid[i-1];
    }

    requests.resize(numPartitionsISendTo);
    Array<double> gidOffsetsForPartitionsIContributeTo(numPartitionsISendTo);

    // Post receives on contributing PIDs.
    for (int i=0; i< numPartitionsISendTo; ++i) {
      int msgTag = partitionsIContributeTo[i];
      MPI_Irecv((void*)&(gidOffsetsForPartitionsIContributeTo[i]), 1, MPI_DOUBLE, partitionOwnersISendTo[i], msgTag, *rawMpiComm, &requests[i]);
    }

    // Do sends by partition owners.
    for (int i=0; i<howManyPidsSendToThisPartition; ++i) {
      int msgTag = myPartitionNumber;
      MPI_Send((void*)&gidOffsets[i] , 1 , MPI_DOUBLE , pidsIReceiveFrom[i], msgTag, *rawMpiComm);
    }

    // Do waits.
    status.resize(numPartitionsISendTo);
    for (int i=0; i<numPartitionsISendTo; ++i)
      MPI_Wait(&requests[i], &status[i]);

    comm->barrier();

    // =================================================================================================
    // Set up a synthetic GID scheme that is the same for both the original unpermuted system and the permuted system.
    // This scheme is *not* the final DOF numbering, but is just for the importer we need to transfer column IDs.
    // =================================================================================================

    // Synthetic GIDS for original unpermuted system.

    // store offsets for easy random access
    hashTable = rcp(new Teuchos::Hashtable<GO, GO>(partitionsIContributeTo.size() + partitionsIContributeTo.size()/2));
    for (int i=0; i<partitionsIContributeTo.size(); ++i) {
      hashTable->put(partitionsIContributeTo[i], (GO)gidOffsetsForPartitionsIContributeTo[i]);
    }
    //store gid offset for those dofs that will remain with me
    if (myPartitionNumber > -1)
      hashTable->put(myPartitionNumber, ((GO)partitionSizeOffset) + myPartitionSize - numDofsThatStayWithMe);
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    Array<GO> uniqueGIDsBeforePermute;
    for (int i=0; i<decompEntries.size(); ++i) {
      GO gid = hashTable->get(decompEntries[i]);
      uniqueGIDsBeforePermute.push_back(gid);
      gid++;
      hashTable->put(decompEntries[i], gid);
    }
    decompEntries = Teuchos::null;

    // Synthetic GIDS for permuted system.

    Array<GO> uniqueGIDsAfterPermute;
    for (int i=0; i<myPartitionSize; ++i) {
        uniqueGIDsAfterPermute.push_back((GO)partitionSizeOffset+i);
    }

    // =================================================================================================
    // Create and apply an importer to communicate column GIDs for the permutation matrix.
    // =================================================================================================

    //TODO we should really supply the global size as another sanity check
    sourceMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    uniqueGIDsBeforePermute(),
                    indexBase,
                    comm);

    targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    uniqueGIDsAfterPermute(),
                    indexBase,
                    comm);

    RCP<const Import> importer = ImportFactory::Build( sourceMap, targetMap);

    RCP<GOVector> sourceVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(sourceMap);
    ArrayRCP<GO> vectorData;
    if (sourceVec->getLocalLength() > 0) {
      vectorData = sourceVec->getDataNonConst(0);
    }

    //load source vector with unpermuted GIDs.
    RCP<const Map> originalRowMap = A->getRowMap();

    //sanity check
    assert(vectorData.size() == (GO) originalRowMap->getNodeNumElements());


    for (LO i=0; i<vectorData.size(); ++i) {
      GO gid = originalRowMap->getGlobalElement(i);
      if (gid == (GO) Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid())
        throw(Exceptions::RuntimeError("Encountered an unexpected error with A's rowmap."));
      vectorData[i] = gid;
    }
    vectorData = Teuchos::null;

    RCP<GOVector> targetVec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(targetMap);
    targetVec->doImport(*sourceVec, *importer, Xpetra::INSERT);

    // =================================================================================================
    // Create an importer between the original row map and the permuted row map
    // =================================================================================================
    vectorData = targetVec->getDataNonConst(0);
    RCP<Map> newRowMap = MapFactory::Build(decomposition->getMap()->lib(),
                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                    vectorData(),
                    indexBase,
                    comm);

    RCP<const Import> importerForRepartitioning = ImportFactory::Build(A->getRowMap(), newRowMap);
    Set(currentLevel, "Importer", importerForRepartitioning);

  } //Build

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeterminePartitionPlacement(Level & currentLevel, GO &myPartitionNumber,
  Array<int> &partitionOwners) const
  {
    FactoryMonitor m(*this, "DeterminePartitionPlacement", currentLevel);

    GO numPartitions = Get<GO>(currentLevel, "number of partitions");
    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
    int mypid = comm->getRank();

    /*
      useDiffusiveHeuristic_ = 0      ====> always put on procs 0..N
      useDiffusiveHeuristic_ = 1      ====> always use diffusive heuristic
      useDiffusiveHeuristic_ = K > 1  ====> if #partitions is > K, put on procs 0..N
                                            otherwise use diffusive
      useDiffusiveHeuristic_ = -1     ====> put on procs 0..N this time only, then use diffusive in remaining rounds
    */
    if (useDiffusiveHeuristic_ != 1) {

      if (numPartitions > useDiffusiveHeuristic_) {
        if (numPartitions==1)
          GetOStream(Runtime0, 0) << "Placing partitions on proc. 0." << std::endl;
        else
          GetOStream(Runtime0, 0) << "Placing partitions on proc. 0-" << numPartitions-1 << "." << std::endl;
        myPartitionNumber = -1;
        for (int i=0; i<(int)numPartitions; ++i) {
          partitionOwners.push_back(i);
          if (i==mypid) myPartitionNumber = i;
        }
        if (useDiffusiveHeuristic_ == -1)
          useDiffusiveHeuristic_ = 1;
        return;
      }
    }

    GetOStream(Runtime0, 0) << "Using diffusive heuristic for partition placement." << std::endl;

    RCP<SubFactoryMonitor> m1 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: Setup", currentLevel));

//RCP<SubFactoryMonitor> m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: getting 'Partition'", currentLevel));
    RCP<GOVector> decomposition = Get<RCP<GOVector> >(currentLevel, "Partition");
    // Figure out how many nnz there are per row.
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: figuring out nnz per row", currentLevel));
    RCP<GOVector> nnzPerRowVector = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO> nnzPerRow;
    if (nnzPerRowVector->getLocalLength() > 0)
      nnzPerRow = nnzPerRowVector->getDataNonConst(0);
    for (int i=0; i<nnzPerRow.size(); ++i)
      nnzPerRow[i] = A->getNumEntriesInLocalRow(i);

    // Use a hashtable to record how many nonzeros in the local matrix belong to each partition.
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: hashing", currentLevel));
    RCP<Teuchos::Hashtable<GO, GO> > hashTable;
    hashTable = rcp(new Teuchos::Hashtable<GO, GO>(numPartitions + numPartitions/2));
    ArrayRCP<const GO> decompEntries;
    if (decomposition->getLocalLength() > 0)
      decompEntries = decomposition->getData(0);
    bool flag=false;
    assert(decompEntries.size() == nnzPerRow.size());
    for (int i=0; i<decompEntries.size(); ++i) {
      if (decompEntries[i] >= numPartitions) flag = true;
      if (hashTable->containsKey(decompEntries[i])) {
        GO count = hashTable->get(decompEntries[i]);
        count += nnzPerRow[i];
        hashTable->put(decompEntries[i], count);
      } else {
        hashTable->put(decompEntries[i], nnzPerRow[i]);
      }
    }

//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: arrayify", currentLevel));
    int problemPid;
    maxAll(comm, (flag ? mypid : -1), problemPid);
    std::ostringstream buf; buf << problemPid;
    TEUCHOS_TEST_FOR_EXCEPTION(problemPid>-1, Exceptions::RuntimeError, "pid " + buf.str() + " encountered a partition number is that out-of-range");
    decompEntries = Teuchos::null;

    Teuchos::Array<GO> allPartitionsIContributeTo;
    Teuchos::Array<GO> localNnzPerPartition;
    hashTable->arrayify(allPartitionsIContributeTo, localNnzPerPartition);
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build target map", currentLevel));
    //map in which all pids have all partition numbers as GIDs.
    //FIXME this next map ctor can be a real time hog in parallel
    Array<GO> allPartitions;
    for (int i=0; i<numPartitions; ++i) allPartitions.push_back(i);
    RCP<Map> targetMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           numPartitions*comm->getSize(),
                                           allPartitions(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);

//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build vectors", currentLevel));
    RCP<Xpetra::Vector<double, LO, GO, NO> > globalWeightVec = Xpetra::VectorFactory<double, LO, GO, NO>::Build(targetMap);  //TODO why does the compiler grumble about this when I omit template arguments?
    RCP<Xpetra::Vector<LO, LO, GO, NO> > procWinnerVec = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(targetMap);
    ArrayRCP<LO> procWinner;
    if (procWinnerVec->getLocalLength() > 0)
      procWinner = procWinnerVec->getDataNonConst(0);
    for (int i=0; i<procWinner.size(); ++i) procWinner[i] = -1;
    procWinner = Teuchos::null;
    RCP<Xpetra::Vector<SC, LO, GO, NO> > scalarProcWinnerVec = Xpetra::VectorFactory<SC, LO, GO, NO>::Build(targetMap);

//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build unique map", currentLevel));
    Array<GO> myPidArray;
    myPidArray.push_back(mypid);

    // intermediate map required by ArbitrateAndCommunicate
    RCP<Map> uniqueMap = MapFactory::Build(decomposition->getMap()->lib(),
                                           comm->getSize(),
                                           myPidArray(),
                                           decomposition->getMap()->getIndexBase(),
                                           comm);
//m3 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: build comm helper", currentLevel));
    MueLu::UCAggregationCommHelper<LO, GO, NO, LMO> commHelper(uniqueMap, targetMap);
    myPartitionNumber = -1;
    int doArbitrate = 1;
//m3 = Teuchos::null;

    /*
       Use ArbitrateAndCommunicate to determine which process should own each partition.
       This may require multiple rounds because a single process i may be found to be the
       largest contributor for more than one partition (say P1, P2, ...Pn).  If that happens,
       process i is made owner of the partition Pk to which it contributes the most.  If
       two processes end up wanting to own the same partition, it's the job of A&C to break
       the tie.

       The contributions of process i to all other partitions are set to 0 (effectively meaning
       i can't be assigned another partition by A&C), and A&C is called again.

       continues until all partitions have been uniquely assigned.
    */

    m1 = Teuchos::null;

    int numRounds=0;
    RCP<SubFactoryMonitor> m2 = rcp(new SubFactoryMonitor(*this, "DeterminePartitionPlacement: Arbitration phase", currentLevel));
    while (doArbitrate)
    {
      ++numRounds;
      ArrayRCP<double> globalWeightVecData = globalWeightVec->getDataNonConst(0);

      //If this process doesn't yet own a partition, record all its nonzeros per partition as weights
      //If it doesn't contribute to a partition, make the weight small (0.1).  In this way, this pid
      //can become the owner of a partition if no one else can take it.
      if (myPartitionNumber == -1) {
        for (int i=0; i<globalWeightVecData.size(); ++i)
          globalWeightVecData[i] = 0.1;
        for (int i=0; i<allPartitionsIContributeTo.size(); ++i)
          globalWeightVecData[ allPartitionsIContributeTo[i] ] = localNnzPerPartition[i];
      } else {
        //this process already owns a partition, so record only the #nonzeros in that partition as weights
        bool noLocalDofsInMyPartition=true;
        for (int i=0; i<allPartitionsIContributeTo.size(); ++i) {
          if (allPartitionsIContributeTo[i] == myPartitionNumber) {
            globalWeightVecData[ myPartitionNumber ] = localNnzPerPartition[i];
            noLocalDofsInMyPartition = false;
            break;
          }
        }
        //In a previous round I was assigned a leftover partition.
        //Make sure I keep it by making the weight associated with it equal to 1.
        //All other PIDs will assign a weight of either 0 (because they own a partition already)
        //or 0.1 because they don't own a partition yet and don't contribute to this one (otherwise
        //they'd own this partition already).
        if (noLocalDofsInMyPartition) globalWeightVecData[ myPartitionNumber ] = 1;
      }
      globalWeightVecData = Teuchos::null;

      commHelper.ArbitrateAndCommunicate(*globalWeightVec, *procWinnerVec, NULL, true);

      if (procWinnerVec->getLocalLength() > 0)
        procWinner = procWinnerVec->getDataNonConst(0);

      // If this process has tentatively been assigned more than one partition,
      // choose the largest as the partition that this process will own.
      GO partitionsIOwn=0;
      GO largestPartitionSize=0;

      GO myPartitionNumberLastRound = myPartitionNumber;
      for (int i=0; i<procWinner.size(); ++i) {
        if (procWinner[i] == mypid) {
          GO partitionSize;
          //prefer partitions for which this pid has DOFs over partitions for which it doesn't.
          if (hashTable->containsKey(i)) partitionSize = hashTable->get(i);
          else                           partitionSize = 1;
          if (partitionSize > largestPartitionSize) {
            myPartitionNumber = i;
            largestPartitionSize = partitionSize;
          }
          partitionsIOwn++;
        }
      }
      procWinner = Teuchos::null;
      //Check to see if my newly assigned partition is one for which I have no DOFs.
      bool gotALeftoverPartitionThisRound=false;
      if (myPartitionNumber > -1 && (myPartitionNumber != myPartitionNumberLastRound)) {
        if (hashTable->containsKey(myPartitionNumber)) gotALeftoverPartitionThisRound = false;
        else                                           gotALeftoverPartitionThisRound = true;
      }

      // If any pid got a leftover partition this round, or if this pid was tentatively assigned
      // more than one partition, then we need to arbitrate again, because either there are unassigned
      // partitions, or there are more than one pid claiming ownership of the same partition.
      int arbitrateAgain;
      if (partitionsIOwn > 1 || gotALeftoverPartitionThisRound) arbitrateAgain=1;
      else                                                      arbitrateAgain=0;
      sumAll(comm, arbitrateAgain, doArbitrate);

    } //while (doArbitrate)
    m2 = Teuchos::null;

    GetOStream(Statistics0, 0) << "Number of arbitration rounds = " << numRounds << std::endl;

    ArrayRCP<const LO> procWinnerConst;
    if (procWinnerVec->getLocalLength() > 0)
      procWinnerConst = procWinnerVec->getData(0);

    int numPartitionOwners=0;
    for (int i=0; i<procWinnerConst.size(); ++i)
      if (procWinnerConst[i] > -1) ++numPartitionOwners;
    buf << numPartitionOwners;
    TEUCHOS_TEST_FOR_EXCEPTION(numPartitionOwners != numPartitions, Exceptions::RuntimeError,
                               "Number of partition owners (" + buf.str() + ") is not equal to number of partitions");
    partitionOwners = procWinnerConst(); //only works if procWinner is const ...

    // print the grid of processors
    GetOStream(Statistics0, 0) << "Partition distribution over cores, + indicates partition ownership" << std::endl;
    int numProc = comm->getSize();
    ArrayRCP<char> grid(numProc, '.');
    for (int i=0; i<partitionOwners.size(); ++i) grid[ partitionOwners[i] ] = '+';
    int sizeOfARow = (int) sqrt(numProc);
    int numRows = numProc / sizeOfARow;
    int leftOvers = numProc  - (numProc/sizeOfARow)*sizeOfARow;
    int ctr=0;
    int pidCtr=0;
    for (int i=0; i<numRows; ++i) {
      for (int j=0; j<sizeOfARow; ++j)
        GetOStream(Statistics0, 0) << grid[ctr++];
      GetOStream(Statistics0, 0) << "      " << pidCtr << ":" << pidCtr+sizeOfARow-1 << std::endl;;
      pidCtr += sizeOfARow;
    }
    if (leftOvers > 0) {
      for (int i=0; i<leftOvers; ++i)
        GetOStream(Statistics0, 0) << grid[ctr++];

      Array<char> aos(sizeOfARow-leftOvers, ' ');
      std::string spaces(aos.begin(), aos.end());
      GetOStream(Statistics0, 0) << spaces << "      " << pidCtr << ":" << pidCtr+leftOvers-1 << std::endl;;
    }

  } //DeterminePartitionPlacement

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetStartLevel(int startLevel) {
    startLevel_ = startLevel;
  }

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetImbalanceThreshold(double threshold) {
    nnzMaxMinRatio_ = threshold;
  }

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMinRowsPerProcessor(GO threshold) {
    minRowsPerProcessor_ = threshold;
  }

  //----------------------------------------------------------------------

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMinNnzPerProcessor(GO threshold) {
    minNnzPerProcessor_ = threshold;
  }

} // namespace MueLu

#endif //ifdef HAVE_MPI

#endif // MUELU_REPARTITIONFACTORY_DEF_HPP
