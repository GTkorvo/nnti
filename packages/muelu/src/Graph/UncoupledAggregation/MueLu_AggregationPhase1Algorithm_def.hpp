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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_

#include <queue>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase1Algorithm_decl.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat,
                  LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    AggOptions::Ordering ordering    = params.get<AggOptions::Ordering>("Ordering");
    LO maxNeighAlreadySelected       = params.get<LO>                  ("MaxNeighAlreadySelected");
    LO minNodesPerAggregate          = params.get<LO>                  ("MinNodesPerAggregate");
    LO maxNodesPerAggregate          = params.get<LO>                  ("MaxNodesPerAggregate");

    TEUCHOS_TEST_FOR_EXCEPTION(maxNodesPerAggregate < minNodesPerAggregate, Exceptions::RuntimeError,
                               "MueLu::UncoupledAggregationAlgorithm::BuildAggregates: minNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

    if (ordering != NATURAL && ordering != RANDOM && ordering != GRAPH)
      throw Exceptions::RuntimeError("UncoupledAggregation::BuildAggregates : bad aggregation ordering option");

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);

    LO numLocalAggregates = aggregates.GetNumAggregates();

    ArrayRCP<LO> randomVector;
    if (ordering == RANDOM) {
      randomVector = arcp<LO>(numRows);
      for (LO i = 0; i < numRows; i++)
        randomVector[i] = i;
      RandomReorder(randomVector);
    }

    int              aggIndex = -1;
    size_t           aggSize  =  0;
    std::vector<int> aggList(graph.getNodeMaxNumRowEntries());

    std::queue<LO>   graphOrderQueue;

    // Main loop over all local rows of graph(A)
    for (LO i = 0; i < numRows; i++) {
      // Step 1: pick the next node to aggregate
      LO rootCandidate = 0;
      if      (ordering == NATURAL) rootCandidate = i;
      else if (ordering == RANDOM)  rootCandidate = randomVector[i];
      else if (ordering == GRAPH) {

        if (graphOrderQueue.size() == 0) {
          // Current queue is empty for GRAPH ordering, populate with one READY node
          for (LO jnode = 0; jnode < numRows; jnode++)
            if (aggStat[jnode] == READY) {
              graphOrderQueue.push(jnode);
              break;
            }
        }
        if (graphOrderQueue.size() == 0) {
          // There are no more ready nodes, end the phase
          break;
        }
        rootCandidate = graphOrderQueue.front();   // take next node from graph ordering queue
        graphOrderQueue.pop();                     // delete this node in list
      }

      if (aggStat[rootCandidate] != READY)
        continue;

      // Step 2: build tentative aggregate
      aggSize = 0;
      aggList[aggSize++] = rootCandidate;

      ArrayView<const LO> neighOfINode = graph.getNeighborVertices(rootCandidate);

      // If the number of neighbors is less than the minimum number of nodes
      // per aggregate, we know this is not going to be a valid root, and we
      // may skip it, but only for NATURAL and RANDOM (for GRAPH we still need
      // to fetch the list of local neighbors to continue)
      if ((ordering == NATURAL || ordering == RANDOM) &&
          neighOfINode.size() < minNodesPerAggregate) {
        continue;
      }

      LO numAggregatedNeighbours = 0;

      for (int j = 0; j < neighOfINode.size(); j++) {
        LO neigh = neighOfINode[j];

        if (neigh != rootCandidate && graph.isLocalNeighborVertex(neigh)) {

          if (aggStat[neigh] == READY || aggStat[neigh] == NOTSEL) {
            // If aggregate size does not exceed max size, add node to the
            // tentative aggregate
            // NOTE: We do not exit the loop over all neighbours since we have
            // still to count all aggregated neighbour nodes for the
            // aggregation criteria
            // NOTE: We check here for the maximum aggregation size. If we
            // would do it below with all the other check too big aggregates
            // would not be accepted at all.
            if (aggSize < as<size_t>(maxNodesPerAggregate))
              aggList[aggSize++] = neigh;

          } else {
            numAggregatedNeighbours++;
          }
        }
      }

      // Step 3: check if tentative aggregate is acceptable
      if ((numAggregatedNeighbours <= maxNeighAlreadySelected) &&   // too many connections to other aggregates
          (aggSize                 >= minNodesPerAggregate)) {      // too few nodes in the tentative aggregate
        // Accept new aggregate
        // rootCandidate becomes the root of the newly formed aggregate
        aggregates.SetIsRoot(rootCandidate);
        aggIndex = numLocalAggregates++;

        for (size_t k = 0; k < aggSize; k++) {
          aggStat     [aggList[k]] = AGGREGATED;
          vertex2AggId[aggList[k]] = aggIndex;
          procWinner  [aggList[k]] = myRank;

          if (ordering == GRAPH) {
            ArrayView<const LO> neighOfJNode = graph.getNeighborVertices(aggList[k]);
            for (int j = 0; j < neighOfJNode.size(); j++) {
              LO neigh = neighOfJNode[j];

              if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == READY)
                graphOrderQueue.push(neigh);
            }
          }
        }

        numNonAggregatedNodes -= aggSize;

      } else {
        // Aggregate is not accepted
        aggStat[rootCandidate] = NOTSEL;

        if (ordering == GRAPH) {
          // Even though the aggregate around rootCandidate is not perfect, we want to try
          // the neighbor nodes of rootCandidate
          for (int j = 0; j < neighOfINode.size(); j++) {
            LO neigh = neighOfINode[j];

            if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == READY)
              graphOrderQueue.push(neigh);
          }
        }
      }
    }

    // Reset all NOTSEL vertices to READY
    // This simplifies other algorithms
    for (LO i = 0; i < numRows; i++)
      if (aggStat[i] == NOTSEL)
        aggStat[i] = READY;

    // update aggregate object
    aggregates.SetNumAggregates(numLocalAggregates);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomReorder(ArrayRCP<LO> list) const {
    //TODO: replace int
    int n = list.size();
    for(int i = 0; i < n-1; i++)
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomOrdinal(int min, int max) const {
    return min + as<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

} // end namespace


#endif /* MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_ */
