/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"

#ifdef IFPACK_NODE_AWARE_CODE
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_Hashtable.hpp"
#include "Teuchos_Array.hpp"
#include "EpetraExt_OperatorOut.h"
extern int ML_NODE_ID;
#endif

using namespace Teuchos;

// ====================================================================== 
// Constructor for the case of one core per subdomain
Ifpack_OverlappingRowMatrix::
Ifpack_OverlappingRowMatrix(const RCP<const Epetra_RowMatrix>& Matrix_in,
                            int OverlapLevel_in)  :
  Matrix_(Matrix_in),
  OverlapLevel_(OverlapLevel_in)
{
  // should not be here if no overlap
  if (OverlapLevel_in == 0)
    IFPACK_CHK_ERRV(-1);

  // nothing to do as well with one process
  if (Comm().NumProc() == 1)
    IFPACK_CHK_ERRV(-1);
  
  NumMyRowsA_ = A().NumMyRows();

  // construct the external matrix
  vector<int> ExtElements; 

  RCP<Epetra_Map> TmpMap;
  RCP<Epetra_CrsMatrix> TmpMatrix; 
  RCP<Epetra_Import> TmpImporter;

  // importing rows corresponding to elements that are 
  // in ColMap, but not in RowMap 
  const Epetra_Map *RowMap; 
  const Epetra_Map *ColMap; 

  for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap) {
    if (TmpMatrix != Teuchos::null) {
      RowMap = &(TmpMatrix->RowMatrixRowMap()); 
      ColMap = &(TmpMatrix->RowMatrixColMap()); 
    }
    else {
      RowMap = &(A().RowMatrixRowMap()); 
      ColMap = &(A().RowMatrixColMap()); 
    }

    int size = ColMap->NumMyElements() - RowMap->NumMyElements(); 
    vector<int> list(size); 

    int count = 0; 

    // define the set of rows that are in ColMap but not in RowMap 
    for (int i = 0 ; i < ColMap->NumMyElements() ; ++i) { 
      int GID = ColMap->GID(i); 
      if (A().RowMatrixRowMap().LID(GID) == -1) { 
        vector<int>::iterator pos 
          = find(ExtElements.begin(),ExtElements.end(),GID); 
        if (pos == ExtElements.end()) { 
          ExtElements.push_back(GID);
          list[count] = GID; 
          ++count; 
        } 
      } 
    } 

    TmpMap = rcp( new Epetra_Map(-1,count, &list[0],0,Comm()) ); 

    TmpMatrix = rcp( new Epetra_CrsMatrix(Copy,*TmpMap,0) ); 

    TmpImporter = rcp( new Epetra_Import(*TmpMap,A().RowMatrixRowMap()) ); 

    TmpMatrix->Import(A(),*TmpImporter,Insert); 
    TmpMatrix->FillComplete(A().OperatorDomainMap(),*TmpMap); 

  }

  // build the map containing all the nodes (original
  // matrix + extended matrix)
  vector<int> list(NumMyRowsA_ + ExtElements.size());
  for (int i = 0 ; i < NumMyRowsA_ ; ++i)
    list[i] = A().RowMatrixRowMap().GID(i);
  for (int i = 0 ; i < (int)ExtElements.size() ; ++i)
    list[i + NumMyRowsA_] = ExtElements[i];

  Map_ = rcp( new Epetra_Map(-1, NumMyRowsA_ + ExtElements.size(),
	  			      &list[0], 0, Comm()) );
# ifdef IFPACK_NODE_AWARE_CODE
  colMap_ = Map_;
# endif
  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = rcp( new Epetra_Map(-1,ExtElements.size(),
					 &ExtElements[0],0,A().Comm()) );
  ExtMatrix_ = rcp( new Epetra_CrsMatrix(Copy,*ExtMap_,*Map_,0) ); 

  ExtImporter_ = rcp( new Epetra_Import(*ExtMap_,A().RowMatrixRowMap()) ); 
  ExtMatrix_->Import(A(),*ExtImporter_,Insert); 
  ExtMatrix_->FillComplete(A().OperatorDomainMap(),*Map_);

  Importer_ = rcp( new Epetra_Import(*Map_,A().RowMatrixRowMap()) );

  // fix indices for overlapping matrix
  NumMyRowsB_ = B().NumMyRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;
  NumMyCols_ = NumMyRows_;
  
  NumMyDiagonals_ = A().NumMyDiagonals() + B().NumMyDiagonals();
  
  NumMyNonzeros_ = A().NumMyNonzeros() + B().NumMyNonzeros();
  Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = A().MaxNumEntries();
  
  if (MaxNumEntries_ < B().MaxNumEntries())
    MaxNumEntries_ = B().MaxNumEntries();

}

#ifdef IFPACK_NODE_AWARE_CODE
// ====================================================================== 
// Constructor for the case of two or more cores per subdomain
Ifpack_OverlappingRowMatrix::
Ifpack_OverlappingRowMatrix(const RCP<const Epetra_RowMatrix>& Matrix_in,
                            int OverlapLevel_in, int nodeID)  :
  Matrix_(Matrix_in),
  OverlapLevel_(OverlapLevel_in)
{

  // should not be here if no overlap
  if (OverlapLevel_in == 0)
    IFPACK_CHK_ERRV(-1);

  // nothing to do as well with one process
  if (Comm().NumProc() == 1)
    IFPACK_CHK_ERRV(-1);

  // nodeID is the node (aka socket) number, and so is system dependent
  // nodeComm is the communicator for all the processes on a particular node
  // these processes will have the same nodeID.
# ifdef HAVE_MPI
  const Epetra_MpiComm *pComm = dynamic_cast<const Epetra_MpiComm*>( &Comm() );
  assert(pComm != NULL);
  MPI_Comm nodeMPIComm;
  MPI_Comm_split(pComm->Comm(),nodeID,pComm->MyPID(),&nodeMPIComm);
  Epetra_MpiComm *nodeComm = new Epetra_MpiComm(nodeMPIComm);
# else
  Epetra_SerialComm *nodeComm =  dynamic_cast<Epetra_MpiComm*>( &(Matrix_in->RowMatrixRowMap().Comm()) );
# endif
  
  NumMyRowsA_ = A().NumMyRows();

/*
  if (A().Comm().MyPID() == 0) printf("====== globA rowmap ========\n");
  cout << A().RowMatrixRowMap() << endl;
  if (A().Comm().MyPID() == 0) printf("====== globA colmap ========\n");
  cout << A().RowMatrixColMap() << endl;
*/

  // off-node GIDs that will be in the overlap
  vector<int> ghostElements; 

  // GIDs that will be in the overlapped matrix's column map
  vector<int> colMapElements; 

  // ghostTable holds off-node GIDs that are connected to on-node rows and can potentially be this PID's overlap
  // TODO hopefully 3 times the # column entries is a reasonable table size
  Teuchos::Hashtable<int,int> ghostTable;

  // nbTable holds node buddy GIDs that are connected to current PID's rows, i.e., GIDs that should be in the overlapped
  // matrix's column map
  Teuchos::Hashtable<int,int> colMapTable;

  RCP<Epetra_Map> TmpMap;
  RCP<Epetra_CrsMatrix> TmpMatrix; 
  RCP<Epetra_Import> TmpImporter;

  // importing rows corresponding to elements that are 
  // in ColMap, but not in RowMap 
  const Epetra_Map *RowMap; 
  const Epetra_Map *ColMap; 
  const Epetra_Map *DomainMap;

  int mypid = Comm().MyPID();

  // TODO Count #connections from nodes I own to each ghost node

/*
  if (!mypid) printf("original matrix colmap\n===========================\n");
  fflush(stdout);
  sleep(1);
  cout << A().RowMatrixColMap() << endl;
  if (!mypid) printf("===========================\n");
  sleep(1);
  Comm().Barrier();
*/

  /* ** ************************************************************************** ** */
  /* ** ********************** start of main overlap loop ************************ ** */
  /* ** ************************************************************************** ** */
  for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap)
  {
    if (TmpMatrix != Teuchos::null) {
      RowMap = &(TmpMatrix->RowMatrixRowMap()); 
      ColMap = &(TmpMatrix->RowMatrixColMap()); 
      DomainMap = &(TmpMatrix->OperatorDomainMap());
    }
    else {
      RowMap = &(A().RowMatrixRowMap()); 
      ColMap = &(A().RowMatrixColMap()); 
      DomainMap = &(A().OperatorDomainMap());
    }

    ghostTable = Teuchos::Hashtable<int,int>(3 * A().RowMatrixColMap().NumMyElements() );
    colMapTable = Teuchos::Hashtable<int,int>(3 * A().RowMatrixColMap().NumMyElements() );

    // For each column ID, determine the owning node (as opposed to core)
    // ID of the corresponding row.
    Epetra_IntVector colIdList( *ColMap );
    Epetra_IntVector rowIdList(*DomainMap);
    rowIdList.PutValue(nodeID);  
    Teuchos::RCP<Epetra_Import> nodeIdImporter = rcp(new Epetra_Import( *ColMap, *DomainMap ));
    colIdList.Import(rowIdList,*nodeIdImporter,Insert);

    int size = ColMap->NumMyElements() - RowMap->NumMyElements(); 
    vector<int> list(size); 
    int count = 0; 

    // define the set of off-node rows that are in ColMap but not in RowMap
    // This naturally does not include off-core rows that are on the same node as me, i.e., node buddy rows.
    for (int i = 0 ; i < ColMap->NumMyElements() ; ++i) {
      int GID = ColMap->GID(i); 
      if ( colIdList[i] != nodeID )
      {
        int votes;
        if (ghostTable.containsKey(GID)) {
          votes = ghostTable.get(GID);
          votes++;
          ghostTable.put(GID,votes);
        } else {
          ghostTable.put(GID,1);
        }
      }
    } //for (int i = 0 ; i < ColMap->NumMyElements() ; ++i)

    Teuchos::Array<int> gidsarray,votesarray;
    ghostTable.arrayify(gidsarray,votesarray);
    int *gids = gidsarray.getRawPtr();
    int *votes = votesarray.getRawPtr();

    /*
       This next bit of code decides which node buddy (NB) gets which ghost points.  Everyone sends their
       list of ghost points to pid 0 of the local subcommunicator.  Pid 0 decides who gets what:

          - if a ghost point is touched by only one NB, that NB gets the ghost point
          - if two or more NBs share a ghost point, whichever NB has the most connections to the ghost
            point gets it.
    */

#   ifdef HAVE_MPI  //FIXME What if we build in serial?!  This file will likely not compile.
    int lengths[nodeComm->NumProc()+1];
    /*
      lengths[i]                    starting position for pid i's entries in ghosts, round, and owningpid
      lengths[i+1] - lengths[i]     #entries belonging to pid i in said vectors
    */
    int *cullList;
    int ncull;
    int mypid = nodeComm->MyPID();

    if (nodeComm->MyPID() == 0)
    {
      // Figure out how much pid 0 is to receive
      MPI_Status status;
      lengths[0] = 0;
      lengths[1] = ghostTable.size();
      for (int i=1; i<nodeComm->NumProc(); i++) {
        int leng;
        MPI_Recv( &leng, 1, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
        lengths[i+1] = lengths[i] + leng;
      }
      int total = lengths[nodeComm->NumProc()];

      int* ghosts = new int[total];
      for (int i=0; i<total; i++) ghosts[i] = -9;
      int *round  = new int[total];
      int *owningpid  = new int[total];

      for (int i=1; i<nodeComm->NumProc(); i++) {
        int count = lengths[i+1] - lengths[i];
        MPI_Recv( ghosts+lengths[i], count, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
        MPI_Recv( round+lengths[i], count, MPI_INT, i, MPI_ANY_TAG, nodeComm->Comm(), &status);
      }

      // put in pid 0's info
      for (int i=0; i<lengths[1]; i++) {
        ghosts[i] = gids[i];
        round[i] = votes[i];
        owningpid[i] = 0;
      }

      // put in the pid associated with each ghost
      for (int j=1; j<nodeComm->NumProc(); j++) {
        for (int i=lengths[j]; i<lengths[j+1]; i++) {
          owningpid[i] = j;
        }
      }

      // sort everything based on the ghost gids
      int* roundpid[2];
      roundpid[0] = round; roundpid[1] = owningpid;
      Epetra_Util epetraUtil;
      epetraUtil.Sort(true,total,ghosts,0,0,2,roundpid);

      //set up arrays that get sent back to node buddies and that tell them what ghosts to cull
      int nlosers[nodeComm->NumProc()];
      int* losers[nodeComm->NumProc()];
      for (int i=0; i<nodeComm->NumProc(); i++) {
        nlosers[i] = 0;
        losers[i] = new int[ lengths[i+1]-lengths[i] ];
      }

      // Walk through ghosts array and and for each sequence of duplicate ghost GIDs, choose just one NB to keep it.
      // The logic is pretty simple.  The ghost list is sorted, so all duplicate PIDs are together.
      // The list is traversed.  As duplicates are found, node pid 0 keeps track of the current "winning"
      // pid.  When a pid is determined to have "lost" (less votes/connections to the current GID), the
      // GID is added to that pid's list of GIDs to be culled.  At the end of the repeated sequence, we have
      // a winner, and other NBs know whether they need to delete it from their import list.
      int max=0;   //for duplicated ghosts, index of pid with most votes and who hence keeps the ghost.
                   // TODO to break ties randomly

      for (int i=1; i<total; i++) {
        if (ghosts[i] == ghosts[i-1]) {
          int idx = i; // pid associated with idx is current "loser"
          if (round[i] > round[max]) {
            idx = max;
            max=i;
          }
          int j = owningpid[idx];
          losers[j][nlosers[j]++] = ghosts[idx];
        } else {
          max=i;
        }
      } //for (int i=1; i<total; i++)

      delete [] round;
      delete [] ghosts;
      delete [] owningpid;

      // send the arrays of ghost GIDs to be culled back to the respective node buddies
      for (int i=1; i<nodeComm->NumProc(); i++) {
        MPI_Send( nlosers+i, 1, MPI_INT, i, 8675, nodeComm->Comm());
        MPI_Send( losers[i], nlosers[i], MPI_INT, i, 8675, nodeComm->Comm());
      }

      //FIXME Unnecessary memory allocation and copying, but makes culling code cleaner
      //Could we stick this info into gids and votes, since neither is being used anymore?
      //TODO Instead of using "losers" arrays, just use "cullList" as in the else clause
      ncull = nlosers[0];
      cullList = new int[ncull+1];
      for (int i=0; i<nlosers[0]; i++)
        cullList[i] = losers[0][i];

      for (int i=0; i<nodeComm->NumProc(); i++)
        delete [] losers[i];

    } else { //everyone but pid 0

      // send to node pid 0 all ghosts that this pid could potentially import
      int hashsize = ghostTable.size();
      MPI_Send( &hashsize, 1, MPI_INT, 0, 8675, nodeComm->Comm());
      MPI_Send( gids, hashsize, MPI_INT, 0, 8675, nodeComm->Comm());
      MPI_Send( votes, hashsize, MPI_INT, 0, 8675, nodeComm->Comm());

      // receive the ghost GIDs that should not be imported (subset of the list sent off just above)
      MPI_Status status;
      MPI_Recv( &ncull, 1, MPI_INT, 0, 8675, nodeComm->Comm(), &status);
      cullList = new int[ncull+1];
      MPI_Recv( cullList, ncull, MPI_INT, 0, 8675, nodeComm->Comm(), &status);
    }

    //TODO clean out hash table after each time through overlap loop   4/1/07 JJH done moved both hash tables to local scope

    // Remove from my row map all off-node ghosts that will be imported by a node buddy.
    for (int i=0; i<ncull; i++) {
      try{ghostTable.remove(cullList[i]);}

      catch(...) {
        printf("pid %d: In OverlappingRowMatrix ctr, problem removing ghost elt %d from ghostTable\n",
               Comm().MyPID(),cullList[i]);
        fflush(stdout);
      }
    }//for

    delete [] cullList;

    // Save off the remaining ghost GIDs from the current overlap round.
    // These are off-node GIDs (rows) that I will import.
    gidsarray.clear(); votesarray.clear();
    ghostTable.arrayify(gidsarray,votesarray);

    count=0;
    for (int i=0; i<ghostTable.size(); i++) {
      ghostElements.push_back(gidsarray[i]);
      list[i] = gidsarray[i];  //FIXME this won't work for >1 level of overlap. list should only contain
                               //FIXME the *current* level of overlap, not *all* the overlap
                               // JJH 4/1/09 this is ok -- list is local scope, so it can *only* have current level of overlap
      count++;
    }

#   endif //ifdef HAVE_MPI

    TmpMap = rcp( new Epetra_Map(-1,count, &list[0],0,Comm()) );

    TmpMatrix = rcp( new Epetra_CrsMatrix(Copy,*TmpMap,0) ); 

    TmpImporter = rcp( new Epetra_Import(*TmpMap,A().RowMatrixRowMap()) ); 

    TmpMatrix->Import(A(),*TmpImporter,Insert); 
    TmpMatrix->FillComplete(A().OperatorDomainMap(),*TmpMap); 

    // These next two imports get the GIDs that need to go into the column map of the overlapped matrix.

    // For each column ID in the overlap, determine the owning node (as opposed to core)
    // ID of the corresponding row.  Save those column IDs whose owning node is the current one.
    // This should get all the imported ghost GIDs.
    Epetra_IntVector ov_colIdList( TmpMatrix->ColMap() );
    ov_colIdList.PutValue(-1);
    Epetra_IntVector ov_rowIdList( TmpMatrix->RowMap() );
    ov_rowIdList.PutValue(nodeID);  
    Teuchos::RCP<Epetra_Import> ov_nodeIdImporter = rcp(new Epetra_Import( TmpMatrix->ColMap(), TmpMatrix->RowMap()));
    ov_colIdList.Import(ov_rowIdList,*ov_nodeIdImporter,Insert);

    for (int i=0 ; i<ov_colIdList.MyLength(); i++) {
       if (ov_colIdList[i] == nodeID) {
         int GID = (ov_colIdList.Map()).GID(i);
         colMapTable.put(GID,1);
       }
    }

    // Do a second import of the owning node ID from A's rowmap to TmpMat's column map.  This ensures that
    // all GIDs that belong to a node buddy and are in a ghost row's sparsity pattern will be in the final
    // overlapped matrix's column map.
    ov_colIdList.PutValue(-1);
    Epetra_IntVector ArowIdList( A().RowMatrixRowMap() );
    ArowIdList.PutValue(nodeID);
    nodeIdImporter = rcp(new Epetra_Import( TmpMatrix->ColMap(), A().RowMatrixRowMap() ));
    ov_colIdList.Import(ArowIdList,*nodeIdImporter,Insert);

    for (int i=0 ; i<ov_colIdList.MyLength(); i++) {
       if (ov_colIdList[i] == nodeID) {
         int GID = (ov_colIdList.Map()).GID(i);
         colMapTable.put(GID,1);
       }
    }

    for (int i = 0 ; i < A().RowMatrixColMap().NumMyElements() ; ++i) {
      int GID = ColMap->GID(i); 
      // Remove any entries that are in A's original column map
      if (colMapTable.containsKey(GID)) {
        try{colMapTable.remove(GID);}
        catch(...) {
          printf("pid %d: In OverlappingRowMatrix ctr, problem removing entry %d from colMapTable\n", Comm().MyPID(),GID);
          fflush(stdout);
        }
      }
    }

    gidsarray.clear(); votesarray.clear();
    colMapTable.arrayify(gidsarray,votesarray);
    for (int i=0; i<colMapTable.size(); i++)
      colMapElements.push_back(gidsarray[i]);

  } //for (int overlap = 0 ; overlap < OverlapLevel_in ; ++overlap)

  /* ** ************************************************************************ ** */
  /* ** ********************** end of main overlap loop ************************ ** */
  /* ** ************************************************************************ ** */

/*
   We need two maps here.  The first is the row map, which we've got by using my original rows
   plus whatever I've picked up in ghostElements.

   The second is a column map.  This map should include my rows, plus any columns that could come from node buddies.
   These GIDs come from the std:array colMapElements, which in turn comes from colMapTable.
   This map should *not* omit ghosts that have been imported by my node buddies, i.e., for any row that I own,
   the stencil should include all column GIDs (including imported ghosts) that are on the node.
*/

  // build the row map containing all the nodes (original matrix + off-node matrix)
  vector<int> rowList(NumMyRowsA_ + ghostElements.size());
  for (int i = 0 ; i < NumMyRowsA_ ; ++i)
    rowList[i] = A().RowMatrixRowMap().GID(i);
  for (int i = 0 ; i < (int)ghostElements.size() ; ++i)
    rowList[i + NumMyRowsA_] = ghostElements[i];

  // row map for the overlapped matrix (local + overlap)
  Map_ = rcp( new Epetra_Map(-1, NumMyRowsA_ + ghostElements.size(), &rowList[0], 0, Comm()) );

  // build the column map for the overlapping matrix
  //vector<int> colList(colMapElements.size());
  // column map for the overlapped matrix (local + overlap)
  //colMap_ = rcp( new Epetra_Map(-1, colMapElements.size(), &colList[0], 0, Comm()) );
  //for (int i = 0 ; i < (int)colMapElements.size() ; i++)
  //  colList[i] = colMapElements[i];
  vector<int> colList(A().RowMatrixColMap().NumMyElements() + colMapElements.size());
  int nc = A().RowMatrixColMap().NumMyElements();
  for (int i = 0 ; i < nc; i++)
    colList[i] = A().RowMatrixColMap().GID(i);
  for (int i = 0 ; i < (int)colMapElements.size() ; i++)
    colList[nc+i] = colMapElements[i];
  // column map for the overlapped matrix (local + overlap)
  colMap_ = rcp( new Epetra_Map(-1, A().RowMatrixColMap().NumMyElements() + colMapElements.size(), &colList[0], 0, Comm()) );

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++ start of sort
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // build the column map, but don't use a copy constructor, b/c local communicator SubComm_ is
  // different from that of Matrix.
  try {
    // build row map based on the local communicator.  We need this temporarily to build the column map.
    RCP<Epetra_Map> nodeMap_ = rcp( new Epetra_Map(-1,NumMyRowsA_ + ghostElements.size(),&rowList[0],0,*nodeComm) );
    //Epetra_Map* nodeMap_ = new Epetra_Map(-1,NumMyRowsA_ + ghostElements.size(),&rowList[0],0,*nodeComm) ;
    int numMyElts = colMap_->NumMyElements();

    Teuchos::RCP<int> myGlobalElts = rcp( new int[numMyElts] );
    colMap_->MyGlobalElements(&*myGlobalElts);
    Teuchos::RCP<int> pidList = rcp( new int[numMyElts] );
    nodeMap_->RemoteIDList(numMyElts, &*myGlobalElts, &*pidList, 0);

   /* The column map *must* be sorted: first locals, then ghosts.
      The ghosts must be further sorted so that they are contiguous by owning processor.  */

    // first sort on the owning pid in the *local* communicator
    Epetra_Util Util;
    int *tt[1];
    tt[0] = &*myGlobalElts;
    Util.Sort(true, numMyElts, &*pidList, 0, (double**)0, 1, tt);

    // for each remote pid, sort the entries in ascending order
    // don't sort the local pid's entries
    int localStart=0;
    while (localStart<numMyElts) {
      int currPID = (&*pidList)[localStart];
      int i=localStart;
      while (i<numMyElts && (&*pidList)[i] == currPID) i++;
      if (currPID != nodeComm->MyPID())
        Util.Sort(true, i-localStart, (&*myGlobalElts)+localStart, 0, 0, 0, 0);
      localStart = i;
    }

    // now move the local entries to the front of the list
    localStart=0;
    while (localStart<numMyElts && (&*pidList)[localStart] != nodeComm->MyPID()) localStart++;
    assert(localStart != numMyElts);
    int localEnd=localStart;
    while (localEnd<numMyElts && (&*pidList)[localEnd] == nodeComm->MyPID()) localEnd++;
    assert(numMyElts!=0);
    Teuchos::RCP<int> mySortedGlobalElts = rcp( new int[numMyElts] );
    //Teuchos::RCP<int> mySortedPidList = rcp( new int[numMyElts] );
    RCP<int> mySortedPidList = rcp( new int[numMyElts] );
    int leng = localEnd - localStart;
    /* This is a little gotcha.  It appears that the ordering of the column map's local entries
       must be the same as that of the domain map's local entries.  See the comment in method
       MakeColMap() in Epetra_CrsGraph.cpp, line 1072. */
    int *rowGlobalElts =  nodeMap_->MyGlobalElements();
    int numRowElts = nodeMap_->NumMyElements();
    //printf("lpid %d: numRowElts = %d, leng = %d\n",Comm().MyPID(),numRowElts,leng); fflush(stdout);
    //assert(numRowElts==leng);
    for (int i=0; i<leng; i++) {
      (&*mySortedGlobalElts)[i] = rowGlobalElts[i];
      (&*mySortedPidList)[i] = nodeComm->MyPID();
/*
      (&*mySortedGlobalElts)[i] = (&*myGlobalElts)[localStart+i];
      (&*mySortedPidList)[i] = (&*pidList)[localStart+i];
*/
    }
    for (int i=leng; i< localEnd; i++) {
      (&*mySortedGlobalElts)[i] = (&*myGlobalElts)[i-leng];
      (&*mySortedPidList)[i] = (&*pidList)[i-leng];
    }
    for (int i=localEnd; i<numMyElts; i++) {
      (&*mySortedGlobalElts)[i] = (&*myGlobalElts)[i];
      (&*mySortedPidList)[i] = (&*pidList)[i];
    }

    int indexBase = colMap_->IndexBase();
    colMap_ = Teuchos::null;
    colMap_ = rcp( new Epetra_Map(-1,numMyElts,&*mySortedGlobalElts,indexBase,Comm()) );

/*
    if (SubComm_->MyPID()==0)
      printf(">>>>>> node %d local column map <<<<<<<\n",ML_NODE_ID);
    cout << *(&*colMap_) << endl;

    sleep(1);

*/
//EPETRA_CHK_ERR(EpetraExt::OperatorToMatlabFile("Ainv.dat", AInvOp));

  }
  catch(...) {
    printf("** * Ifpack_OverlappingRowmatrix ctor: problem creating column map * **\n\n");
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++ end of sort
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Teuchos::RCP<int> pidList = rcp( new int[colMapElements.size()] );

  Map_->RemoteIDList(colMapElements.size(), &colList[0], &*pidList, 0);

/*
  printf("++++++++\ngpid %d\n++++++++\n",Comm().MyPID()); 
  for (int i=0; i<colMapElements.size(); i++) printf("   %d:  %d\n",colList[i],(&*pidList)[i]);
  fflush(stdout);
*/

/*

   FIXME
   Does the column map need to be sorted for the overlapping matrix?

   The column map *must* be sorted:

        first locals
        then ghosts

   The ghosts must be further sorted so that they are contiguous by owning processor

  int* RemoteSizeList = 0
  int* RemoteColIndices = ColIndices.Values() + NumLocalColGIDs; // Points to back end of ColIndices

  EPETRA_CHK_ERR(DomainMap.RemoteIDList(NumRemoteColGIDs, RemoteColIndices, PIDList.Values(), 0, RemoteSizeList));
  Epetra_Util epetraUtil;
  SortLists[0] = RemoteColIndices;
  SortLists[1] = RemoteSizeList;
  epetraUtil.Sort(true, NumRemoteColGIDs, PIDList.Values(), 0, 0, NLists, SortLists);
*/

/*
  if (!mypid)
    printf("row map\n==========================\n");
  fflush(stdout);sleep(1);
  cout << *Map_ << endl;
  sleep(1);
  Comm().Barrier();

  if (!mypid)
    printf("overlapped col map\n===================\n");
  fflush(stdout);
  sleep(1);
  cout << *colMap_ << endl;
  sleep(3);
*/

  // now build the map corresponding to all the external nodes
  // (with respect to A().RowMatrixRowMap().
  ExtMap_ = rcp( new Epetra_Map(-1,ghostElements.size(), &ghostElements[0],0,Comm()) );
  ExtMatrix_ = rcp( new Epetra_CrsMatrix(Copy,*ExtMap_,*colMap_,0) ); 

  ExtImporter_ = rcp( new Epetra_Import(*ExtMap_,A().RowMatrixRowMap()) ); 
  ExtMatrix_->Import(A(),*ExtImporter_,Insert); 

/*
  Notes to self:    In FillComplete, the range map does not have to be 1-1 as long as
                    (row map == range map).  Ditto for the domain map not being 1-1
                    if (col map == domain map).
                    
*/

  ExtMatrix_->FillComplete( *colMap_ , *Map_ ); //FIXME wrong

  // Note: B() = *ExtMatrix_ .

  Importer_ = rcp( new Epetra_Import(*Map_,A().RowMatrixRowMap()) ); //FIXME is this right?!

  // fix indices for overlapping matrix
  NumMyRowsB_ = B().NumMyRows();
  NumMyRows_ = NumMyRowsA_ + NumMyRowsB_;  //TODO is this wrong for a subdomain on >1 processor? // should be ok
  //NumMyCols_ = NumMyRows_;  //TODO is this wrong for a subdomain on >1 processor?  // YES!!!
  //NumMyCols_ = A().NumMyCols() + B().NumMyCols();
  NumMyCols_ = B().NumMyCols();

  /*FIXME*/ //somehow B's NumMyCols is the entire subdomain (local + overlap)

  NumMyDiagonals_ = A().NumMyDiagonals() + B().NumMyDiagonals();
  
  NumMyNonzeros_ = A().NumMyNonzeros() + B().NumMyNonzeros();
  Comm().SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  MaxNumEntries_ = A().MaxNumEntries();
  
  if (MaxNumEntries_ < B().MaxNumEntries())
    MaxNumEntries_ = B().MaxNumEntries();

/*
  if (!Comm().MyPID()) printf("========\nA matrix\n========\n"); fflush(stdout);
  EpetraExt::OperatorToMatlabFile("Aorig.dat",A());
  sleep(3);
  Comm().Barrier();
*/

/*
  if (!Comm().MyPID()) printf("========\nB matrix\n========\n"); fflush(stdout);
  cout << *ExtMatrix_ << endl;
  //EpetraExt::OperatorToMatlabFile("Bmatrix.dat",B());
  sleep(3);
  Comm().Barrier();

  if (Comm().MyPID()==0)
    printf("=======================================\nIfpack_OverlappingRowMatrix Importer_\n=======================================\n"); fflush(stdout);
  cout << *(A().RowMatrixImporter()) << endl;

  sleep(2);
*/

# ifdef HAVE_MPI
  delete nodeComm;
# endif

} //Ifpack_OverlappingRowMatrix() ctor for more than one core
#endif //ifdef IFPACK_NODE_AWARE_CODE

// ======================================================================
int Ifpack_OverlappingRowMatrix::
NumMyRowEntries(int MyRow, int & NumEntries) const
{
  if (MyRow < NumMyRowsA_)
    return(A().NumMyRowEntries(MyRow,NumEntries));
  else
    return(B().NumMyRowEntries(MyRow - NumMyRowsA_, NumEntries));
}

#ifdef IFPACK_NODE_AWARE_CODE
// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractMyRowCopy(int LocRow, int Length, int & NumEntries, double *Values, 
                 int * Indices) const
{
  assert(1==0);
  int ierr;
  const Epetra_Map *Themap;
  if (LocRow < NumMyRowsA_) {
    ierr = A().ExtractMyRowCopy(LocRow,Length,NumEntries,Values,Indices);
    Themap=&A().RowMatrixColMap();
  }
  else {
    ierr = B().ExtractMyRowCopy(LocRow-NumMyRowsA_,Length,NumEntries,Values,Indices);
    Themap=&B().RowMatrixColMap();
  }

  IFPACK_RETURN(ierr);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractGlobalRowCopy(int GlobRow, int Length, int & NumEntries, double *Values, 
                 int * Indices) const
{
  int ierr;
  const Epetra_Map *Themap;
  int LocRow = A().RowMatrixRowMap().LID(GlobRow);
  if (LocRow < NumMyRowsA_ && LocRow != -1) { //TODO don't need to check less than nummyrows
    ierr = A().ExtractMyRowCopy(LocRow,Length,NumEntries,Values,Indices);
    Themap=&A().RowMatrixColMap();
  }
  else {
    LocRow = B().RowMatrixRowMap().LID(GlobRow);
    assert(LocRow!=-1);
    //ierr = B().ExtractMyRowCopy(LocRow-NumMyRowsA_,Length,NumEntries,Values,Indices);
    ierr = B().ExtractMyRowCopy(LocRow,Length,NumEntries,Values,Indices);
    Themap=&B().RowMatrixColMap();
  }

  for (int i=0; i<NumEntries; i++) {
    Indices[i]=Themap->GID(Indices[i]);
    assert(Indices[i] != -1);
  }

  IFPACK_RETURN(ierr);
}
#else

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values,
                 int * Indices) const
{
  int ierr;
  if (MyRow < NumMyRowsA_)
    ierr = A().ExtractMyRowCopy(MyRow,Length,NumEntries,Values,Indices);
  else
    ierr = B().ExtractMyRowCopy(MyRow - NumMyRowsA_,Length,NumEntries,
                                Values,Indices);

  IFPACK_RETURN(ierr);
}
#endif

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  IFPACK_CHK_ERR(-1);
}


// ======================================================================
int Ifpack_OverlappingRowMatrix::
Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();
  vector<int> Ind(MaxNumEntries_);
  vector<double> Val(MaxNumEntries_);

  Y.PutScalar(0.0);

  // matvec with A (local rows)
  for (int i = 0 ; i < NumMyRowsA_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(A().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i] += Val[j] * X[k][Ind[j]];
      }
    }
  }

  // matvec with B (overlapping rows)
  for (int i = 0 ; i < NumMyRowsB_ ; ++i) {
    for (int k = 0 ; k < NumVectors ; ++k) {
      int Nnz;
      IFPACK_CHK_ERR(B().ExtractMyRowCopy(i,MaxNumEntries_,Nnz, 
                                          &Val[0], &Ind[0]));
      for (int j = 0 ; j < Nnz ; ++j) {
        Y[k][i + NumMyRowsA_] += Val[j] * X[k][Ind[j]];
      }
    }
  }
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Multiply(UseTranspose(),X,Y));
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-1);
}

// ======================================================================
#ifndef IFPACK_NODE_AWARE_CODE
Epetra_RowMatrix& Ifpack_OverlappingRowMatrix::B() const
{
  return(*ExtMatrix_);
}
#endif
// ======================================================================
const Epetra_BlockMap& Ifpack_OverlappingRowMatrix::Map() const
{
  return(*Map_);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ImportMultiVector(const Epetra_MultiVector& X, Epetra_MultiVector& OvX,
                  Epetra_CombineMode CM)
{
  OvX.Import(X,*Importer_,CM);
  return(0);
}

// ======================================================================
int Ifpack_OverlappingRowMatrix::
ExportMultiVector(const Epetra_MultiVector& OvX, Epetra_MultiVector& X,
                  Epetra_CombineMode CM)
{
  X.Export(OvX,*Importer_,CM);
  return(0);
}

