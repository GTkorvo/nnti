// @HEADER
// ***********************************************************************
//
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

package Jpetra;

import java.util.ArrayList;
import java.util.Iterator;
import java.io.Serializable;

/**
 * <code>BasicDirectory</code> is not used directly by the user.
 *
 * @author Jason Cross
 */
public class BasicDirectory extends JpetraObject implements Directory {
    private VectorSpace vectorSpace;
    private VectorSpace directoryVectorSpace;
    private int[] directoryVnodeIds;
    private int[] directoryLids;
    private int[] allMinGids;  // for linear nonuniform distributions
    private int[] numGidsOnVnodes; // for linear nonuniform distributions
    
    public BasicDirectory(VectorSpace vectorSpace) {
        //this.outputStreams.put("DIRECTORY", new Output("BasicDirectory: ", true, System.out, false, System.out));
        
        this.vectorSpace = vectorSpace;
        
        if (!vectorSpace.isDistributedGlobally()) {
            return;  // nothing to setup
        }
        
        if (vectorSpace.isDistributedUniformly()) {
            return; // nothing to setup
        }
        
        if (vectorSpace.isDistributedLinearly()) {
            int[] tmp = new int[]{vectorSpace.getMyMinGlobalIndex(), vectorSpace.getNumMyGlobalEntries()};
            //this.println("DIRECTORY", "Submitting min gid: " + tmp[0] + "to the gatherAll()");
            tmp = vectorSpace.getComm().gatherAll(tmp);
            this.allMinGids = new int[tmp.length + 1];
            System.arraycopy(tmp, 0, this.allMinGids, 0, tmp.length);
            this.allMinGids[allMinGids.length - 1] = vectorSpace.getMaxGlobalEntryId(); // set max cap for GID searching algorithm below in getDirectoryEntries
            //this.println("DIRECTORY", "vectorSpace.getMyMinGlobalIndex(): " + vectorSpace.getMyMinGlobalIndex());
            for(int i=0; i < this.allMinGids.length-1; i+=2) {
                //this.println("DIRECTORY", "vnode: " + i/2 + " minGid: " + this.allMinGids[i] + " numGids: " + this.allMinGids[i+1]);
            }
            
            return; // nothing else left to do
        }
        
        // setup arbitrary distribution
        generateDirectoryTable();
    }
    
    private void generateDirectoryTable() {
        int minGlobalGid = vectorSpace.getMinGlobalEntryId();
        int maxGlobalGid = vectorSpace.getMaxGlobalEntryId();
        int numDirectoryGlobalEntries = maxGlobalGid - minGlobalGid + 1;
        this.directoryVectorSpace = new VectorSpace(new ElementSpace(numDirectoryGlobalEntries, minGlobalGid, vectorSpace.getComm()));
        
        // debug stuff
        //this.println("DIRECTORY", "numGlobalEntries: " + directoryVectorSpace.getNumGlobalEntries() + " myMin: " + this.directoryVectorSpace.getMyMinGlobalIndex() + " myMax: " + this.directoryVectorSpace.getMyMaxGlobalIndex());
        //this.println("DIRECTORY", "minGlobalEntryId: " + directoryVectorSpace.getMinGlobalEntryId() + " maxGlobalEntryId: " + directoryVectorSpace.getMaxGlobalEntryId());
        int[] myGidsTemp = directoryVectorSpace.getMyGlobalEntryIds();
        for(int i=0; i < myGidsTemp.length; i++) {
            //this.println("DIRECTORY", "Gid: " + myGidsTemp[i]);
        }
        //end debug stuff
        
        int numMyGlobalDirectoryEntries = directoryVectorSpace.getNumMyGlobalEntries();
        
        // Get list of processors owning the directory entries for the Map GIDs
        int[] myGids = vectorSpace.getMyGlobalEntryIds();
        int[][] tmp = directoryVectorSpace.getRemoteVnodeIdList(myGids);  // get remote vnodeIds
        int[] sendGidsToVnodes = tmp[0];
        
        // use distributor to send out my Gids to those vnodes who own them in the directory
        Distributor distributor = vectorSpace.getComm().createDistributor();
        distributor.createFromSends(sendGidsToVnodes, vectorSpace.getComm());
        // now pack up the gids/lids that we have to send to the owner gid directory vnodes
        Serializable[] toSendData = new Serializable[sendGidsToVnodes.length];
        for(int i=0; i < myGids.length; i++) {
            //this.println("DIRECTORY", "packing my gid: " + myGids[i]);
            toSendData[i] = new int[]{myGids[i], vectorSpace.getLocalIndex(myGids[i])};
        }
        // data is packed, so send it off
        Serializable[] receivedData = distributor.distribute(toSendData, false);
        int[] senders = distributor.getSenders();
        Serializable[] gidsLids;
        int directoryLid;
        int[] gidLid;
        // now unpack data into the tables that support directory lookup
        this.directoryVnodeIds = new int[directoryVectorSpace.getNumMyGlobalEntries()];
        this.directoryLids = new int[directoryVectorSpace.getNumMyGlobalEntries()];
        for(int i=0; i < senders.length; i++) {
            if (senders[i] == 1) {
                gidsLids = (Serializable[]) receivedData[i];
                for(int j=0; j < gidsLids.length; j++) {
                    gidLid = (int[]) gidsLids[j];
                    directoryLid = directoryVectorSpace.getLocalIndex(gidLid[0]);
                    if (directoryLid == -1) {
                        this.println("ERR", "I don't keep track of gid " + gidLid[0] + " in my directory.");
                    } else {
                        //this.println("DIRECTORY", "Adding gid location " + gidLid[0] + " to my directory.");
                        this.directoryVnodeIds[directoryLid] = i;
                        this.directoryLids[directoryLid] = gidLid[1];
                    }
                }
            }
        }
        // the directory tables are now built
        
        //debug code
        int currentGid = directoryVectorSpace.getMyMinGlobalIndex();
        for(int i=0; i < directoryVectorSpace.getNumMyGlobalEntries(); i++) {
            //this.println("DIRECTORY", "LID: " + i + " GID: " + currentGid++ + " Owner: " + directoryVnodeIds[i] + " OwnerLid: " + directoryLids[i]);
        }
        //this.println("DIRECTORY", "BasicDirectory.generate has finished.");
    }
    
    /**
     *
     * @return vnodeIds are in int[0] and Lids in int[1]
     */
    public int[][] getDirectoryEntries(int[] globalElements) {
        int[][] vnodeIdsLids = new int[2][globalElements.length];
        
        // for serial
        if (!vectorSpace.isDistributedGlobally()) {
            // since Java initializes all elements to zero, and the vnodeId=0, we can leave vnodeIdsLids[0] alone
            for(int i=0; i < vnodeIdsLids[1].length; i++) {
                vnodeIdsLids[1][i] = vectorSpace.getLocalIndex(globalElements[i]);
            }
            return vnodeIdsLids;
        }
        
        // for parallel
        // for linear continous distribution
        if (vectorSpace.isDistributedUniformly()) {
            int numRemainderIndicies = vectorSpace.getNumRemainderIndices();
            int numIndiciesPerVnode = vectorSpace.getNumIndicesPerVnode();
            int minGlobalElementId = vectorSpace.getMinGlobalEntryId();
            int maxGlobalElementId = vectorSpace.getMaxGlobalEntryId();
            int numGlobalElements = vectorSpace.getNumGlobalEntries();
            int maxRemainderIndex = numRemainderIndicies + (numIndiciesPerVnode * numRemainderIndicies);
            int gid;
            for(int i=0; i < globalElements.length; i++) {
                gid = globalElements[i];
                if (gid >= maxRemainderIndex) {
                    vnodeIdsLids[0][i] = ((gid - maxRemainderIndex) / numIndiciesPerVnode) + numRemainderIndicies;
                    vnodeIdsLids[1][i] = ((gid - maxRemainderIndex) % numIndiciesPerVnode);
                }
                else {
                    vnodeIdsLids[0][i] = gid / (numIndiciesPerVnode + 1);
                    vnodeIdsLids[1][i] = gid % (numIndiciesPerVnode + 1);
                    
                }
            }
            
            return vnodeIdsLids;
        }
        
        if (vectorSpace.isDistributedLinearly()) {
            //this.println("DIRECTORY", "trying to find Gids using Linear Distribution...");
            int minAllGids = this.vectorSpace.getMinGlobalEntryId();
            int maxAllGids = this.vectorSpace.getMaxGlobalEntryId();
            int numVnodes = this.vectorSpace.getComm().getNumVnodes();
            int n_over_p = numVnodes / this.vectorSpace.getNumGlobalEntries();
            
            for (int i=0; i < globalElements.length; i++) {
                int LID = -1; // Assume not found
                int Proc = -1;
                int GID = globalElements[i];
                boolean found;
                boolean wrap;
                boolean doneWrap;
                int Proc1;
                //this.println("DIRECTORY", "Checking for gid: " + GID);
                if (GID < minAllGids) {
                    //this.println("DIRECTORY", "gid: " + GID + " is less than minAllGids: " + minAllGids);
                    // error
                }
                else if (GID > maxAllGids) {
                    //this.println("DIRECTORY", "gid: " + GID + " is greater than maxAllGids: " + maxAllGids);
                    // error
                }
                else {
                    // Guess uniform distribution and start a little above it
                    Proc1 = Util.min(GID/Util.max(n_over_p,1) + 2, numVnodes-1);
                    Proc1 = Proc1 * 2;
                    found = false;
                    wrap = true;
                    doneWrap = false;
                    while(wrap) {
                        while (Proc1 >= 0 && Proc1 < numVnodes*2) {
                            ////this.println("DIRECTORY", "Proc1: " + Proc1/2 + " GID: " + GID);
                            if (allMinGids[Proc1] <= GID) {
                                //if (GID < allMinGids[Proc1+2]) {
                                if ((allMinGids[Proc1 + 1] + allMinGids[Proc1]) > GID) { // allMinGids[Proc1] is the minGid of that vnode, allMinGids[Proc1 + 1] is the number of Gids on that vnode
                                    found = true;
                                    break;
                                }
                                else Proc1+=2;
                            }
                            else Proc1-=2;
                        }
                        
                        if (found) {
                            wrap = false;
                        }
                        else if(allMinGids[Proc1-2] != allMinGids[Proc1]) {
                            wrap = false;
                        }
                        else {
                            if (doneWrap) {
                                wrap = false;
                            }
                            else {
                                doneWrap = true;
                                Proc1 = 0;
                            }
                        }
                    }
                    if (found) {
                        Proc = Proc1;
                        LID = GID - allMinGids[Proc];
                        //this.println("DIRECTORY", "gid: " + GID + "found on Proc: " + Proc);
                    }
                }
                vnodeIdsLids[0][i] = Proc;
                vnodeIdsLids[1][i] = LID;
            }
            
            return vnodeIdsLids;
        }
        
        /*//this.println("STD", "global elements this vnode needs:");
        for(int i=0; i < globalElements.length; i++) {
            //this.println("STD", globalElements[i] + "");
        }*/
        
        // for arbitrary distribution
        int[][] tmp = this.directoryVectorSpace.getRemoteVnodeIdList(globalElements);
        Distributor distributor = vectorSpace.getComm().createDistributor();
        int[] gidsToSend = distributor.createFromReceives(globalElements, tmp[0], vectorSpace.getComm());  // tmp[0] is the directoryVnodes array
        Serializable[] packedGidsToSend = new Serializable[gidsToSend.length];
        int directoryLid;
        for(int i=0; i < gidsToSend.length; i++) {
            //this.println("DIRECTORY", "I need to send gid owner vnode/lid: " + gidsToSend[i]);
            directoryLid = directoryVectorSpace.getLocalIndex(gidsToSend[i]);
            packedGidsToSend[i] = new int[]{this.directoryVnodeIds[directoryLid], gidsToSend[i], this.directoryLids[directoryLid]};
        }
        
        Serializable[] receives = distributor.distribute(packedGidsToSend, false);
        int[] tmp1;
        Serializable[] tmp2;
        for(int i=0; i < receives.length; i++) {
            if (receives[i] != null) {
                tmp2 = (Serializable[]) receives[i];
                for(int j=0; j < tmp2.length; j++) {
                    tmp1 = (int[]) tmp2[j];
                    //this.println("DIRECTORY", "dirVnode: " + i + " ownwerVnode: " + tmp1[0] + " gid: " + tmp1[1] + " lid: " + tmp1[2]);
                }
            }
        }
        
        // insert the receives into the output array in the same order as the gids in the import array
        Serializable[] intArrays;
        int[] intArray;
        for(int i=0; i < receives.length; i++) {
            if (receives[i] != null) {
                intArrays = (Serializable[]) receives[i];
                for(int j=0; j < intArrays.length; j++) {
                    intArray = (int[]) intArrays[j];
                    for(int k=0; k < globalElements.length; k++) {
                        if (globalElements[k] == intArray[1]) {
                            vnodeIdsLids[0][k] = intArray[0];
                            vnodeIdsLids[1][k] = intArray[2];
                            break;
                        }
                    }
                }
            }
        }
        
        
        return vnodeIdsLids;
    }
}
