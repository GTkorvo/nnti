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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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
 *
 * @author  Jason Cross
 */
public class BasicDirectory extends JpetraObject implements Directory {
    private VectorSpace vectorSpace;
    //private int numVnodes;
    //private int[] minVnodeGids;
    //private int minMyGid;
    
    public BasicDirectory(VectorSpace vectorSpace) {
        this.vectorSpace = vectorSpace;
        
        if (!vectorSpace.isDistributedGlobally()) {
            return;  // nothing to setup
        }
        
        if (vectorSpace.isDistributedLinearly()) {
            //numVnodes = vectorSpace.getComm().getNumVnodes();
            //minVnodeGids = new int[numVnodes+1];
            //minMyGid = vectorSpace.getMinGlobalEntryId();
            //vectorSpace.getComm().gatherAll(vectorSpace.getMinGlobalEntryId());
            //minVnodeGids[numVnodes] = 1 + vectorSpace.getMaxGlobalEntryId(); // Set max cap
            
            return; // nothing to setup
        }
        
        // setup arbitrary distribution
        generateDirectoryTable();
    }
    
    private void generateDirectoryTable() {
        int minGlobalGid = vectorSpace.getMinGlobalEntryId();
        int maxGlobalGid = vectorSpace.getMaxGlobalEntryId();
        int numDirectoryGlobalEntries = maxGlobalGid - minGlobalGid + 1;
        VectorSpace directoryVectorSpace = new VectorSpace(new ElementSpace(numDirectoryGlobalEntries, minGlobalGid, vectorSpace.getComm()));
        
        int numMyGlobalDirectoryEntries = directoryVectorSpace.getNumMyGlobalEntries();
        
        // Get list of processors owning the directory entries for the Map GIDs
        int[][] tmp = directoryVectorSpace.getRemoteVnodeIdList(vectorSpace.getMyGlobalEntryIds());  // get remote vnodeIds
        int[] sendGidsToVnodes = tmp[0];
        
        // use distributor to send out my Gids to those vnodes who own them in the directory
        Distributor distributor = vectorSpace.getComm().createDistributor();
        distributor.createFromSends(sendGidsToVnodes, vectorSpace.getComm());
        
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
        if (vectorSpace.isDistributedLinearly()) {
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
        
        /*this.println("STD", "global elements this vnode needs:");
        for(int i=0; i < globalElements.length; i++) {
            this.println("STD", globalElements[i] + "");
        }*/
        
        // for arbitrary distribution
        Comm comm = vectorSpace.getComm();
        int[][] neededGlobalElements = comm.gatherAll2dArray(globalElements);
        
        this.println("STD", "Gids to find...");
        for (int i=0; i < neededGlobalElements.length; i++) {
            this.println("STD", "vnode " + i + " needs Gids:");
            for(int j=0; j < neededGlobalElements[i].length; j++) {
                this.println("STD", neededGlobalElements[i][j] + "");
            }
        }
        // set index of vnodeId to 1 if I will send to that vnode
        int[] response = new int[neededGlobalElements.length];
        ArrayList[] responseGids = new ArrayList[neededGlobalElements.length];
        // go through neededGlobalElements, find those that belong to this vnode, and register the Gid to be sent to the neccessary vnode
        for (int i=0; i < neededGlobalElements.length; i++) {
            responseGids[i]= new ArrayList();
            for(int j=0; j < neededGlobalElements[i].length; j++) {
                if(vectorSpace.isMyGlobalIndex(neededGlobalElements[i][j])) {
                    response[i] = 1;
                    responseGids[i].add(new Integer(neededGlobalElements[i][j]));
                }
            }
        }
        
        // collect all the responses on the root node
        int[][] totalResponse = comm.gather(response);
        // transpose the responses to correspond to each receiving vnode
        int[][] totalResponseTranspose = new int[comm.getNumVnodes()][comm.getNumVnodes()];
        if (comm.getVnodeId() == 0) {
            for(int i=0; i < totalResponse.length; i++) {
                for(int j=0; j < totalResponse[i].length; j++) {
                    totalResponseTranspose[j][i] = totalResponse[i][j];
                }
            }
        }
        // scatter the corresponding response arrays to each vnode from the root node
        int[] senders = comm.scatter2dArray(totalResponseTranspose);
        // setup for receiving is done, now we delay receiving until after we have done all our own sends
        
        // do async_sends for all vnodes this vnode needs to send to
        this.println("STD", "Doing sends...");
        int[] individualResponse;  // the container for the Gids sent to other vnodes
        Iterator iterator;
        for(int i=0; i < responseGids.length; i++) {
            if (responseGids[i].size() > 0) {
                individualResponse = new int[responseGids[i].size()];  // output an int[] from the ArrayList
                individualResponse[0] = comm.getVnodeId();  // the first element will hold the vnode id of the sending vnode
                iterator = responseGids[i].iterator();
                for(int j=0; j < individualResponse.length; j++) {
                    individualResponse[j] = ((Integer) iterator.next()).intValue();
                }
                this.println("STD", "Sending to vnode " + i);
                comm.send(individualResponse, i);
            }
        }
        
        this.println("STD", "Doing barrier...");
        comm.barrier();  // this barrier may not be needed, for testing purposes
        
        // do blocking receives
        this.println("STD", "I'm vnode " + comm.getVnodeId() + "  Doing receives...");
        comm.barrier();  // this barrier is probably not needed, for testing purposes
        // find how many receives we have to do
        int numSenders = 0;
        for(int i=0; i < senders.length; i++) {
            if (senders[i] == 1) {
                numSenders++;
            }
        }
        // allocate space to receive the Gids owned by each respective vnode
        int[] senderVnodeIds = new int[numSenders];
        int[][] receivedGids = new int[numSenders][];
        int senderIndex = 0;
        // since each index in the senders array maps to a vnode, i is the vnode ID of the sending vnode
        for(int i=0; i < senders.length; i++) {
            if (senders[i] == 1) {
                senderVnodeIds[senderIndex] = i;
                receivedGids[senderIndex++] = (int[]) comm.receive(i);
            }
        }
        
        // unpack and process receives
        int[] gids;
        for(int i=0; i < receivedGids.length; i++) {
            this.println("STD", "Processing message...");
            gids = receivedGids[i];
            this.println("STD", "vnode " + senderVnodeIds[i] + " has Gids:");
            for(int j=0; j < gids.length; j++) {
                this.println("STD", gids[j] + "");
            }
        }
        
        this.println("STD", "Doing barrier...");
        comm.barrier();
        this.println("STD", "Done doing receives.");
        //this.println("ERR", "A non-serial non-linear continous vectorSpace was passed to BasicDirectory.getDirectoryEntries().  This is not supported yet!");
        return null; //temporary so class will compile
    }
    
}
