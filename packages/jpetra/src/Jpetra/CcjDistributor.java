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

import java.io.Serializable;
import java.util.LinkedHashSet;

/**
 *
 * @author  Jason Cross
 */
public class CcjDistributor extends JpetraObject implements Distributor {
    int numReceives;
    Comm comm;
    boolean vnodesInOrder;
    int[] startIndices;
    int[] numSends;
    int[] nextIndex;
    int[] exportVnodeIds;
    int[] senders;
    
    int[] packedGidsToSend; //used for createFromReceives
    
    
    // for reverse op.
    int[] reverseSenders;
    private boolean doneForwardOp;
    private int[][] reverseExportVnodeIdsGidsLids;
    
    public CcjDistributor() {
        //this.outputStreams.put("DISTRIBUTOR", new Output("CcjDistributor: ", true, System.out, false, System.out));
        // empty
    }
    
    public int[] createFromReceives(int[] remoteGids, int[] remoteVnodeIds, Comm comm) {
        this.comm = comm;
        
        createFromSends(computeSends(remoteGids, remoteVnodeIds), this.comm);
        return this.packedGidsToSend;
    }
    
    public int[] computeSends(int[] remoteGids, int[] remoteVnodeIds) {
        Distributor distributor = comm.createDistributor();
        distributor.createFromSends(remoteVnodeIds, comm);
        int[][] gidsToSend = distributor.distribute(remoteGids);
        
        int totalGidsToSend = 0;
        for(int i=0; i < gidsToSend.length; i++) {
            if (gidsToSend[i] == null) {
                gidsToSend[i] = new int[0];
            }
            totalGidsToSend += gidsToSend[i].length;
        }
        int[] sendToVnodeIds = new int[totalGidsToSend];
        this.packedGidsToSend = new int[totalGidsToSend];
        int index = 0;
        for(int i=0; i < gidsToSend.length; i++) {
            for(int j=0; j < gidsToSend[i].length; j++) {
                sendToVnodeIds[index] = i;
                this.packedGidsToSend[index++] = gidsToSend[i][j];
            }
        }
        
        return sendToVnodeIds;
    }
    
    public void createFromSends(int[] exportVnodeIds, Comm comm) {
        this.exportVnodeIds = exportVnodeIds;
        this.comm = comm;
        
        // figure out if the exportVnodeIds are linearly blocked together by vnodeId or if the vnodeIds are out of order
        boolean vnodesInOrder = true;
        this.startIndices = new int[comm.getNumVnodes()];
        this.numSends = new int[comm.getNumVnodes()];
        
        if (exportVnodeIds.length != 0) {
            int iTemp;
            for(iTemp=0; iTemp < this.exportVnodeIds.length-1; iTemp++) {
                if (this.vnodesInOrder && (this.exportVnodeIds[iTemp] > this.exportVnodeIds[iTemp+1])) {
                    this.vnodesInOrder = false;
                }
                this.numSends[this.exportVnodeIds[iTemp]]++;  // count up how many elements to send to the vnode
            }
            this.numSends[this.exportVnodeIds[iTemp]]++;
            
            if (this.vnodesInOrder) {
                // easy setup since the exportVnodeIds are in linear order...
                for(int i=0; i < this.numSends.length-1; i++) {
                    this.startIndices[i+1] = this.numSends[i] + this.startIndices[i];
                }
            }
            else {
                // vnodes are out of order
                this.nextIndex = new int[exportVnodeIds.length];
                int[] previousIndex = new int[comm.getNumVnodes()];
                int[] count = new int[comm.getNumVnodes()];  // counts how many gid will be sent to the respective vnode
                for(int i=0; i < this.startIndices.length; i++) {
                    this.startIndices[i] = -1;
                    previousIndex[i] = -1;
                }
                
                int lastIndex;
                for(int i=0; i < this.exportVnodeIds.length; i++) {
                    lastIndex = previousIndex[this.exportVnodeIds[i]];
                    if (lastIndex == -1) {
                        this.startIndices[this.exportVnodeIds[i]] = i;
                    }
                    else {
                        this.nextIndex[lastIndex] = i;
                    }
                    previousIndex[exportVnodeIds[i]] = i;
                    count[this.exportVnodeIds[i]]++;
                }
            }
        }
        
        ComputeReceives();
    }
    
    public void ComputeReceives() {
        int[] receivingVnodes = new int[comm.getNumVnodes()];  //remember that Java 0's all elements in new int arrays
        for(int i=0; i < receivingVnodes.length; i++) {
            if (numSends[i] > 0) {
                receivingVnodes[i] = 1;
            }
        }
        
        // collect all the sending to receiving vnode mappings on the root node
        int[][] allReceivingVnodes = comm.gather(receivingVnodes);
        // transpose the receives to correspond to each receiving vnode
        int[][] allReceivingTranspose = new int[comm.getNumVnodes()][comm.getNumVnodes()];
        if (comm.getVnodeId() == 0) {
            for(int i=0; i < allReceivingVnodes.length; i++) {
                for(int j=0; j < allReceivingVnodes[i].length; j++) {
                    allReceivingTranspose[j][i] = allReceivingVnodes[i][j];
                }
            }
        }
        // scatter the corresponding receiving arrays to each vnode from the root node
        this.senders = comm.scatter2dArray(allReceivingTranspose);
        // setup for receiving is done
    }
    
    public int[][] distribute(int[] toSendData) {
        // do async_sends
        int[] buffer;
        int dataIndex;
        for(int i=0; i < numSends.length; i++) {
            if(numSends[i] > 0) {
                //this.println("DISTRIBUTOR", "Sending " + numSends[i] + " objects to vnode " + i);
                // we're going to send data to vnode i
                // so buffer up all send objects
                buffer = new int[numSends[i]];
                dataIndex = this.startIndices[i];
                for(int j=0; j < numSends[i]; j++) {
                    buffer[j] = toSendData[dataIndex];
                    dataIndex = nextIndex[dataIndex];
                }
                // buffer object filled so send it off to vnode i
                comm.send(buffer, i);
            }
        }
        // done doing sends, now do receives from known list of senders
        int[][] receivedData = new int[senders.length][];
        for(int i=0; i < senders.length; i++) {
            if (senders[i] == 1) {
                //this.println("DISTRIBUTOR", "Receiving from vnode " + i);
                receivedData[i] = (int[]) comm.receive(i);
            }
        }
        // done doing receives, that means we're all finished
        
        return receivedData;
    }
    
    public Serializable[] distribute(Serializable[] exportObjects, boolean doReverse) {
        // do async_sends
        Serializable[] buffer;
        int dataIndex;
        if (!doReverse) {
            this.reverseSenders = new int[comm.getNumVnodes()];
            for(int i=0; i < numSends.length; i++) {
                if(numSends[i] > 0) {
                    this.reverseSenders[i] = 1;  // if we reverse this op. then we will receive from vnode i
                    //this.println("STD", "Sending " + numSends[i] + " objects to vnode " + i);
                    
                /*this.println("STD", "startIndices.length: " + startIndices.length);
                for(int i2=0; i2 < this.startIndices.length; i2++) {
                    this.println("STD", "startIndices[" + i2 + "]=" + this.startIndices[i2]);
                }*/
                    
                    // we're going to send data to vnode i
                    // so buffer up all send objects
                    buffer = new Serializable[numSends[i]];
                    dataIndex = this.startIndices[i];
                    for(int j=0; j < numSends[i]; j++) {
                        buffer[j] = exportObjects[dataIndex];
                        dataIndex = nextIndex[dataIndex];
                    }
                    // buffer object filled so send it off to vnode i
                    comm.send(buffer, i);
                }
            }
        }
        else {
            dataIndex = 0;
            for(int i=0; i < reverseExportVnodeIdsGidsLids[0].length; i++) {
                if (reverseExportVnodeIdsGidsLids[0][i] > 0) {
                    //this.println("DISTRIBUTOR", "Sending " + reverseExportVnodeIdsGidsLids[0][i] + " objects to vnode " + i);
                    
                    // we're going to send data to vnode i
                    // so buffer up all send objects
                    buffer = new Serializable[reverseExportVnodeIdsGidsLids[0][i]];
                    for(int j=0; j < reverseExportVnodeIdsGidsLids[0][i]; j++) {
                        buffer[j] = exportObjects[dataIndex++];
                    }
                    // buffer object filled so send it off to vnode i
                    comm.send(buffer, i);
                }
            }
        }
        
        // done doing sends, now do receives from known list of senders
        Serializable[] receivedData = new Serializable[comm.getNumVnodes()];
        if (!doReverse) {
            for(int i=0; i < senders.length; i++) {
                if (senders[i] == 1) {
                    //("DISTRIBUTOR", "Receiving from vnode " + i);
                    receivedData[i] = comm.receive(i);
                }
            }
        } else {
            for(int i=0; i < reverseSenders.length; i++) {
                if (reverseSenders[i] == 1) {
                    //this.println("DISTRIBUTOR", "Receiving from vnode " + i);
                    receivedData[i] = comm.receive(i);
                }
            }
        }
        // done doing receives, that means we're all finished
        
        return receivedData;
    }
    
    public int[] getSenders() {
        return this.senders;
    }
    
    public int[] getExportVnodeIds() {
        return this.exportVnodeIds;
    }
    
    public void setReverseExportVnodeIdsGidsLids(int[][] reverseExportVnodeIdsGidsLids) {
        this.reverseExportVnodeIdsGidsLids = reverseExportVnodeIdsGidsLids;
    }
    
    public int[] getReverseExportVnodeIds() {
        return this.reverseExportVnodeIdsGidsLids[1];
    }
    
    public int[] getReverseExportGids() {
        return this.reverseExportVnodeIdsGidsLids[2];
    }
    
    public int[] getReverseExportLids() {
        return this.reverseExportVnodeIdsGidsLids[3];
    }
    
    public boolean doneForwardOp() {
        return this.doneForwardOp;
    }
    
    public void setDoneForwardOp(boolean doneForwardOp) {
        this.doneForwardOp = doneForwardOp;
    }
}
