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
    private VectorSpace directoryVectorSpace;
    private int[] directoryVnodeIds;
    private int[] directoryLids;
    
    public BasicDirectory(VectorSpace vectorSpace) {
        this.vectorSpace = vectorSpace;
        
        if (!vectorSpace.isDistributedGlobally()) {
            return;  // nothing to setup
        }
        
        if (vectorSpace.isDistributedLinearly()) {
            return; // nothing to setup
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
        this.println("STD", "numGlobalEntries: " + directoryVectorSpace.getNumGlobalEntries() + " myMin: " + this.directoryVectorSpace.getMyMinGlobalIndex() + " myMax: " + this.directoryVectorSpace.getMyMaxGlobalIndex());
        this.println("STD", "minGlobalEntryId: " + directoryVectorSpace.getMinGlobalEntryId() + " maxGlobalEntryId: " + directoryVectorSpace.getMaxGlobalEntryId());
        int[] myGidsTemp = directoryVectorSpace.getMyGlobalEntryIds();
        for(int i=0; i < myGidsTemp.length; i++) {
            this.println("STD", "Gid: " + myGidsTemp[i]);
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
            this.println("STD", "packing my gid: " + myGids[i]);
            toSendData[i] = new int[]{myGids[i], vectorSpace.getLocalIndex(myGids[i])};
        }
        // data is packed, so send it off
        Serializable[] receivedData = distributor.distribute(toSendData);
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
                        this.println("STD", "Adding gid location " + gidLid[0] + " to my directory.");
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
            this.println("STD", "LID: " + i + " GID: " + currentGid++ + " Owner: " + directoryVnodeIds[i] + " OwnerLid: " + directoryLids[i]);
        }
        this.println("STD", "BasicDirectory.generate has finished.");
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
        int[][] tmp = this.directoryVectorSpace.getRemoteVnodeIdList(globalElements);
        Distributor distributor = vectorSpace.getComm().createDistributor();
        distributor.createFromReceives(globalElements, tmp[0], vectorSpace.getComm());  // tmp[0] is the directorVnodes array

    }  
}
