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

import java.io.Serializable;
import java.util.LinkedHashSet;

/**
 *
 * @author  Jason Cross
 */
public class CcjDistributor extends JpetraObject implements Distributor {
    int numReceives;
    public CcjDistributor() {
    }
    
    /*public int[] createFromRecieves(int[] remoteGlobalElementIds, int[] remoteGlobalVnodeIds, int[] exportElementIds, int[] exportVnodeIds) {
        return null;  // !! not implemented
    }*/
    public void createFromReceives(int[] importVnodeIds) {
        //TreeSet treeSet = new TreeSet();
        this.numReceives = 0;
        int current = -1;
        for(int i=0; i < importVnodeIds.length; i++) {
            if (current != importVnodeIds[i]) {
                current = importVnodeIds[i];
                this.numReceives++;
                //treeSet.add(new Integer(importVnodeIds[i]));
            }
        }
        
    }
    
    public void createFromSends(int[] exportVnodeIds, Comm comm) {
        // figure out if the exportVnodeIds are linearly blocked together by vnodeId or if the vnodeIds are out of order
        boolean vnodesInOrder = true;
        int[] startIndices = new int[comm.getNumVnodes()];
        int[] numSends = new int[comm.getNumVnodes()];
        for(int i=0; i < exportVnodeIds.length-1; i++) {
            if (vnodesInOrder && (exportVnodeIds[i] > exportVnodeIds[i+1])) {
                vnodesInOrder = false;
            }
            numSends[exportVnodeIds[i]]++;  // count up how many elements to send to the vnode
        }
        
        if (vnodesInOrder) {
            // easy setup since the exportVnodeIds are in linear order...
            for(int i=0; i < numSends.length-1; i++) {
                startIndices[i+1] = numSends[i] + startIndices[i];
            }
        }
        else {
            // vnodes are out of order
            int[] nextIndex = new int[exportVnodeIds.length];
            int[] previousIndex = new int[comm.getNumVnodes()];
            int[] count = new int[comm.getNumVnodes()];
            for(int i=0; i < startIndices.length; i++) {
                startIndices[i] = -1;
                previousIndex[i] = -1;
            }
            
            int lastIndex;
            for(int i=0; i < exportVnodeIds.length; i++) {
                lastIndex = previousIndex[exportVnodeIds[i]];
                if (lastIndex == -1) {
                    startIndices[exportVnodeIds[i]] = i;
                }
                else {
                    nextIndex[lastIndex] = i;
                }
                previousIndex[exportVnodeIds[i]] = i;
                count[i]++;
            }
        }
    }
    
    public void distribute(Comm comm, int[] exportVnodeIds, int[] exportGids, int[] exportLids, Serializable[] exportObjects) {
        Serializable[] importObjects = new Serializable[this.numReceives];
        
        // do async sends
        int current = exportVnodeIds[0];
        int start = 0;
        int numSend = 0;
        for(int i=0; i < exportVnodeIds.length; i++) {
            if (current != exportVnodeIds[i]) {
                // copy gids and objects from start to start+numSend
                // do async_send to vnode current
                this.println("STD", "sending to " + current);
                // setup next array slice
                start = i;
                numSend = 1;
                current = exportVnodeIds[i];
            }
            numSend++;
        }
        // do blocking receives
        for(int i=0; i < this.numReceives; i++) {
            // do receive
            this.println("STD", "Waiting to receive... " + 1);
        }
    }
}
