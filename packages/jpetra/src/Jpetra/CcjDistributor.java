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
//import java.util.TreeSet;

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
    
    public int createFromSends(int[] exportVnodeIds) {
        return 0;  // should change to public void
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
