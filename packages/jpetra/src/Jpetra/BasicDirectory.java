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

/**
 *
 * @author  Jason Cross
 */
public class BasicDirectory extends JpetraObject implements Directory {
    private VectorSpace vectorSpace;
    private int numVnodes;
    private int[] minVnodeGids;
    private int minMyGid;
    
    public BasicDirectory(VectorSpace vectorSpace) {
        this.vectorSpace = vectorSpace;
        
        if (vectorSpace.isDistributedGlobally()) {
            return;  // nothing to setup
        }
        
        if (vectorSpace.isDistributedLinearly()) {
            numVnodes = vectorSpace.getComm().getNumVnodes();
            minVnodeGids = new int[numVnodes+1];
            minMyGid = vectorSpace.getMinGlobalEntryId();
            vectorSpace.getComm().gatherAll(vectorSpace.getMinGlobalEntryId());
            minVnodeGids[numVnodes] = 1 + vectorSpace.getMaxGlobalEntryId(); // Set max cap
        }
        
        this.println("ERR", "The vector space is neither distributed globally or linearly, thus it is not supported.");
    }
    
    public int[][] getDirectoryEntries(int[] globalElements) {
        return null; //temporary so class will compile
    }
    
}
