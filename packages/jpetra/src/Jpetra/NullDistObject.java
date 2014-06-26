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

/**
 *
 * @author  Jason Cross
 */
public class NullDistObject extends DistObject {
    private VectorSpace vectorSpace;
    
    /** Creates a new instance of NullDistObject */
    public NullDistObject(VectorSpace vectorSpace) {
        this.vectorSpace = vectorSpace;
    }
    
    public void copyAndPermute(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids, int combineMode) {
        // empty
    }
    
    public Serializable[] packAndPrepare(DistObject distObjectSource, int[] exportGids, int[] exportLids) {
        Serializable[] exportData = new Serializable[exportGids.length];  // the object to be exported by Distributor.distribute
        
        // we just need to packup the gids
        for(int i=0; i < exportGids.length; i++) {
            exportData[i] = new Integer(exportGids[i]);
        }
        
        return exportData;
    }
    
    public int[][] unpackAndCombine(Serializable[] importData, int combineMode) {
        Serializable[] entry;
        Serializable[] entryData;
        int gid;
        int lid;
        
        // for reverse op
        int sumEntries = 0;
        for(int i=0; i < importData.length; i++) {
            // if a vnode didn't send us any data, then importData[vnodeId] == null
            if (importData[i] != null) {
                sumEntries += ((Serializable[]) importData[i]).length;
            }
        }
        int[][] reverseExportVnodeIdsGidsLids = new int[4][];
        reverseExportVnodeIdsGidsLids[0] = new int[vectorSpace.getComm().getNumVnodes()];
        reverseExportVnodeIdsGidsLids[1] = new int[sumEntries];
        reverseExportVnodeIdsGidsLids[2] = new int[sumEntries];
        reverseExportVnodeIdsGidsLids[3] = new int[sumEntries];
        int revCount = 0;
        for(int i=0; i < importData.length; i++) {
            // if a vnode didn't send us any data, then importData[vnodeId] == null
            if (importData[i] != null) {
                entry = (Serializable[]) importData[i];  // get the array of elements send to us by the vnode i
                // unpack each element
                for(int j=0; j < entry.length; j++) {
                    gid = ((Integer) entry[j]).intValue();
                    lid = vectorSpace.getLocalIndex(gid);
                    reverseExportVnodeIdsGidsLids[0][i]++;
                    reverseExportVnodeIdsGidsLids[1][revCount] = i;
                    reverseExportVnodeIdsGidsLids[2][revCount] = gid;
                    reverseExportVnodeIdsGidsLids[3][revCount++] = lid;
                }
            }
        }
        
        return reverseExportVnodeIdsGidsLids;
    }
    
    public VectorSpace getVectorSpace() {
        return this.vectorSpace;
    }
}
