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
public class Import extends JpetraObject {
    VectorSpace sourceVectorSpace;
    VectorSpace targetVectorSpace;
    
    int numSameGids;
    int[] remoteLids;
    int[] remoteGids;
    int[] permuteToLids;
    int[] permuteFromLids;
    
    public Import(VectorSpace sourceVectorSpace, VectorSpace targetVectorSpace) {
        this.sourceVectorSpace = sourceVectorSpace;
        this.targetVectorSpace = targetVectorSpace;
        
        int[] sourceGids;
        if (sourceVectorSpace.getNumMyGlobalEntries() > 0) {
            sourceGids = sourceVectorSpace.getMyGlobalEntryIds();
        }
        else {
            sourceGids = new int[0];
        }
        
        int[] targetGids;
        if (targetVectorSpace.getNumMyGlobalEntries() > 0) {
            targetGids = targetVectorSpace.getMyGlobalEntryIds();
        }
        else {
            targetGids = new int[0];
        }
        
        //int numSameGids;
        int minNumIds = Util.min(sourceGids.length, targetGids.length);
        this.println("STD", "sourceGids.length: " + sourceGids.length + " targetGids.length: " + targetGids.length);
        for (numSameGids = 0; numSameGids < minNumIds; numSameGids++) {
            if (sourceGids[numSameGids] != targetGids[numSameGids]) {
                break;
            }
        }
        this.println("STD", "numSameGids: " + numSameGids);
        
        int numPermuteGids = 0;
        int numRemoteGids = 0;
        for (int i=numSameGids; i < targetGids.length; i++) {
            if (sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                this.println("STD", "Going to permute gid: " + targetGids[i]);
                numPermuteGids++;
            }
            else {
                numRemoteGids++;
            }
        }
        
        //int[] remoteLids;
        //int[] remoteGids;
        //int[] permuteToLids;
        //int[] permuteFromLids;
        if (numRemoteGids > 0) {
            remoteLids = new int[numRemoteGids];
            remoteGids = new int[numRemoteGids];
        }
        else {
            remoteLids = new int[0];
            remoteGids = new int[0];
        }
        if (numPermuteGids > 0)  {
            permuteToLids = new int[numPermuteGids];
            permuteFromLids = new int[numPermuteGids];
        }
        else {
            permuteToLids = new int[0];
            permuteFromLids = new int[0];
        }
        
        numPermuteGids = 0;
        numRemoteGids = 0;
        for (int i=numSameGids; i < targetGids.length; i++) {
            if (sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                this.println("STD", "targetGids[i]" + targetGids[i]);
                permuteFromLids[numPermuteGids] = i;
                permuteToLids[numPermuteGids++] = sourceVectorSpace.getLocalIndex(targetGids[i]);
            }
            else {
                remoteGids[numRemoteGids] = targetGids[i];
                remoteLids[numRemoteGids++] = i;
            }
        }
        
        if ((numRemoteGids > 0) && !sourceVectorSpace.isDistributedGlobally()) {
            this.println("ERR", "A non-distributed globally vector space has remote elements.");
        }
        
        this.println("STD", "targetGids.length: " + targetGids.length);
        this.println("STD", "permuteToLids.length: " + permuteToLids.length);
        
        
        /*int[] importGids = new int[numImports];
        int importGidsIndex = 0;
        this.println("STD", "importGids.length: " + importGids.length);
        for (int i=numSameGids; i < targetGids.length; i++) {
            if (!sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                importGids[importGidsIndex++] = targetGids[i];
            }
        }*/
        
        //int[] exportVnodeIds;
        int[] remoteVnodeIds;
        if (sourceVectorSpace.isDistributedGlobally()) {
            int[][] tmp = targetVectorSpace.getRemoteVnodeIdList(remoteGids);  // get remote vnodeIds
            remoteVnodeIds = tmp[0];
            
            //Get rid of IDs that don't exist in SourceMap
            if(remoteGids.length > 0) {
                int count = 0;
                for(int i = 0; i < remoteGids.length; i++)
                    if( remoteVnodeIds[i] == -1 ) count++;
                if( count > 0 ) {
                    if( remoteGids.length - count > 0 ) {
                        int[] newRemoteGids = new int[remoteGids.length-count];
                        int[] newRemoteVnodeIds = new int[remoteGids.length-count];
                        count = 0;
                        for(int i = 0; i < remoteGids.length; i++)
                            if( remoteVnodeIds[i] != -1 ) {
                                newRemoteGids[count] = remoteGids[i];
                                newRemoteVnodeIds[count++] = remoteVnodeIds[i];
                            }
                        remoteGids = newRemoteGids;
                        remoteVnodeIds = newRemoteVnodeIds;
                        this.println("ERR", "Warning in Import: Target IDs not found in Source VectorSpace (Do you want to import to subset of Target VectorSpace?)");
                    }
                    else { //valid remoteGids empty
                        remoteGids = new int[0];
                        remoteVnodeIds = new int[0];
                    }
                }
            }
            
            tmp[0] = remoteLids;
            tmp[1] = remoteGids;
            Util.sort(true, remoteVnodeIds, new double[0][0], tmp);
        }
        else {
            remoteVnodeIds = new int[0];
        }
        
    }
}