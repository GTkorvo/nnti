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

/**
 *
 * @author  Jason Cross
 */
public class Import extends JpetraObject {
    VectorSpace sourceVectorSpace;
    VectorSpace targetVectorSpace;
    
    // for data I will receive
    private int numSameGids;
    private int[] remoteLids;  // another vnode's lid for the gid I need
    private int[] remoteGids;  // gid I need not on my vnode
    private int[] permuteToLids;
    private int[] permuteFromLids;
    
    // for data I will send
    private int[] exportLids;
    private int[] exportGids;
    private int[] exportVnodeIds;
    
    private Distributor distributor;
    
    public Import(VectorSpace sourceVectorSpace, VectorSpace targetVectorSpace) {
        //this.outputStreams.put("IMPORT", new Output("Import: ", true, System.out, false, System.out));
        
        this.sourceVectorSpace = sourceVectorSpace;
        this.targetVectorSpace = targetVectorSpace;
        
        // first get all the global IDs from the source and target VectorSpaces and put them into arrays we can use
        int[] sourceGids;
        if (sourceVectorSpace.getNumMyGlobalEntries() > 0) {
            sourceGids = sourceVectorSpace.getMyGlobalEntryIds();
        }
        else {
            // setting sourceGids to a 0 length array allows for loops to work properly
            sourceGids = new int[0];
        }
        
        int[] targetGids;
        if (targetVectorSpace.getNumMyGlobalEntries() > 0) {
            targetGids = targetVectorSpace.getMyGlobalEntryIds();
        }
        else {
            targetGids = new int[0];
        }
        
        // now figure out how many global ids are the same between the source and target VectorSpaces
        int minNumIds = Util.min(sourceGids.length, targetGids.length);
        //this.println("IMPORT", "sourceGids.length: " + sourceGids.length + " targetGids.length: " + targetGids.length);
        for (numSameGids = 0; numSameGids < minNumIds; this.numSameGids++) {
            if (sourceGids[this.numSameGids] != targetGids[this.numSameGids]) {
                break;
            }
        }
        //this.println("IMPORT", "numSameGids: " + this.numSameGids);
        
        // next figure out how many global ids are in both the target VectorSpace and in the source VectorSpace
        int numPermuteGids = 0;
        int numRemoteGids = 0;
        for (int i=this.numSameGids; i < targetGids.length; i++) {
            if (sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                //this.println("IMPORT", "Going to permute gid: " + targetGids[i]);
                numPermuteGids++;
            }
            else {
                // if the target global id isn't found in the source VectorSpace then it must be owned by a remote vnode
                numRemoteGids++;
            }
        }
        
        // pre-allocate arrays
        if (numRemoteGids > 0) {
            this.remoteLids = new int[numRemoteGids];
            this.remoteGids = new int[numRemoteGids];
        }
        else {
            this.remoteLids = new int[0];
            this.remoteGids = new int[0];
        }
        if (numPermuteGids > 0)  {
            this.permuteToLids = new int[numPermuteGids];
            this.permuteFromLids = new int[numPermuteGids];
        }
        else {
            this.permuteToLids = new int[0];
            this.permuteFromLids = new int[0];
        }
        
        // now put all target global ids into an array of either those that are local and will be permuted or those that are on other vnodes
        numPermuteGids = 0;
        numRemoteGids = 0;
        // i is my local id
        for (int i=this.numSameGids; i < targetGids.length; i++) {
            if (sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                permuteToLids[numPermuteGids] = i;
                permuteFromLids[numPermuteGids++] = sourceVectorSpace.getLocalIndex(targetGids[i]);
            }
            else {
                remoteGids[numRemoteGids] = targetGids[i];
                remoteLids[numRemoteGids++] = i;
            }
        }
        
        if ((numRemoteGids > 0) && !sourceVectorSpace.isDistributedGlobally()) {
            this.println("WRN", "A non-distributed globally vector space has remote elements.");
        }
        
        //this.println("IMPORT", "targetGids.length: " + targetGids.length);
        //this.println("IMPORT", "permuteToLids.length: " + permuteToLids.length);
        
        // if the sourceVectorSpace is distributed globally then we have to figure out which remote vnodes own the global ids we need
        int[] remoteVnodeIds;
        if (sourceVectorSpace.isDistributedGlobally()) {
            int[][] tmp = sourceVectorSpace.getRemoteVnodeIdList(remoteGids);  // get remote vnodeIds corresponding to the passed in Gids
            remoteVnodeIds = tmp[0];  // tmp[0] contains the vnode  Ids, tmp[1] contains those vnodes' Lids (which we don't need)
            
            //Get rid of IDs that don't exist in the source VectorSpace
            if(remoteGids.length > 0) {
                int count = 0;
                for(int i = 0; i < remoteGids.length; i++) {
                    // remoteVnodeIds are set to -1 if no vnode owns the corresponding Gid
                    if( remoteVnodeIds[i] == -1 ) {
                        count++;
                    }
                }
                if( count > 0 ) {
                    if( remoteGids.length - count > 0 ) {
                        // create new arrays without the unfound Gids
                        int[] newRemoteGids = new int[remoteGids.length-count];
                        int[] newRemoteVnodeIds = new int[remoteGids.length-count];
                        count = 0;
                        for(int i = 0; i < remoteGids.length; i++) {
                            if( remoteVnodeIds[i] != -1 ) {
                                newRemoteGids[count] = remoteGids[i];
                                newRemoteVnodeIds[count++] = remoteVnodeIds[i];
                            }
                        }
                        remoteGids = newRemoteGids;
                        remoteVnodeIds = newRemoteVnodeIds;
                    }
                    else {
                        // make remoteGids and remoteVnodesIds 0 length arrays since no vnodes own any of the remote Gids we need
                        remoteGids = new int[0];
                        remoteVnodeIds = new int[0];
                    }
                    this.println("WRN", "In Import: Target IDs not found in Source VectorSpace (Do you want to import to subset of Target VectorSpace?)");
                }
            }
            
            // pack up remoteLids and remoteGids to be sorted the same way as remoteVnodeIds is sorted
            tmp[0] = remoteLids;
            tmp[1] = remoteGids;
            // sort the Gids by vnode IDs
            Util.sort(true, remoteVnodeIds, new double[0][0], tmp);
            // we will need to figure out what Gids to export and which vnodes have the Gids we need
            // so we'll use a distributor to do so
            this.distributor = sourceVectorSpace.getComm().createDistributor();
            // construct a list of exports that calling vnode needs to send as a result
            // of everyone asking for what it needs to receive.
            this.exportGids = distributor.createFromReceives(remoteGids, remoteVnodeIds, sourceVectorSpace.getComm());
            this.exportLids = new int[this.exportGids.length];
            this.exportVnodeIds = distributor.getExportVnodeIds();
            // find the export local ids for all the corresponding export global ids
            for(int i=0; i < exportLids.length; i++) {
                this.exportLids[i] = sourceVectorSpace.getLocalIndex(this.exportGids[i]);
                //this.println("IMPORT", "sending gid: " + this.exportGids[i] +" to vnode: " + this.exportVnodeIds[i] + " with myLid: " + this.exportLids[i]);
            }
        }  // end if (sourceVectorSpace.isDistributedGlobally())
    }
    
    
    public int getNumSameGids() {
        return this.numSameGids;
    }
    public int[] getRemoteLids() {
        return this.remoteLids;  // another vnode's lids for the gids I need
    }
    public int[] getRemoteGids() {
        return this.remoteGids;  // gids I need not on my vnode
    }
    public int[] getPermuteToLids() {
        return this. permuteToLids;
    }
    public int[] getPermuteFromLids() {
        return this.permuteFromLids;
    }
    
    
    public int[] getExportLids() {
        return this.exportLids;
    }
    public int[] getExportGids() {
        return this.exportGids;
    }
    public int[] getExportVnodeIds() {
        return this.exportVnodeIds;
    }
    public Distributor getDistributor() {
        return this.distributor;
    }
}