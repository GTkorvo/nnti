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
public class Export extends JpetraObject {
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
    
    public Export(VectorSpace sourceVectorSpace, VectorSpace targetVectorSpace) {
        //this.outputStreams.put("EXPORT", new Output("Export: ", true, System.out, false, System.out));
        
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
        int temp = sourceGids.length;
        temp = targetGids.length;
        int minNumIds = Util.min(sourceGids.length, targetGids.length);
        //this.println("EXPORT", "sourceGids.length: " + sourceGids.length + " targetGids.length: " + targetGids.length);
        for (numSameGids = 0; numSameGids < minNumIds; this.numSameGids++) {
            if (sourceGids[this.numSameGids] != targetGids[this.numSameGids]) {
                break;
            }
        }
        //this.println("EXPORT", "numSameGids: " + this.numSameGids);
        
        // next figure out how many global ids are in both the target VectorSpace and in the source VectorSpace
        int numPermuteGids = 0;
        int numExportGids = 0;
        for (int i=this.numSameGids; i < sourceGids.length; i++) {
            if (targetVectorSpace.isMyGlobalIndex(sourceGids[i])) {
                //this.println("EXPORT", "Going to permute gid: " + targetGids[i]);
                numPermuteGids++;
            }
            else {
                // if the source global id isn't found in the target VectorSpace then it must be owned by a remote vnode so we need to export it
                numExportGids++;
            }
        }
        
        // pre-allocate arrays
        if (numExportGids > 0) {
            this.exportLids = new int[numExportGids];
            this.exportGids = new int[numExportGids];
        }
        else {
            this.exportLids = new int[0];
            this.exportGids = new int[0];
        }
        if (numPermuteGids > 0)  {
            this.permuteToLids = new int[numPermuteGids];
            this.permuteFromLids = new int[numPermuteGids];
        }
        else {
            this.permuteToLids = new int[0];
            this.permuteFromLids = new int[0];
        }
        
        // now put all source global ids into an array of either those that are local and will be permuted or those that are on other vnodes
        numPermuteGids = 0;
        numExportGids = 0;
        // i is my local id
        for (int i=this.numSameGids; i < sourceGids.length; i++) {
            if (targetVectorSpace.isMyGlobalIndex(sourceGids[i])) {
                //this.println("EXPORT", "sourceGids[" + i + "]" + sourceGids[i]);
                permuteFromLids[numPermuteGids] = i;
                permuteToLids[numPermuteGids++] = targetVectorSpace.getLocalIndex(sourceGids[i]);
            }
            else {
                exportGids[numExportGids] = sourceGids[i];
                exportLids[numExportGids++] = i;
            }
        }
        
        if ((numExportGids > 0) && !sourceVectorSpace.isDistributedGlobally()) {
            this.println("WRN", "In Export: A non-distributed globally vector space has remote elements.");
        }
        
        //this.println("EXPORT", "targetGids.length: " + targetGids.length);
        //this.println("EXPORT", "permuteToLids.length: " + permuteToLids.length);
        
        // if the sourceVectorSpace is distributed globally then we need to figure out which vnodes to send out exportGids to
        int[] exportVnodeIds;
        if (sourceVectorSpace.isDistributedGlobally()) {
            int[][] tmp = targetVectorSpace.getRemoteVnodeIdList(exportGids);  // get export vnodeIds corresponding to the passed in Gids
            exportVnodeIds = tmp[0];  // tmp[0] contains the vnode  Ids, tmp[1] contains those vnodes' Lids (which we don't need)
            
            //Get rid of IDs that don't exist in the target VectorSpace
            if(exportGids.length > 0) {
                int count = 0;
                for(int i = 0; i < exportGids.length; i++) {
                    // exportVnodeIds are set to -1 if no vnode owns the corresponding Gid
                    if(exportVnodeIds[i] == -1) {
                        //this.println("EXPORT", "Gid " + exportGids[i] + " not found.");
                        count++;
                    }
                }
                if(count > 0) {
                    if(exportGids.length - count > 0) {
                        // create new arrays without the unfound Gids
                        int[] newExportGids = new int[exportGids.length-count];
                        int[] newExportVnodeIds = new int[exportGids.length-count];
                        count = 0;
                        for(int i = 0; i < exportGids.length; i++) {
                            if( exportVnodeIds[i] != -1 ) {
                                newExportGids[count] = exportGids[i];
                                newExportVnodeIds[count++] = exportVnodeIds[i];
                            }
                        }
                        exportGids = newExportGids;
                        exportVnodeIds = newExportVnodeIds;
                    }
                    else {
                        // make remoteGids and remoteVnodesIds 0 length arrays since no vnodes own any of the remote Gids we need
                        exportGids = new int[0];
                        exportVnodeIds = new int[0];
                    }
                    this.println("WRN", "In Export: Target IDs not found in Source VectorSpace (Do you want to import to subset of Target VectorSpace?)");
                }
            }
            
            // pack up exportLids and exportGids to be sorted the same way as exportVnodeIds is sorted
            tmp[0] = exportLids;
            tmp[1] = exportGids;
            // sort the Gids by vnode IDs
            Util.sort(true, exportVnodeIds, new double[0][0], tmp);
            // we will need to figure out which Gids we will be receiving and from who
            this.distributor = sourceVectorSpace.getComm().createDistributor();
            // construct a list of imports based on all the vnodes' exports by using a distributor
            distributor.createFromSends(exportVnodeIds, sourceVectorSpace.getComm());
            // pack up Gids for export
            Serializable[] packedExportGids = new Serializable[exportGids.length];
            for(int i=0; i < packedExportGids.length; i++) {
                packedExportGids[i] = new Integer(exportGids[i]);
            }
            // distribute export Gids to all vnodes that we export to and receive back a import Gids
            // from other vnodes sending to us
            Serializable[] packedRemoteGids = distributor.distribute(packedExportGids, false);
            // unpack the import Gids we received
            int numPackedRemoteGids = 0;
            for(int i=0; i < packedRemoteGids.length; i++) {
                if (packedRemoteGids[i] != null) {
                    numPackedRemoteGids += ((Serializable[]) packedRemoteGids[i]).length;
                }
            }
            this.remoteGids = new int[numPackedRemoteGids];
            Serializable[] tmpGidArray;
            int gidIndex = 0;
            for(int i=0; i < packedRemoteGids.length; i++) {
                if (packedRemoteGids[i] != null) {
                    tmpGidArray = (Serializable[]) packedRemoteGids[i];
                    for(int j=0; j < tmpGidArray.length; j++) {
                        this.remoteGids[gidIndex++] = ((Integer) tmpGidArray[j]).intValue();
                    }
                }
            }
            
            // find the export local ids for all the corresponding export global ids
            this.remoteLids = new int[this.remoteGids.length];
            for(int i=0; i < remoteLids.length; i++) {
                this.remoteLids[i] = targetVectorSpace.getLocalIndex(this.remoteGids[i]);
                //this.println("STD", "sending gid: " + this.exportGids[i] +" to vnode: " + this.exportVnodeIds[i] + " with myLid: " + this.exportLids[i]);
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