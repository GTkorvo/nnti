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

import java.util.TreeSet;

/**
 *
 * @author  Jason Cross
 */
public class Export extends JpetraObject {
    VectorSpace sourceVectorSpace;
    VectorSpace targetVectorSpace;
    int numSameGids;
    int numImports;
    int[] exportLids;
    int[] exportGids;
    int[] permuteToLids;
    int[] permuteFromLids;
    int[] exportVnodeIds;
    Distributor distributor;
    
    public Export(VectorSpace sourceVectorSpace, VectorSpace targetVectorSpace) {
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
        int numExportGids = 0;
        for (int i=numSameGids; i < sourceGids.length; i++) {
            if (targetVectorSpace.isMyGlobalIndex(sourceGids[i])) {
                numPermuteGids++;
            }
            else {
                numExportGids++;
            }
        }
        
        //int[] exportLids;
        //int[] exportGids;
        //int[] permuteToLids;
        //int[] permuteFromLids;
        if (numExportGids > 0) {
            exportLids = new int[numExportGids];
            exportGids = new int[numExportGids];
        }
        else {
            exportLids = new int[0];
            exportGids = new int[0];
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
        numExportGids = 0;
        for (int i=numSameGids; i < sourceGids.length; i++) {
            if (targetVectorSpace.isMyGlobalIndex(sourceGids[i])) {
                this.println("STD", "sourceGids[i]" + sourceGids[i]);
                permuteFromLids[numPermuteGids] = i;
                permuteToLids[numPermuteGids++] = targetVectorSpace.getLocalIndex(sourceGids[i]);
            }
            else {
                //NumSend_ +=SourceMap.ElementSize(i); // Count total number of entries to send
                //NumSend_ +=SourceMap.MaxElementSize(); // Count total number of entries to send (currently need max)
                exportGids[numExportGids] = sourceGids[i];
                exportLids[numExportGids++] = i;
            }
        }
        
        //TreeSet importVnodeIds = new TreeSet();
        this.println("STD", "targetGids.length: " + targetGids.length);
	this.println("STD", "permuteToLids.length: " + permuteToLids.length);
        numImports = targetGids.length - numSameGids - permuteToLids.length;
        /*for (int i=numSameGids; i < targetGids.length; i++) {
            if (!sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                numImports++;
            }
        }*/
        
        int[] importGids = new int[numImports];
        int importGidsIndex = 0;
        this.println("STD", "importGids.length: " + importGids.length);
        for (int i=numSameGids; i < targetGids.length; i++) {
            if (!sourceVectorSpace.isMyGlobalIndex(targetGids[i])) {
                importGids[importGidsIndex++] = targetGids[i];
            }
        }
        
        //int[] exportVnodeIds;
        int[] importVnodeIds;
        if (sourceVectorSpace.isDistributedGlobally()) {
            exportVnodeIds = new int[exportGids.length];
            importVnodeIds = new int[importGids.length];
            int[][] tmp = targetVectorSpace.getRemoteVnodeIdList(exportGids);  // finds vnodeIds for exportGids
            exportVnodeIds = tmp[0];
            tmp = sourceVectorSpace.getRemoteVnodeIdList(importGids);  // finds vnodeIds for importGids
            importVnodeIds = tmp[0];
        }
        else {
            exportVnodeIds = new int[0];
            importVnodeIds = new int[0];
        }
        
        //Get rid of IDs not in Target Map
        if(exportVnodeIds.length > 0) {
            int count = 0;
            for(int i = 0; i < exportVnodeIds.length; i++ ) {
                if( exportVnodeIds[i] == -1 ) {
                    count++;
                }
            }
            if (count > 0) {
                int[] newExportGids = null;
                int[] newExportVnodeIds = null;
                int count1 = exportVnodeIds.length - count;
                if (count1 > 0) {
                    newExportGids = new int[count1];
                    newExportVnodeIds = new int[count1];
                }
                count = 0;
                for(int i = 0; i < exportVnodeIds.length; i++)
                    if(exportVnodeIds[i] != -1) {
                        newExportGids[count] = exportGids[i];
                        newExportVnodeIds[count] = exportVnodeIds[i];
                        count++;
                    }
                exportGids = newExportGids;
                exportVnodeIds = newExportVnodeIds;
                this.println("STD", "Warning in Export: Source IDs not found in Target Map (Do you want to export from subset of Source Map?");
            }
        }
        
        int[][] tmp = new int[2][];
        tmp[0] = exportLids;
        tmp[1] = exportGids;
        Util.sort(true, exportVnodeIds, new double[0][0], tmp);
        
        tmp[0] = new int[0];
        tmp[1] = new int[0];
        Util.sort(true, importVnodeIds, new double[0][0], tmp);
        
        
        this.println("STD", "numImports: " + numImports);
        
        distributor = sourceVectorSpace.getComm().createDistributor();
        distributor.createFromReceives(importVnodeIds);
        
        // need to do distributor.do
        
        //for (i=0; i< remoteGids.length; i++) {
        //    remoteLids[i] = targetVectorSpace.getLocalIndex(remoteGids[i]);
        //}
        
        this.println("STD", "LIDS  GIDS");
        this.println("STD", "----------");
        for(int i=0; i < exportLids.length; i++) {
            this.println("STD", exportLids[i] + " " + exportGids[i]);
        }
        
        this.println("STD", "ToLIDS  FromLIDS");
        this.println("STD", "----------------");
        for(int i=0; i < permuteToLids.length; i++) {
            this.println("STD", permuteToLids[i] + " " + permuteFromLids[i]);
        }
        
    }
    
    public Distributor getDistributor() {
        return this.distributor;
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
}