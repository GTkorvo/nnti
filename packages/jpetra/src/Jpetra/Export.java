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
public class Export extends JpetraObject {
    VectorSpace sourceVectorSpace;
    VectorSpace targetVectorSpace;
    
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
        
        int numSameGids;
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
        for (int i=numSameGids; i< sourceGids.length; i++) {
            if (targetVectorSpace.isMyGlobalIndex(sourceGids[i])) {
                numPermuteGids++;
            }
            else {
                numExportGids++;
            }
        }
        
        int[] exportLids;
        int[] exportGids;
        int[] permuteToLids;
        int[] permuteFromLids;
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
    
}