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

package Jpetra.CcjSupport;

import java.io.Serializable;

import CCJ.*;

/**
 * Used by CCJ to collect all the <code>allElements</code> together during a <code>allGather</code> call.
 * Implements the CCJ interface Partitionable.
 *
 * @author Jason Cross
 */
 
public class CcjGatherSerializableArray implements Partitionable {
    
    /**
     * holds each process' elements indexed by process ID
     */
    private Serializable[][] allElements;
    
    /**
     * Creates a new CcjGatherIntArray.
     * 
     * @param groupSize number of processes in the group
     */
    public CcjGatherSerializableArray (int groupSize) {
        allElements = new Serializable[groupSize][];
    }
    
    /**
     * Used internally by CCJ.
     *
     * @param index process ID or rank
     * @param groupSize number of processes in the group
     * @param object object being added to <code>allElements[index]</code>
     */
    public void setElementAt(int index, int groupSize, Serializable object) {
        // casts object into something usefull
        allElements[index] = (Serializable[]) object;
    }
    
    /**
     * Used internally by CCJ.
     *
     * @param index process ID or rank
     * @param groupSize number of processes in the group
     * @return ArrayList from process with this <code>index</code>
     */
    public Serializable elementAt(int index, int groupSize) {
        return allElements[index];
    }

    /**
     * Takes all second degree elements of <code>allElements</code> and joins them into one int array.
     *
     * @return all arrays from all proceses put in order by process into one int array
     */
    public Serializable[] getAllElements () {
        /**
         * the size of the int array to be returned
         */
        int resultSize=0;
        
        /**
         * tracks the position of <code>returnResults</code> when performing the joining of <code>allElements</code>
         */
        int offSet = 0;
        
        // determines the size of the array so that array will be proper size
        // for hopefully somewhat faster Java array coping
        for(int i = 0; i < allElements.length; i++) {
            resultSize += allElements[i].length;
        }
        
        /**
         * contains the contiguous list of values to be returned
         */
        Serializable[] returnResults = new Serializable[resultSize];

        
        // copies each processes' int array into one contigous array <code>returnResults</code>
        for(int i = 0; i < allElements.length; i++) {
            System.arraycopy(allElements[i], 0, returnResults, offSet, allElements[i].length);
            offSet += allElements[i].length;
        }
        
        return returnResults;
    }
}
