/*
 * CcjGatherIntArray.java
 *
 * Created on June 20, 2003, 5:17 AM
 */

package Jpetra;

import java.io.Serializable;

import CCJ.*;

/**
 * Used by CCJ to collect all the <code>allElements</code> together during a <code>allGather</code> call.
 * Implements the CCJ interface Partitionable.
 *
 * @author Jason Cross
 */
 
public class CcjGatherIntArray implements Partitionable {
    
    /**
     * holds each process' elements indexed by process ID
     */
    private int[][] allElements;
    
    /**
     * Creates a new CcjGatherIntArray.
     * 
     * @param groupSize number of processes in the group
     */
    public CcjGatherIntArray (int groupSize) {
        allElements = new int[groupSize][];
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
        allElements[index] = (int[]) object;
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
    public int[] getAllElements () {
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
        int[] returnResults = new int[resultSize];

        
        // copies each processes' int array into one contigous array <code>returnResults</code>
        for(int i = 0; i < allElements.length; i++) {
            System.arraycopy(allElements[i], 0, returnResults, offSet, allElements[i].length);
            offSet += allElements[i].length;
        }
        
        return returnResults;
    }
    
    /**
     * Starts at <code>allElements[index]</code> and decraments down allElements adding all second degree values of
     * <code>allElements</code> to the second degree values of <code>allElements[index]</code>.
     
     * @param index equal to the rank of the calling process
     * @return array of scan sums
     */
    public int[] scanSums (int index) {
        for(int i=index-1; i >= 0; i--) {
            for(int j=0; j < allElements[index].length; j++) {
                allElements[index][j] += allElements[i][j];
            }
        }
    return allElements[index];
    }
}
