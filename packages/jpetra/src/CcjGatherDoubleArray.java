/*
 * CcjGatherDoubleArray.java
 *
 * Created on June 20, 2003, 5:17 AM
 */

package Jpetra;

import java.io.Serializable;

import CCJ.*;

/**
 * Used by CCJ to collect all the <code>allElements</code> together during a CCJ <code>allGather<code> call.
 * Implements the CCJ interface Partitionable.
 *
 * @author Jason Cross
 */
 
public class CcjGatherDoubleArray implements Partitionable {
    
    /**
     * holds each process' elements indexed by process ID
     */
    private double[][] allElements;
    
    /**
     * Creates a new CcjGatherDoubleArray.
     * 
     * @param groupSize number of processes in the group
     */
    public CcjGatherDoubleArray (int groupSize) {
        allElements = new double[groupSize][];
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
        allElements[index] = (double[]) object;
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
    public double[] returnAllElements () {
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
        double[] returnResults = new double[resultSize];

        
        // copies each processes' int array into one contigous array <code>returnResults</code>
        for(int i = 0; i < allElements.length; i++) {
            System.arraycopy(allElements[i], 0, returnResults, offSet, allElements[i].length);
            offSet += allElements[i].length;
        }
        
        return returnResults;
    }
    
    /**
     * Copies <code>allElements[0]</code> and then compares the rest of the elements of<code>allElements</code>
     * to the first element array to find the max values for each second degree element.
     *
     * @return array with all the second degree max values
     */
    public double[] getMaxs () {
        //copy allElements[0]
        double[] doubleMaxs = new double[allElements[0].length];
        System.arraycopy(allElements[0], 0, doubleMaxs, 0, allElements[0].length);
        
        //find the maxs
        //step through each element of allElemenents
        for(int i=1; i < allElements.length; i++) {
            //step through each second degree element of current allElements element
            for(int j=0; j < allElements[i].length; j++) {
                //whichever value is the max is stored to doubleMaxs
                if (doubleMaxs[j] < allElements[i][j]) {
                    doubleMaxs[j] = allElements[i][j];
                }
            }
        }
        
        return doubleMaxs;
    }

    /**
     * Copies <code>allElements[0]</code> and then compares the rest of the elements of<code>allElements</code>
     * to the first element array to find the min values for each second degree element.
     *
     * @return array with all the second degree min values
     */
    //see double[] getMaxs() for code comments    
    public double[] getMins () {
        double[] doubleMins = new double[allElements[0].length];
        System.arraycopy(allElements[0], 0, doubleMins, 0, allElements[0].length);
        for(int i=1; i < allElements.length; i++) {
            for(int j=0; j < allElements[i].length; j++) {
                if (doubleMins[j] > allElements[i][j]) {
                    doubleMins[j] = allElements[i][j];
                }
            }
        }
        
        return doubleMins;
    } 
    
    /**
     * Starts at <code>allElements[index]</code> and decraments down allElements adding all second degree values of
     * <code>allElements</code> to the second degree values of <code>allElements[index]</code>.
     
     * @param index equal to the rank of the calling process
     * @return array of scan sums
     */
    public double[] scanSums (int index) {
        for(int i=index-1; i >= 0; i--) {
            for(int j=0; j < allElements[index].length; j++) {
                allElements[index][j] += allElements[i][j];
            }
        }
    return allElements[index];
    }   
}
