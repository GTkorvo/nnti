/*
 * CcjGatherDoubleArray.java
 *
 * Created on June 20, 2003, 7:58 AM
 */

package Jpetra;

import CCJ.*;
import java.io.Serializable;

/**
 * This class defines the algorithm used by by CCJ to sum all the <code>partialSums</code> together during a CCJ
 * <code>allReduce</code> call.
 * Implements the CCJ interface Reducible.
 *
 * @author Jason Cross
 */

class CcjReduceAllSumIntArray implements Reducible {

    public Serializable reduce(Serializable partialSums1, Serializable partialSums2) {

        //cast the partialSums into a useful object type 
        int [] partialSumsA = (int []) partialSums1;
        int [] partialSumsB = (int []) partialSums2;
    
        //add the partial sums together
        for(int i=0; i < partialSumsA.length; i++) {
            partialSumsB[i] += partialSumsA[i];
        }
        
        //returned the summed partial sums
        return partialSumsB;
    }
}