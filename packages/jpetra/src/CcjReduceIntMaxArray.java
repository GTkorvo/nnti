package Jpetra;

import CCJ.*;
import java.io.Serializable;

/**
 * This class defines the algorithm used by by CCJ to compare all the <code>partialMaxs</code> during a CCJ
 * <code>allReduce</code> call.
 * Implements the CCJ interface <code>Reducible</code>.
 *
 * @author Jason Cross
 */

class CcjReduceIntMaxArray implements Reducible {

    public Serializable reduce(Serializable partialMaxs1, Serializable partialMaxs2) {
        int [] partialMaxsA = (int []) partialMaxs1;
        int [] partialMaxsB = (int []) partialMaxs2;
    
        //the order here is very important!
        //if partialMaxsA is modified so will be the array passed
        //which may break things if the programmer is not aware of
        //that consequene
        //therefore for saftey, do not modify partialMaxsA
        for(int i=0; i < partialMaxsA.length; i++) {
            if (partialMaxsB[i] < partialMaxsA[i]) {
                partialMaxsB[i] = partialMaxsA[i];
            }
        }
        
        return partialMaxsB;
    }
}