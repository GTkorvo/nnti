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

class CcjReduceDoubleMaxArray implements Reducible {

    public Serializable reduce(Serializable partialMaxs1, Serializable partialMaxs2) {
        double [] partialMaxsA = (double []) partialMaxs1;
        double [] partialMaxsB = (double []) partialMaxs2;
    
        for(int i=0; i < partialMaxsA.length; i++) {
            //the order here is very important!
            //if partialMaxsA is modified so will be the array passed
            //which may break things if the programmer is not aware of
            //that consequene
            //therefore for saftey, do not modify partialMaxsA
            if (partialMaxsB[i] < partialMaxsA[i]) {
                partialMaxsB[i] = partialMaxsA[i];
            }
        }
        
        return partialMaxsB;
    }
}