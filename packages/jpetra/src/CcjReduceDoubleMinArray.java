package Jpetra;

import CCJ.*;
import java.io.Serializable;

/**
 * This class defines the algorithm used by by CCJ to compare all the <code>partialMins</code> during a CCJ
 * <code>allReduce</code> call.
 * Implements the CCJ interface <code>Reducible</code>.
 *
 * @author Jason Cross
 */

class CcjReduceDoubleMinArray implements Reducible {

    public Serializable reduce(Serializable partialMins1, Serializable partialMins2) {
        double [] partialMinsA = (double []) partialMins1;
        double [] partialMinsB = (double []) partialMins2;

        //the order here is very important!
        //if partialMinsA is modified so will be the array passed
        //which may break things if the programmer is not aware of
        //that consequene
        //therefore for saftey, do not modify partialMinsA    
        for(int i=0; i < partialMinsA.length; i++) {
            if (partialMinsB[i] > partialMinsA[i]) {
                partialMinsB[i] = partialMinsA[i];
            }
        }
        
        return partialMinsB;
    }
}