package Jpetra;

import CCJ.*;
import java.io.Serializable;

class CcjReduceDoubleMaxArray implements Reducible {

    public Serializable reduce(Serializable partialMaxs1, Serializable partialMaxs2) {
        double [] partialMaxsA = (double []) partialMaxs1;
        double [] partialMaxsB = (double []) partialMaxs2;
    
        for(int i=0; i < partialMaxsA.length; i++) {
            if (partialMaxsB[i] < partialMaxsA[i]) {
                partialMaxsB[i] = partialMaxsA[i];
            }
        }
        
        return partialMaxsB;
    }
}