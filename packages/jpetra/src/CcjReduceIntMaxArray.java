package Jpetra;

import CCJ.*;
import java.io.Serializable;

class CcjReduceIntMaxArray implements Reducible {

    public Serializable reduce(Serializable partialMaxs1, Serializable partialMaxs2) {
        int [] partialMaxsA = (int []) partialMaxs1;
        int [] partialMaxsB = (int []) partialMaxs2;
    
        for(int i=0; i < partialMaxsA.length; i++) {
            if (partialMaxsB[i] < partialMaxsA[i]) {
                partialMaxsB[i] = partialMaxsA[i];
            }
        }
        
        return partialMaxsB;
    }
}