package Jpetra;

import CCJ.*;
import java.io.Serializable;

class CcjReduceIntMinArray implements Reducible {

    public Serializable reduce(Serializable partialMins1, Serializable partialMins2) {
        int [] partialMinsA = (int []) partialMins1;
        int [] partialMinsB = (int []) partialMins2;
    
        for(int i=0; i < partialMinsA.length; i++) {
            if (partialMinsB[i] > partialMinsA[i]) {
                partialMinsB[i] = partialMinsA[i];
            }
        }
        
        return partialMinsB;
    }
}