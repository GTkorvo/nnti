package Jpetra;

import CCJ.*;
import java.io.Serializable;

class CcjReduceDoubleMinArray implements Reducible {

    public Serializable reduce(Serializable partialMins1, Serializable partialMins2) {
        double [] partialMinsA = (double []) partialMins1;
        double [] partialMinsB = (double []) partialMins2;
    
        for(int i=0; i < partialMinsA.length; i++) {
            if (partialMinsB[i] > partialMinsA[i]) {
                partialMinsB[i] = partialMinsA[i];
            }
        }
        
        return partialMinsB;
    }
}