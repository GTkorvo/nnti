package Jpetra;
/*
 * Util.java
 *
 * Created on May 29, 2001, 2:21 PM
 */

/**
 * @author  Michael William Boldt
 * @version 
 */
public class Util extends JpetraObject {
    public Util() {
    }

    public Util(Util util) {
    }

    public void sort(boolean sortAscending, int numKeys, int [] keys, int numDoubleCompanions, 
                     double [][] doubleCompanions, int numIntCompanions, int [][] intCompanions) {
                     
        int i;
        int n = numKeys;
        final int [] list = keys;
        int m = n / 2;
    
        while(m > 0) {
            int max = n-m;
            for(int j=0; j<max; j++) {
                for(int k=j; k>=0; k-=m) {
                    if(list[k+m] >= list[k])
                    break;
                    int temp = list[k+m];
                    list[k+m] = list[k];
                    list[k] = temp;
                    
                    for(i=0; i<numDoubleCompanions; i++) {
                        double dTemp = doubleCompanions[i][k+m];
                        doubleCompanions[i][k+m] = doubleCompanions[i][k];
                        doubleCompanions[i][k] = dTemp;
                    }
                    for(i=0; i<numIntCompanions; i++) {
                        int iTemp = intCompanions[i][k+m];
                        intCompanions[i][k+m] = intCompanions[i][k];
                        intCompanions[i][k] = iTemp;
                    }
                }
            }
            m = m / 2;
        }
    }
}

