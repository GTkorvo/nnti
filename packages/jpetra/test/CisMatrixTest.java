package test;

import Jpetra.*;

/**
 *
 * @author  Jason Cross
 */
public class CisMatrixTest extends JpetraObject {
    public static void main(String[] args) {
        new CisMatrixTest();
    }
    
    public CisMatrixTest() {
        this.initializeOutput();
        this.println("STD", "Creating SerialComm...");
        Comm myComm = new SerialComm();
        this.println("STD", "Creating ElementSpace...");
        ElementSpace rowElementSpace = new ElementSpace(5, 0, myComm);
        this.println("STD", "Creating VectorSpace...");
        VectorSpace rowVectorSpace = new VectorSpace(rowElementSpace);
        this.println("STD", "Creating CisMatrix...");
        CisMatrix rowMatrix = new CisMatrix(rowVectorSpace, CisMatrix.ROW_ORIENTED);
        
        this.println("STD", "Filling the CisMatrix...");
        
        // completely fill the matrix
        double[] temp = new double[5];
        for(int i=0; i < 5; i++) {
            for (int j=0; j < 5; j++) {
                temp[j] = i*5+j;
            }
            rowMatrix.insertEntries(i, new int[]{0,1,2,3,4}, temp, DistObject.ADD);
        }
        
        
        /*
        // user defined matrix
        rowMatrix.insertEntries(0, new int[]{1,2,3}, new double[]{2,3,4}, DistObject.ADD);
        rowMatrix.insertEntries(0, new int[]{1,2,3}, new double[]{-1,-2,-3}, DistObject.ADD);
        rowMatrix.insertEntries(1, new int[]{1,3}, new double[]{4,5}, DistObject.ADD);
        rowMatrix.insertEntries(2, new int[]{0,1,2,3,4}, new double[]{6,7,8,9,10}, DistObject.ADD);
        rowMatrix.insertEntries(3, new int[]{4}, new double[]{11}, DistObject.ADD);
        rowMatrix.insertEntries(3, new int[]{2}, new double[]{1}, DistObject.ADD);
        rowMatrix.insertEntries(4, new int[]{0,4}, new double[]{12,13}, DistObject.ADD);
        */
        
        this.println("STD", "Calling fillComplete() on CisMatrix...");
        rowMatrix.fillComplete();
        
        this.println("STD", "Calling scale(10) on CisMatrix...");
        rowMatrix.scale(10);
        this.println("STD", "Calling printOut() on CisMatrix...");
        rowMatrix.printOut("STD");
        this.println("STD", "Done.");
        
    }
}
