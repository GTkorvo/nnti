package test;

import Jpetra.*;
/**
 *
 * @author  Jason Cross
 */
public class MultiVectorTest extends JpetraObject {
    public static void main(String[] args) {
        new MultiVectorTest();
    }
    
    public MultiVectorTest() {
        this.initializeOutput();
        this.println("STD", "Creating SerialComm...");
        Comm myComm = new SerialComm();
        this.println("STD", "Creating ElementSpace...");
        ElementSpace rowElementSpace = new ElementSpace(3, myComm);
        this.println("STD", "Creating VectorSpace...");
        VectorSpace rowVectorSpace = new VectorSpace(rowElementSpace);
        this.println("STD", "Creating and filling two MultiVectors...");
        double[][] xData = new double[3][];
        xData[0] = new double[]{0,1,2,3,4,5};
        xData[1] = new double[]{0,1,2,3,4,5};
        xData[2] = new double[]{1,2,3,4,5,6};
        
        double[][] yData = new double[3][];
        yData[0] = new double[]{0,-1,-2,-3,-4,-5};
        yData[1] = new double[]{0,1,2,3,4,5};
        yData[2] = new double[]{1,2,3,4,5,0''};
        
        MultiVector x = new MultiVector(rowVectorSpace, xData);
        MultiVector y = new MultiVector(rowVectorSpace, yData);
        
        this.println("STD", "Doing x.dot(y)...");
        double[] result = x.dot(y);
        this.println("STD", "The results is:\n" + result[0] + "\n" + result[1] + "\n" + result[2]);
    }
    
}
