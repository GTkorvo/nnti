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
        //this.println("STD", "Creating SerialComm...");
        //Comm myComm = new SerialComm();
        this.println("STD", "Creating CcjComm...");
        Comm myComm = new CcjComm("test/ccjhosts.txt");
        int myVnodeId = myComm.getVnodeId();
        this.println("STD", "I'm vnode #" + myVnodeId);
        this.println("STD", "Creating ElementSpace...");
        ElementSpace rowElementSpace = new ElementSpace(12, 0, myComm);
        this.println("STD", "Creating VectorSpace...");
        VectorSpace rowVectorSpace = new VectorSpace(rowElementSpace);
        this.println("STD", "Creating and filling two MultiVectors...");
        double[][] xData = new double[3][];
        xData[0] = new double[]{0,1,5,3,4,-5};
        xData[1] = new double[]{0,1,10,3,4,5};
        xData[2] = new double[]{1,2,3,4,5,6};
        
        double[][] yData = new double[3][];
        yData[0] = new double[]{0,-1,-2,-3,-4,-5};
        yData[1] = new double[]{0,1,2,3,4,5};
        yData[2] = new double[]{1,2,15,4,5,0};
        
        MultiVector x = null;
        if (myVnodeId == 0) {
            x = new MultiVector(rowVectorSpace, xData);
        }
        else if (myVnodeId == 1) {
            x = new MultiVector(rowVectorSpace, yData);
        }
        //MultiVector x = new MultiVector(rowVectorSpace, xData);
        //MultiVector y = new MultiVector(rowVectorSpace, yData);
        
        this.println("STD", "Doing x.scale(10.0)");
        x.scale(10.0);
        x.printOut("STD");
        
        //this.println("STD", "Doing x.dot(y)...");
        //double[] result = x.dot(y);
        //this.println("STD", "The results is:\n" + result[0] + "\n" + result[1] + "\n" + result[2]);
        //this.println("STD", "Doing x.norm2()");
        //double[] result = x.norm2();
        //this.println("STD", "The results is:\n" + result[0] + "\n" + result[1] + "\n" + result[2]);
        //this.println("STD", "Doing x.maxValue()");
        //double[] result = x.norm2();
        //this.println("STD", "The results is:\n" + result[0] + "\n" + result[1] + "\n" + result[2]);
        
        System.exit(0);  // must call explicitly or else program will not end if using CcjComm
    }
    
}
