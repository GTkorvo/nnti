/*
 * PerformanceTest1.java
 *
 * Created on July 12, 2004, 11:58 AM
 */

package test;

import Jpetra.*;

import java.lang.Runtime;

/**
 *
 * @author  Jason Cross
 */
public class PerformanceTest1 extends JpetraObject {
    
    /** Creates a new instance of PerformanceTest1 */
    public PerformanceTest1(int numElementsPerVnode, int numColumnsA, int numVectorsInX) {
        initializeOutput();
        Comm comm = new CcjComm("test/ccjhosts.txt");
        setAmIroot(comm);
        //setVnodesPrint("STD", true);
        this.outputStreams.put("MAIN", new Output("-->", true, System.out, false, System.out));
        this.println("MAIN", "Comm is setup, creating objects for test.");
        
        int myVnodeId = comm.getVnodeId();
        
        ElementSpace aElementSpace = new ElementSpace(numElementsPerVnode * comm.getNumVnodes(), numElementsPerVnode, 0, comm);
        CisMatrix A = new CisMatrix(new VectorSpace(aElementSpace), CisMatrix.ROW_ORIENTED);
        
        int startValue = numElementsPerVnode * myVnodeId;
        double[] values = new double[numColumnsA];
        int[] indices = new int[numColumnsA];
        for(int col=0; col < numColumnsA; col++) {
            indices[col] = col;
        }
        for(int row=0; row < numElementsPerVnode; row++) {
            for(int col=0; col < numColumnsA; col++) {
                values[col] = startValue++;
            }
            
            A.insertEntries(row, indices, values, DistObject.REPLACE);
        }
        
        A.fillComplete();
        this.println("MAIN", "A construction is complete.");
        
        double[][] xValues = new double[numVectorsInX][numColumnsA];
        for(int i=0; i < xValues.length; i++) {
            for(int j=0; j < xValues[i].length; j++) {
                xValues[i][j] = 3;
            }
        }
        
        ElementSpace xElementSpace = new ElementSpace (numColumnsA, numColumnsA, 0, comm);
        MultiVector x = new MultiVector(new VectorSpace(xElementSpace), xValues);
        this.println("MAIN", "x construction is complete.");
        
        int numYelements = numElementsPerVnode * comm.getNumVnodes();
        double[][] yValues = new double[numVectorsInX][numElementsPerVnode];
        //ElementSpace yElementSpace = new ElementSpace(numElementsPerVnode, numYelements, 0, comm);
        MultiVector y = new MultiVector(new VectorSpace(aElementSpace), yValues);
        this.println("MAIN", "y construction is complete.");
        
        A.setFlopCounter(new FlopCounter());
        this.println("MAIN", "A has been given a FlopCounter.");
        this.println("MAIN", "do y=A*x");
        Time myTime = new Time(comm);
        
        A.multiply(CisMatrix.USE_A, x, y);
        
        double elapsedTime = myTime.getElapsedTime();
        double[] flopsDone = comm.gatherAll(new double[]{A.getFlops()});
        Runtime runtime = java.lang.Runtime.getRuntime();
        double totalMemory = runtime.totalMemory();
        double freeMemory = runtime.freeMemory();
        double[] memoryInfo = comm.gatherAll(new double[]{totalMemory/(1024*1024), freeMemory/(1024*1024)});
        
        this.println("MAIN", "y=A*x complete.");
        this.println("MAIN", "operation took " + elapsedTime + " seconds and did " + flopsDone[0] + " flops.");
        int count = 0;
        for(int i=0; i < comm.getNumVnodes(); i++) {
            this.println("MAIN", "mem info for vnode " + i);
            this.println("MAIN", "  totalMemory=" + memoryInfo[count] + " freeMemory=" + memoryInfo[count+1] + " usedMemory=" + (memoryInfo[count] - memoryInfo[count+1]));
            count += 2;
        }
        this.println("MAIN", "DONE");
        System.exit(0);
    }
    
    public static void main (String[] args) {
        if (args.length != 3) {
            System.out.println("Usage: java test/PerformanceTest1 numRowsA numColsA numVectorsX");
            System.exit(-1);
        }
        new PerformanceTest1(Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]));
    }
}
