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
    static int numElementsPerVnode;
    static int numColumnsA;
    static int numVectorsInX;
    static String inFileName;
    
    /** Creates a new instance of PerformanceTest1 */
    public PerformanceTest1() {
        Runtime runtime = java.lang.Runtime.getRuntime();
        initializeOutput();
        Comm comm = new CcjComm("test/ccjhosts.txt");
        setAmIroot(comm);
        //setVnodesPrint("STD", true);
        setRootPrint("STD", false);
        this.outputStreams.put("MAIN", new Output("-->", true, System.out, true, System.out));
        this.println("MAIN", "Comm is setup, creating objects for test.");
        
        int myVnodeId = comm.getVnodeId();
        
        //dense fill
        /*ElementSpace aElementSpace = new ElementSpace(numElementsPerVnode * comm.getNumVnodes(), numElementsPerVnode, 0, comm);
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
            //this.println("MAIN", "freeMemory=" + (runtime.freeMemory()/(1024*1024)) + " totalMemory=" + (runtime.totalMemory()/(1024*1024)) + " usedMemory=" + ((runtime.totalMemory() - runtime.freeMemory())/(1024*1024)));
            //System.out.flush();
            A.insertEntries(row, indices, values, DistObject.REPLACE);
        }*/
        
        //matrix market reader
        CisMatrix inMatrix = null;
        try {
            this.println("MAIN", "reading in inMatrix from " + inFileName + "...");
            inMatrix = Jpetra.MatrixMarketIO.CisMatrixReader.read(inFileName, CisMatrix.ROW_ORIENTED, comm);
        } catch (java.io.IOException e) {
            this.println("FATALERR", e.toString());
            System.exit(1);
        }
        this.println("MAIN", "inMatrix.getNumMyNonZeros()=" + inMatrix.getNumMyNonZeros());
        this.println("MAIN", "inMatrix has been read in.");
        VectorSpace globalVectorSpace = new VectorSpace(new ElementSpace(inMatrix.getPrimaryVectorSpace().getNumGlobalEntries(), 0, comm));
        this.println("MAIN", "globalVectorSpace has been created.");
        CisMatrix A = new CisMatrix(globalVectorSpace, CisMatrix.ROW_ORIENTED);
        Export exporter = new Export(inMatrix.getPrimaryVectorSpace(), globalVectorSpace);
        A.exportValues(inMatrix, exporter, DistObject.REPLACE);
        this.println("MAIN", "inMatrix to A export is done.");
        // make sure the GC knows that inMatrix is no longer needed
        inMatrix = null;
        A.fillComplete();
        
        //A.printOutAllVnodes("MAIN");
        this.println("MAIN", "A construction is complete.");
        
        //double[][] xValues = new double[numVectorsInX][numColumnsA];
        double[][] xValues = new double[numVectorsInX][A.getNumGlobalColumns()];
        for(int i=0; i < xValues.length; i++) {
            for(int j=0; j < xValues[i].length; j++) {
                xValues[i][j] = 3;
            }
        }
        
        //ElementSpace xElementSpace = new ElementSpace (numColumnsA, numColumnsA, 0, comm);
        ElementSpace xElementSpace = new ElementSpace(A.getNumGlobalColumns(), A.getNumGlobalColumns(), 0, comm);
        MultiVector x = new MultiVector(new VectorSpace(xElementSpace), xValues);
        this.println("MAIN", "x construction is complete.");
        
        //int numYelements = numElementsPerVnode * comm.getNumVnodes();
        //double[][] yValues = new double[numVectorsInX][numElementsPerVnode];
        double[][] yValues = new double[numVectorsInX][A.getPrimaryVectorSpace().getNumMyGlobalEntries()];
        //ElementSpace yElementSpace = new ElementSpace(numElementsPerVnode, numYelements, 0, comm);
        //MultiVector y = new MultiVector(new VectorSpace(aElementSpace), yValues);
        MultiVector y = new MultiVector(globalVectorSpace, yValues);
        this.println("MAIN", "y construction is complete.");
        
        A.setFlopCounter(new FlopCounter());
        this.println("MAIN", "A has been given a FlopCounter.");
        this.println("MAIN", "do y=A*x");
        Time myTime = new Time(comm);
        
        A.multiply(CisMatrix.USE_A, x, y);
        
        double elapsedTime = myTime.getElapsedTime();
        double flops = A.getFlops();
        this.println("MAIN", myVnodeId + ": numFlops=" + flops);
        double[] flopsDone = comm.sumAll(new double[]{flops});
        
        double totalMemory = runtime.totalMemory();
        double freeMemory = runtime.freeMemory();
        double[] memoryInfo = comm.gatherAll(new double[]{totalMemory/(1024*1024), freeMemory/(1024*1024)});
        
        //y.printOutAllVnodes("MAIN");
        this.println("MAIN", "y=A*x complete.");
        this.println("MAIN", "operation took " + elapsedTime + " seconds and did " + (flopsDone[0]/(1000000)) + " Mflops which is " + ((flopsDone[0]/elapsedTime)/1000000) + " Mflops/sec.");
        int count = 0;
        for(int i=0; i < comm.getNumVnodes(); i++) {
            this.println("MAIN", "mem info for vnode " + i);
            this.println("MAIN", "  totalMemory=" + memoryInfo[count] + " freeMemory=" + memoryInfo[count+1] + " usedMemory=" + (memoryInfo[count] - memoryInfo[count+1]));
            count += 2;
        }
        this.println("MAIN", "DONE");
        System.exit(0);
    }
    
    public static void main(String[] args) {
        boolean DO_FROM_FILE = true;
        boolean DO_DENSE = false;
        
        if (DO_FROM_FILE) {
            if (args.length != 2) {
                System.out.println("Usage: java test/PerformanceTest1 matrixMarketFormatFile.mtx numVectorsX");
                System.exit(-1);
            }
            else {
                inFileName = args[0];
                numVectorsInX = Integer.parseInt(args[1]);
            }
        }
        
        if (DO_DENSE) {
            if (args.length != 3) {
                System.out.println("Usage: java test/PerformanceTest1 numRowsA numColsA numVectorsX");
                System.exit(-1);
            } else {
                numElementsPerVnode = Integer.parseInt(args[0]);
                numColumnsA = Integer.parseInt(args[1]);
                numVectorsInX = Integer.parseInt(args[2]);
            }
        }
        
        new PerformanceTest1();
    }
}
