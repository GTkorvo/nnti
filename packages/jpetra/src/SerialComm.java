package Jpetra;

/*
 * SerialComm.java
 *
 * Created on May 29, 2001, 1:32 PM
 */


/**
 * SerialComm is the implementatin of the Comm interface
 * for machines with a single processor.
 *
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @version 
 */
public class SerialComm extends JpetraObject implements Comm {

    private int myRank;
    private int myNode;
    private int myThread;
    private int numNodes;
    private boolean isSerial;
    
    /** Creates new <code>SerialComm</code> */
    public SerialComm() {
        myNode = 0;
        numNodes = 1;
    }
    
    /**
     * Creates a new <code>SerialComm</code> with the same attributes as the argument.
     *
     * @param comm  the <code>SerialComm</code> of which to make a copy
     */
    public SerialComm(SerialComm comm) {
        myRank = comm.myRank;
        myNode = comm.myNode;
        myThread = comm.myThread;
        numNodes = comm.numNodes;
        isSerial = comm.isSerial;
    }

    /**
     * Creates and returns a clone of the calling <code>SerialComm</code>.
     *
     * @return  the clone of the calling <code>SerialComm</code>
     */
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }  
    
    public boolean getIsSerial() {
        return isSerial;
    }      
    
    public int getMyNode() {
        return myNode;
    }
    
    public int getNumNodes() {
        return numNodes;
    }
    
    public void barrier() {}
    
    public void nodeBarrier() {}
        
    public int broadcast(int numElements,double[] elements,int root) {
        return 0;
    }
    
    public int broadcast(int numElements,int[] elements,int root) {
        return 0;
    }
    
    public int gatherAll(int numElements,double[] myElements,double[] allElements) {
        for (int i=0; i<numElements; i++) allElements[i] = myElements[i];
        return 0;
    }
    
    public int gatherAll(int numElements,int[] myElements,int[] allElements) {
        for (int i=0; i<numElements; i++) allElements[i] = myElements[i];
        return 0;
    }
    
    public int sumAll(int numElements, double [] partialSums, double [] finalSums) {
        System.arraycopy(partialSums, 0, finalSums, 0, numElements);
        return 0;
    }
    
    public int sumAll(int numElements,int[] partialSums,int[] finalSums) {
        System.arraycopy(partialSums, 0, finalSums, 0, numElements);
        return 0;
    }
    
    public int maxAll(int numElements,double[] partialMaxs,double[] globalMaxs) {
        for (int i=0; i<numElements; i++) globalMaxs[i] = partialMaxs[i];
        return 0;
    }
    
    public int maxAll(int numElements,int[] partialMaxs,int[] globalMaxs) {
        for (int i=0; i<numElements; i++) globalMaxs[i] = partialMaxs[i];
        return 0;
    }
    
    public int minAll(int numElements,double[] partialMins,double[] globalMins) {
        for (int i=0; i<numElements; i++) globalMins[i] = partialMins[i];
        return 0;
    }
    
    public int minAll(int numElements,int[] partialMins,int[] globalMins) {
        for (int i=0; i<numElements; i++) globalMins[i] = partialMins[i];
        return 0;
    }
    
    public int scanSums(int numElements,double[] myElements,double[] sums) {
        for (int i=0; i<numElements; i++) sums[i] = myElements[i];
        return 0;
    }
    
    public int scanSums(int numElements,int[] myElements,int[] sums) {
        for (int i=0; i<numElements; i++) sums[i] = myElements[i];
        return 0;
    }
    
    public int getPID() {
        return 0;
    }
    
    public int getThreadID() {
        return myThread;
    }
    
    public int setThreadID(int aThread) {
        myThread = aThread;
        return 0;
    }
    
    public int setMyNodeID(int newNodeID) {
        myNode = newNodeID;
        return 0;
    }
    
    /**
     * Accessor for the number of processors in the communicator.
     */
    public int getNumProc() {
        return 1;
    }
    
}
