/*
 * SerialComm.java
 *
 * Created on May 19, 2004, 5:39 PM
 */

package Jpetra;

/**
 *
 * @author  jacross
 */
public class SerialComm implements Comm {
    
    /** Creates a new instance of SerialComm */
    public SerialComm() {
    }
    
    public void barrier() {
        // blank
    }
    
    public int broadcast(int value, int root) {
        return value;
    }
    
    public double broadcast(double value, int root) {
        return value;
    }
    
    public java.io.Serializable broadcast(java.io.Serializable value, int root) {
        return value;
    }
    
    public double[] gatherAll(double myDouble) {
        double[] temp = {myDouble};
        return temp;
    }
    
    public double[] gatherAll(double[] myElements) {
    }
    
    public java.io.Serializable[] gatherAll(java.io.Serializable[] myElements) {
    }
    
    public int[] gatherAll(int[] myElements) {
    }
    
    public int[] gatherAll(int myInt) {
        int[] temp = {myInt};
        return temp;
    }
    
    public boolean getIsSerial() {
        return true;
    }
    
    public int getNumMyThreads() {
        return 1;
    }
    
    public int getNumThreads() {
        return 1;
    }
    
    public int getNumVnodes() {
        return 1;
    }
    
    public int getThreadID() {
        return 0;
    }
    
    public int getVnodeID() {
        return 0;
    }
    
    public double[] maxAll(double[] partialMaxs) {
        return partialMaxs;
    }
    
    public int[] maxAll(int[] partialMaxs) {
        return partialMaxs;
    }
    
    public int[] minAll(int[] partialMins) {
        return partialMins;
    }
    
    public double[] minAll(double[] partialMins) {
        return partialMins;
    }
    
    public double[] scanSums(double[] myElements) {
    }
    
    public int[] scanSums(int[] myElements) {
    }
    
    public int setMyVnodeID(int newVnodeID) {
        // empty
    }
    
    public int setThreadID(int newThreadID) {
        // empty
    }
    
    public int[] sumAll(int[] partialSums) {
        return partialSums;
    }
    
    public double[] sumAll(double[] partialSums) {
        return partialSums;
    }
    
    public void threadBarrier() {
        // blank
    }
    
}
