// @HEADER
// ***********************************************************************
// 
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

package Jpetra;

/**
 * <code>SerialComm</comm> is a simple class that implements the <code>Comm</code>
 * interface for serial communication.  Most methods either simply return the
 * parameters passed to them or return a hard coded value.
 *
 * @author  Jason Cross
 */
public class SerialComm extends JpetraObject implements Comm {
    
    /** Creates a new instance of SerialComm */
    public SerialComm() {
        // empty
    }
    
    /**
     * Does nothing.
     */
    public void barrier() {
        // empty
    }
    
    /**
     * Does a broadcast of <code>value</code> to the node itself.
     *
     * @param value The int to broadcast to all nodes.
     * @param root The root vnode to gather the ints onto.  This is not used for serial communication and is ignored.
     * @return <code>value</code> that was passed in
     */
    public int broadcast(int value, int root) {
        return value;
    }
    
    /**
     * Does a broadcast of <code>value</code> to the node itself.
     *
     * @param value The double to broadcast to all nodes.
     * @param root The root vnode to gather the doubles onto.  This is not used for serial communication and is ignored.
     * @return <code>value</code> that was passed in
     */
    public double broadcast(double value, int root) {
        return value;
    }
    
    /**
     * Does a broadcast of <code>value</code> to the node itself.
     *
     * @param value The Serializable object to broadcast to all nodes.
     * @param root The root vnode to gather the Serializable objects onto.  This is not used for serial communication and is ignored.
     * @return <code>value</code> that was passed in
     */
    public java.io.Serializable broadcast(java.io.Serializable value, int root) {
        return value;
    }
    
    /**
     * Simply returns an array of size one of <code>myDouble</code>.
     *
     * @param myDouble the double to gather
     * @return <code>myDouble</code>
     */
    public double[] gatherAll(double myDouble) {
        double[] temp = {myDouble};
        return temp;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     *
     * @param myElements the double array to gather
     * @return <code>myElements</code>
     */
    public double[] gatherAll(double[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     *
     * @param myElements the serialized object array to gather
     * @return <code>myElements</code>
     */
    public java.io.Serializable[] gatherAll(java.io.Serializable[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     *
     * @param myElements the int array to gather
     * @return <code>myElements</code>
     */
    public int[] gatherAll(int[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns an array of size one of <code>myInt</code>.
     *
     * @param myInt the int to gather
     * @return <code>myInt</code>
     */
    public int[] gatherAll(int myInt) {
        int[] temp = {myInt};
        return temp;
    }
    
    /**
     * Accessor to determine if <code>this</code> is in serial mode.  Always returns true.
     * @return true
     */
    public boolean isSerial() {
        return true;
    }
    
    /**
     * Threads are not supported yet and thus this value returned will always be 1.
     * @return 1
     */
    public int getNumMyThreads() {
        return 1;
    }
    
    /**
     * Since threads are not supported this value will always be one.
     * @return 1
     */
    public int getNumThreads() {
        return 1;
    }
    
    /**
     * Since this is a serial implementation of <code>comm</code> the value returned will always be 1.
     * @return 1
     */
    public int getNumVnodes() {
        return 1;
    }
    
    /**
     * Always returns 0 since threads are not yet implemented.
     * @return 0
     */
    public int getThreadId() {
        return 0;
    }
    
    /**
     * Since this is a serial implementation of <code>Comm</code> it always returns 0.
     * @return 0
     */
    public int getVnodeId() {
        return 0;
    }
    
    /**
     * Simply returns <code>partialMaxs</code>.
     * @param partialMaxs the partial maximums
     * @return <code>partialMaxs</code>
     */
    public double[] maxAll(double[] partialMaxs) {
        return partialMaxs;
    }
    
    /**
     * Simply returns <code>partialMaxs</code>.
     * @param partialMaxs the partial maximums
     * @return <code>partialMaxs</code>
     */
    public int[] maxAll(int[] partialMaxs) {
        return partialMaxs;
    }
    
    /**
     * Simply returns <code>partialMins</code>.
     * @param partialMins the partial minimums
     * @return <code>partialMins</code>
     */
    public int[] minAll(int[] partialMins) {
        return partialMins;
    }
    
    /**
     * Simply returns <code>partialMins</code>.
     * @param partialMins the partial minimums
     * @return <code>partialMins</code>
     */
    public double[] minAll(double[] partialMins) {
        return partialMins;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     * @param myElements the values to be summmed
     * @return <code>myElements</code>
     */
    public double[] scanSums(double[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     * @param myElements the values to be summmed
     * @return <code>myElements</code>
     */
    public int[] scanSums(int[] myElements) {
        return myElements;
    }
    
    /**
     * Not implemeted.
     * @param newVnodeID new vnode ID for <code>this</code vnode
     */
    public void setMyVnodeID(int newVnodeID) {
        // empty
    }
    
    /**
     * Not implemented.
     * @param newThreadID new thread ID for <code>this</code thread
     */
    public void setThreadID(int newThreadID) {
        // empty
    }
    
    /**
     * Simply returns <code>partialSums</code>.
     *
     * @param partialSums the partial sums
     * @return <code>partialSums</code>
     */
    public int[] sumAll(int[] partialSums) {
        return partialSums;
    }
    
    /**
     * Simply returns <code>partialSums</code>.
     *
     * @param partialSums the partial sums
     * @return <code>partialSums</code>
     */
    public double[] sumAll(double[] partialSums) {
        return partialSums;
    }
    
    /**
     * Not yet implemented.
     */
    public void threadBarrier() {
        // empty
    }
    
    //public Distributor createDistributor() {
    //}
    
    public Directory createDirectory(ElementSpace elementSpace) {
        return new SerialDirectory(elementSpace);
    }
    
    public Distributor createDistributor() {
    }
    
    public int[][] gatherAll2dArray(int[] myElements) {
    }
    
    public void send(Serializable[] exportObject, int destinationVnode) {
    }
    
    public void send(int[] exportObject, int destinationVnode) {
    }
    
    public void setupReceives(int numReceives) {
    }
    
    public Serializable receive(int senderId) {
    }
    
    public int[] scatterIntArray(int[] in) {
    }
    
    public int[] scatter2dArray(int[][] in) {
    }
    
    public int[][] gather(int[] in) {
    }
    
    public void send(java.io.Serializable exportObject, int destinationVnode) {
    }
    
}
