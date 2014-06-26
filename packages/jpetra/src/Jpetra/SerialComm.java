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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

package Jpetra;

import java.io.Serializable;

/**
 * <code>SerialComm</code> is a simple class that implements the <code>Comm</code>
 * interface for serial communication.  Most methods either simply return the
 * parameters passed to them or return a hard coded value.
 *
 * @author  Jason Cross
 */
public class SerialComm extends JpetraObject implements Comm {
    
    /**
     * Creates a new instance of SerialComm.
     */
    public SerialComm() {
        // empty
    }
    
    /**
     * No-op for a serial communicator.
     */
    public void barrier() {
        // empty
    }
    
    /**
     * Does a broadcast of <code>value</code> to the node itself.
     *
     * @param value The int to broadcast to all nodes.
     * @param root The root vnode to gather the ints onto.  This is not used for serial communication and is ignored.
     *
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
     *
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
     *
     * @return <code>value</code> that was passed in
     */
    public Serializable broadcast(java.io.Serializable value, int root) {
        return value;
    }
    
    /**
     * Simply returns an array of size one of <code>myDouble</code>.
     *
     * @param myDouble the double to gather
     *
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
     *
     * @return <code>myElements</code>
     */
    public double[] gatherAll(double[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     *
     * @param myElements the serialized object array to gather
     *
     * @return <code>myElements</code>
     */
    public java.io.Serializable[] gatherAll(java.io.Serializable[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns <code>myElements</code>.
     *
     * @param myElements the int array to gather
     *
     * @return <code>myElements</code>
     */
    public int[] gatherAll(int[] myElements) {
        return myElements;
    }
    
    /**
     * Simply returns an array of size one of <code>myInt</code>.
     *
     * @param myInt the int to gather
     *
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
     * Since this is a serial implementation of <code>comm</code> the value returned will always be 1.
     * @return 1
     */
    public int getNumVnodes() {
        return 1;
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
    public Serializable minAll(Serializable partialMins) {
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
     * No-op for a serial communicator.
     *
     * @param newVnodeID new vnode ID for <code>this</code vnode
     */
    public void setMyVnodeID(int newVnodeID) {
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
    
    public int[][] gatherAll2dArray(int[] myElements) {
        int[][] out = new int[1][];
        out[0] = myElements;
        return out;
    }
    
    /*
     * No-op for a serial communicator.
     */
    public void send(Serializable[] exportObject, int destinationVnode) {
        // empty
    }
    
    /**
     * No-op for a serial communicator.
     */
    public void send(int[] exportObject, int destinationVnode) {
        // empty
    }
    
    /**
     * No-op for a serial communicator.
     */
    public Serializable receive(int senderId) {
        return null;
    }
    
    public int[] scatterIntArray(int[] in) {
        return new int[]{in[0]};
    }
    
    public int[] scatter2dArray(int[][] in) {
        return in[0];
    }
    
    public int[][] gather(int[] in) {
        int[][] out = new int[1][];
        out[0] = in;
        return out;
    }
    
    /**
     * No-op for a serial communicator.
     */
    public void send(Serializable exportObject, int destinationVnode) {
    }
    
    /**
     * No-op for a serial communicator.
     */
    public void send(double[] exportObject, int destinationVnode) {
    }
    
    public Distributor createDistributor() {
        return new SerialDistributor();
    }
    
    public Directory createDirectory(VectorSpace vectorSpace) {
        return new BasicDirectory(vectorSpace);
    }
    
    public Serializable maxAll(Serializable partialMaxs) {
        return partialMaxs;
    }
    
    public double[] minAll(double[] partialMins) {
        return partialMins;
    }
    
}
