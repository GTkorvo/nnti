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

import java.io.Serializable;
/*
 * SerialComm.java
 *
 * Created on May 29, 2001, 1:32 PM
 */


/**
 * SerialComm is the implementatin of the Comm interface
 * for machines with a single vnode.
 *
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @version 
 */
public class SerialComm extends JpetraObject implements Comm {

    private int myRank;
    private int myVnode;
    private int myThread;
    private int numVnodes;
    private boolean isSerial;
    
    /** Creates new <code>SerialComm</code> */
    public SerialComm() {
        myVnode = 0;
        numVnodes = 1;
    }
    
    /**
     * Creates a new <code>SerialComm</code> with the same attributes as the argument.
     *
     * @param comm  the <code>SerialComm</code> of which to make a copy
     */
    public SerialComm(SerialComm comm) {
        myRank = comm.myRank;
        myVnode = comm.myVnode;
        myThread = comm.myThread;
        numVnodes = comm.numVnodes;
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
    
    public void barrier() {}
    
    public void threadBarrier() {}
    
    public Serializable broadcast(Serializable value, int root) {
        return value;
    }
    
    public double[] broadcast(double[] elements,int root) {
        return elements;
    }
    
    public int[] broadcast(int[] elements,int root) {
        return elements;
    }
    
    public int broadcast(int value, int root) {
        return value;
    }
    
    public double broadcast(double value, int root) {
        return value;
    }
    
    public double[] gatherAll(double[] myElements) {
        return myElements;
    }
    
    public int[] gatherAll(int[] myElements) {
        return myElements;
    }
    
    public int[] gatherAll(int myInt) {
        return new int[]{myInt};   
    }
    
    public double[] gatherAll(double myDouble) {
        return new double[]{myDouble};   
    }

    public Serializable[] gatherAll(Serializable [] myElements) {
        return myElements;
    }
    
    public double[] sumAll( double [] partialSums) {
        return partialSums;
    }
    
    public int[] sumAll(int[] partialSums) {
        return partialSums;
    }
    
    public double[] maxAll(double[] partialMaxs) {
        return partialMaxs;
    }
    
    public int[] maxAll(int[] partialMaxs) {
        return partialMaxs;
    }
    
    public double[] minAll(double[] partialMins) {
        return partialMins;
    }
    
    public int[] minAll(int[] partialMins) {
        return partialMins;
    }
    
    public double[] scanSums(double[] myElements) {
        return myElements;
    }
    
    public int[] scanSums(int[] myElements) {
        return myElements;
    }
    
    public int getVnodeID() {
        return 0;
    }
    
    public int getThreadID() {
        return myThread;
    }
    
    public int setThreadID(int aThread) {
        myThread = aThread;
        return 0;
    }
    
    public int getNumThreads() {
        return 1;
    }
    
    public int getNumMyThreads() {
        return 1;
    }    
    
    public int setMyVnodeID(int newVnodeID) {
        myVnode = newVnodeID;
        return 0;
    }
    
        public boolean getIsSerial() {
        return isSerial;
    }      
    
    public int getMyVnode() {
        return myVnode;
    }
    
    /**
     * Accessor for the number of vnodes in the communicator.
     */
    public int getNumVnodes() {
        return 1;
    }
    
    public Distributor createDistributor() {
    	return new SerialDistributor(this);
    }
    
}
