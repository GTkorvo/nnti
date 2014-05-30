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

import CCJ.*;
import java.io.Serializable;

import Jpetra.CcjSupport.*;

/**
 * <code>CcjComm</code> is the implementatin of the <code>Comm</code> interface
 * for vnodes running in a parallel environment.
 * It utilizes the Collective Communication (CCJ) API
 * for MPI like message passing.
 * <p>
 * To use <code>CcjComm</code> first an instance of CcjComm needs to be created.
 * A the path and file name of a ccjhosts file needs to be passed to the <code>CcjComm</code> contructor.
 * The format of this file is one IP or host name and port per line, seperated
 * by a colon (:).  Note that the port number is optional.
 * <p>Example ccjhosts file:
 * <pre>
 * // comments begining with /
 * # and # are allowed; blank lines are also ok
 *
 * //root vnode
 * 192.168.1.8
 * //slave vnodes
 * 192.168.1.8:1
 * 192.168.1.8:24
 * sun1
 * sun2
 * //do not put comments on a host line: 192.168.1.5 // root  <== this comment will break things!
 * </pre>
 *
 * Each vnode needs to be able to access the ccjhosts file.  This could entail copying the
 * file to the same path on each physical node, or a shared file system could be used.
 *<p>
 *
 * Adapted from <code>SerialComm</code> by Mike Heroux and
 * Michael William Boldt.
 *
 * @author  Jason Cross
 */
public class CcjComm extends JpetraObject implements Comm {
    /**
     * this vnode's rank/ID
     */
    private int myVnode;
    
    /**
     * total number of vnodes
     */
    private int numVnodes;
    
    /**
     * this vnode's instance of the wrapper object between Jpetra and CCJ
     */
    private CcjLink myCcjLink;
    
    /**
     * true if there is only one vnode
     */
    private boolean isSerial;
    
    /** Creates new <code>CcjComm</code> */
    public CcjComm(String filePath) {
        
        //create the CcjLink object used for wrapping CcjComm calls
        try {
            ColGroupMaster groupMaster = new ColGroupMaster(filePath);
            this.myCcjLink = new CcjLink(groupMaster);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ Setup: " + e);
        }
        
        //set fields
        this.myVnode = this.myCcjLink.getRank();
        this.numVnodes = this.myCcjLink.getNumVnodes();
        
        //threads are not implemented, so set to 1
        if (this.numVnodes == 1) {
            this.isSerial = true;
        }
        
    }
    
    /**
     * Accessor for <code>isSerial</code>.
     *
     * @return <code>isSerial</code>
     */
    public boolean isSerial() {
        return isSerial;
    }
    
    /**
     * Accessor for <code>myVnode</code>.
     *
     * @return <code>myVnode</code>
     */
    public int getVnodeId() {
        return myVnode;
    }
    
    /**
     * Accessor for <code>numVnodes</code>.
     *
     * @return <code>numVnodes</code>
     */
    public int getNumVnodes() {
        return numVnodes;
    }
    
    /**
     * Makes each CCJ vnode wait until all CCJ vnodes are ready.
     */
    public void barrier() {
        this.myCcjLink.barrier();
    }
    
    /**
     * Makes each thread in this CCJ vnode wait until all threads in this vnode are ready.
     * Not currently implemented.
     */
    public void threadBarrier() {}
    
    /**
     * Wrapper to CCJ broadcast.
     */
    public Serializable broadcast(Serializable value, int root) {
        return this.myCcjLink.broadcast(value, root);
    }
    
    /**
     * Wrapper to CCJ broadcast.
     */
    public int broadcast(int value, int root) {
        Integer valueInteger = (Integer) this.myCcjLink.broadcast(new Integer(value), root);
        return valueInteger.intValue();
    }
    
    /**
     * Wrapper to CCJ broadcast.
     */
    public double broadcast(double value, int root) {
        Double valueDouble = (Double) this.myCcjLink.broadcast(new Double(value), root);
        return valueDouble.doubleValue();
    }
    
    /**
     * Wrapper to CCJ gatherAll.
     */
    public Serializable[] gatherAll(Serializable [] myElements) {
        return this.myCcjLink.gatherAll(myElements);
    }
    
    /**
     * Wrapper to CCJ gatherAll.
     */
    public int[] gatherAll(int [] myElements) {
        return this.myCcjLink.gatherAll(myElements);
    }
    
    public int[][] gatherAll2dArray(int[] myElements) {
        return this.myCcjLink.gatherAll2dArray(myElements);
    }
    
    /**
     * Wrapper to CCJ gatherAll.
     */
    public double[] gatherAll(double [] myElements) {
        return this.myCcjLink.gatherAll(myElements);
    }
    
    /**
     * Wrapper to CCJ gatherAll.
     */
    public int[] gatherAll(int myInt) {
        int[] myIntArray = new int[]{myInt};
        return this.myCcjLink.gatherAll(myIntArray);
    }
    
    /**
     * Wrapper to CCJ gatherAll.
     */
    public double[] gatherAll(double myDouble) {
        double[] myDoubleArray = new double[]{myDouble};
        return this.myCcjLink.gatherAll(myDoubleArray);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>sumAll</code>.
     */
    public double[] sumAll(double [] partialSums) {
        return this.myCcjLink.sumAll(partialSums);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>sumAll</code>.
     */
    public int[] sumAll(int [] partialSums) {
        return this.myCcjLink.sumAll(partialSums);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>maxAll</code>.
     */
    public double[] maxAll(double [] partialMaxs) {
        return this.myCcjLink.maxAll(partialMaxs);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>maxAll</code>.
     */
    public int[] maxAll(int [] partialMaxs) {
        return this.myCcjLink.maxAll(partialMaxs);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>minAll</code>.
     */
    public Serializable minAll(Serializable partialMins) {
        return this.myCcjLink.minAll(partialMins);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>minAll</code>.
     */
    public double[] minAll(double[] partialMins) {
        return this.myCcjLink.minAll(partialMins);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>minAll</code>.
     */
    public int[] minAll(int [] partialMins) {
        return this.myCcjLink.minAll(partialMins);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>scanSums</code>.
     */
    public double[] scanSums(double [] myElements) {
        return this.myCcjLink.scanSums(myElements);
    }
    
    /**
     * Wrapper to <code>CcjLink</code> <code>scanSums</code>.
     */
    public int[] scanSums(int [] myElements) {
        return this.myCcjLink.scanSums(myElements);
    }
    
    /*public send(int[] exportObj, int destinationVnode) {
        this.myCcjLink.send(exportObj, int destinationVnode);
    }
     
    public send(int[] exportOb, int destinationVnode) {
        this.myCcjLink.send(exportObj, int destinationVnode);
    }
     
    public int[] recieve(int sourceVnode) {
        return this.myCcjLink.send(int sourceVnode);
    }
     
    public double[] recieve(int sourceVnode) {
        return this.myCcjLink.send(int sourceVnode);
    }*/
    
    /**
     * Accessor for <code>myVnode</code>.
     *
     * @return <code>myVnode</code>
     */
    public int getVnodeID() {
        return myVnode;
    }
    
    /**
     * Sets the vnode ID.  The ID is also the rank, do not change this
     * unless you know exaclty what you are doing.
     *
     * @param newVnodeID new proccess ID/rank
     */
    public void setMyVnodeID(int newVnodeID) {
        myVnode = newVnodeID;
    }
    
    public Distributor createDistributor() {
        return new CcjDistributor();
    }
    
    public CcjLink getCcjLink() {
        return this.myCcjLink;
    }
    
    public String getGroupName() {
        return "";
    }
    
    public Directory createDirectory(VectorSpace vectorSpace) {
        return new BasicDirectory(vectorSpace);
    }
    
    public void send(double[] exportObject, int destinationVnode) {
        this.myCcjLink.send(exportObject, destinationVnode);
    }
    
    public void send(int[] exportObject, int destinationVnode) {
        this.myCcjLink.send(exportObject, destinationVnode);
    }
    
    public void send(Serializable exportObject, int destinationVnode) {
        this.myCcjLink.send(exportObject, destinationVnode);
    }
    
    public void setupReceives(int numReceives) {
        this.myCcjLink.setupReceives(numReceives);
    }
    
    public Serializable receive(int senderId) {
        return this.myCcjLink.receive(senderId);
    }
    
    public int[] scatter2dArray(int[][] in) {
        return this.myCcjLink.scatter2dArray(in);
    }
    
    public int[] scatterIntArray(int[] in) {
        return this.myCcjLink.scatterIntArray(in);
    }
    
    public int[][] gather(int[] in) {
        return this.myCcjLink.gather(in);
    }
    
    public Serializable[][] gather(Serializable[] in) {
        return this.myCcjLink.gather(in);
    }
    
    public Serializable maxAll(Serializable in) {
        return this.myCcjLink.maxAll(in);
    }
    
    public String getMyHostName() {
        return this.myCcjLink.getMyHostName();
    }
    
}
