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

package Jpetra.CcjSupport;

import CCJ.*;
import java.io.Serializable;

/*
 * CCJCom.java
 *
 * Created on Tue June 17 11:59:00 CDT 2003
 */

/**
 * <code>CcjLink</code> provides the wrappers and data representation objects
 * that <code>CcjComm</code> uses to interact with CCJ.
 *
 * @author  Jason Cross
 */

public class CcjLink extends ColMember {
    
    /**
     * number of vnodes in the <code>group</code>
     */
    private int numVnodes;
    
    /**
     * vnode group; all vnodes belong to the same group
     */
    private ColGroup group;
    
    /**
     * root vnode
     */
    private ColGroupMaster groupMaster;
    
    /**
     * this rank
     */
    private int rank;
    
    /**
     * Contacts the root vnode and joins the vnode
     * <code>group</code> unless <code>this</code> is the root vnode, then it
     * listens for connections.
     */
    public CcjLink(ColGroupMaster groupMaster) throws CCJException {
        super();
        this.groupMaster = groupMaster;
        
        this.numVnodes = groupMaster.getNumberOfCpus();
        
        // notice that the group name is a string
        this.groupMaster.addMember("myGroup", this);
        this.group = groupMaster.getGroup("myGroup", this.numVnodes);
        
        this.rank = group.getRank(this);
        
        
        // calls the run method
        // not used by CcjLink
        // begin();
    }
    
    /**
     * Accessor for <code>numVnodes</code> provided by CCJ.
     *
     * @return <code>numVnodes</code>
     */
    public int getNumVnodes() {
        return this.numVnodes;
    }
    
    /**
     * Accessor for <code>rank</code> provided by CCJ.
     *
     * @return <code>rank</code>
     */
    public int getRank() {
        return this.rank;
    }
    
    public ColGroup getGroup() {
        return this.group;
    }
    
    /**
     * Wrapper to CCJ <code>barrier</code>.  Causes each vnode in <code>group</group> to wait
     * until all vnodees are ready.
     */
    public void barrier() {
        try {
            barrier(this.group);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ barrier: " + e);
        }
    }
    
    /**
     * Wrapper to CCJ <code>broadcast</code>.  Broadcasts <code>value</code> from the <code>root</code>
     * vnode to all other vnodees in <code>group</code>.
     */
    public Serializable broadcast(Serializable value, int root) {
        try {
            value = (Serializable) broadcast(group, value, root);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ broadcast: " + e);
        }
        
        return value;
    }
    
    /**
     * Wrapper to CCJ <code>allGather</code>.  Broadcasts <code>value</code> from the <code>root</code>
     * vnode to all other vnodees in <code>group</code>.
     */
    public Serializable[] gatherAll(Serializable [] myElements) {
        //used internally by CCJ to handle the <code>gatherAll</code> of arrays
        CcjGatherSerializableArray allElements = new CcjGatherSerializableArray(group.size());
        
        try {
            allGather(group, allElements, myElements);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ allGather: " + e);
        }
        
        return allElements.getAllElements();
    }
    
    /**
     * Wrapper to CCJ <code>allGather</code>.  Broadcasts <code>value</code> from the <code>root</code>
     * vnode to all other vnodees in <code>group</code>.
     */
    public int[] gatherAll(int [] myElements) {
        //used internally by CCJ to handle the <code>gatherAll</code> of arrays
        CcjGatherIntArray allElements = new CcjGatherIntArray(group.size());
        
        try {
            allGather(group, allElements, myElements);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ allGather: " + e);
        }
        
        return allElements.getAllElements();
    }
    
    public int[][] gatherAll2dArray(int[] myElements) {
        //used internally by CCJ to handle the <code>gatherAll</code> of arrays
        CcjGatherIntArray allElements = new CcjGatherIntArray(group.size());
        
        try {
            allGather(group, allElements, myElements);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ allGather: " + e);
        }
        
        return allElements.getAllElements2dArray();
    }
    
    /**
     * Wrapper to CCJ <code>allGather</code>.
     */
    public double[] gatherAll(double [] myElements) {
        //used internally by CCJ to handle the <code>gatherAll</code> of arrays
        CcjGatherDoubleArray allElements = new CcjGatherDoubleArray(group.size());
        
        try {
            allGather(group, allElements, myElements);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ allGather: " + e);
        }
        
        return allElements.returnAllElements();
    }
    
    /**
     * Wrapper to CCJ <code>allReduce</code>.
     */
    public double[] sumAll(double[] partialSums) {
        //used internally by CCJ to handle the <code>allReduce</code> of arrays
        CcjReduceAllSumDoubleArray doublesAdder = new CcjReduceAllSumDoubleArray();
        
        double[] toReturn = null;
        try {
            toReturn = (double[]) allReduce(group, partialSums, doublesAdder);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ doubleSumAll: " + e);
        }
        
        return toReturn;
    }
    
    /**
     * Wrapper to CCJ <code>allReduce</code>.
     */
    public int[] sumAll(int[] partialSums) {
        //used internally by CCJ to handle the <code>allReduce</code> of arrays
        CcjReduceAllSumIntArray intsAdder = new CcjReduceAllSumIntArray();
        
        int[] toReturn=null;
        try {
            toReturn = (int[]) allReduce(group, partialSums, intsAdder);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ intSumAll: " + e);
        }
        
        return toReturn;
    }
    
    /**
     * Wrapper to CCJ <code>allReduce</code>.
     */
    public int[] maxAll(int [] partialMaxs) {
        CcjReduceIntMaxArray maxInts = new CcjReduceIntMaxArray();
        
        int[] toReturn=null;
        try {
            toReturn = (int[]) allReduce(group, partialMaxs, maxInts);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ intMaxAll: " + e);
        }
        
        return toReturn;
    }
    
    /**
     * Wrapper to CCJ <code>allReduce</code>.
     */
    public double[] maxAll(double [] partialMaxs) {
        CcjReduceDoubleMaxArray maxDoubles = new CcjReduceDoubleMaxArray();
        
        double[] toReturn=null;
        try {
            toReturn = (double[]) allReduce(group, partialMaxs, maxDoubles);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ doubleMaxAll: " + e);
        }
        
        return toReturn;
    }
    
    /**
     * Wrapper to CCJ <code>allReduce</code>.
     */
    public double[] minAll(double [] partialMins) {
        CcjReduceDoubleMinArray minDoubles = new CcjReduceDoubleMinArray();
        
        double[] toReturn=null;
        try {
            toReturn = (double[]) allReduce(group, partialMins, minDoubles);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ doubleMinAll: " + e);
        }
        
        return toReturn;
    }
    
    /**
     * Wrapper to CCJ <code>allReduce</code>.
     */
    public int[] minAll(int [] partialMins) {
        CcjReduceIntMinArray minInts = new CcjReduceIntMinArray();
        
        int[] toReturn=null;
        try {
            toReturn = (int[]) allReduce(group, partialMins, minInts);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ intMinAll: " + e);
        }
        
        return toReturn;
    }
    
    /**
     * Wrapper to CCJ <code>allGather</code> then calls <code>CcjGatherDoubleArray.scanSums()</code>.
     */
    public double[] scanSums(double [] myElements) {
        //used internally by CCJ to handle the <code>gatherAll</code> of arrays
        CcjGatherDoubleArray globalSums = new CcjGatherDoubleArray(group.size());
        
        try {
            allGather(group, globalSums, myElements);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ scanSums: " + e);
        }
        
        return globalSums.scanSums(rank);
    }
    
    /**
     * Wrapper to CCJ <code>allGather</code> then calls <code>CcjGatherIntArray.scanSums()</code>.
     */
    public int[] scanSums(int [] myElements) {
        //used internally by CCJ to handle the <code>gatherAll</code> of arrays
        CcjGatherIntArray globalSums = new CcjGatherIntArray(group.size());
        
        try {
            allGather(group, globalSums, myElements);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ scanSums: " + e);
        }
        
        return globalSums.scanSums(rank);
    }
    
    public int[] scatter2dArray(int[][] in) {
        CcjScatter2DArray scatter = new CcjScatter2DArray(in);
        
        int[] out = null;
        try {
            out = (int[]) scatter(group, scatter, 0);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ scatter2dArray: " + e);
        }
        
        return out;
    }
    
    public int[] scatterIntArray(int[] in) {
        CcjScatterIntArray scatter = new CcjScatterIntArray(in);
        
        int[] out = null;
        try {
            out = (int[]) scatter(group, scatter, 0);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ scatter2dArray: " + e);
        }
        
        return out;
    }
    
    public int[][] gather(int[] in) {
        CcjGatherIntArray rec = new CcjGatherIntArray(group.size());
        
        try {
            rec = (CcjGatherIntArray) gather(group, rec, in, 0);
            
            if (group.getRank(this) == 0) {
                return rec.getAllElements2dArray();
            }
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ gather: " + e);
        }
        return null;
    }
    
    /**
     * Wrapper to CCJ <code>setupRecords</code> which tells CCJ how many messages to expect to receive.
     *
     * @param numReceives the number of messages that CCJ should expect to receive
     */
    public void setupReceives(int numReceives) {
        setupRecords(numReceives);  // tell CCJ how many messages to wait for
    }
    
    /**
     * Wrapper to CCJ <code>getRecords</code>, <code>receive</code>, and <code>endRecords</code>.
     * <code>setupRecords</code> must be called before <code>getReceives</code>.  <code>getReceives</code>
     * does all the work receiving all expected messages at once and then returns them.  <b>Note</b>: this
     * IS a blocking operation.
     *
     * return the objects that CCJ received
     */
    /*public Serializable[] getReceives() {
        int[] records = getRecords();  // get the message IDs that CCJ has received
        Serializable[] results = new Serializable[records.length];
     
        try {
            for(int i=0; i < records.length; i++) {
                System.out.println("Going to try to receive message " + records[i]);
                results[i] = receive(group, records[i]);
                System.out.println("Message " + records[i] + " has been received.");
            }
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ receive: " + e);
        }
     
        endRecords();  // tell CCJ we're done with messages and don't expect to receive anymore
        return results;
    }*/
    public Serializable receive(int senderId) {
        try {
            return receive(group, senderId);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ receive: " + e);
        }
        return null;  // we should never get this far unless their is an error
    }
    
    /**
     * Sends an int arry to a single specified vnode.  Wrapper to CCJ <code>send_async</code>
     * <b>Note</b>: this is NOT a blocking operation.
     *
     * @param exportObject the int arrray to be sent
     * @param destinationVnode the vnode ID of the receiving vnode
     */
    public void send(int[] exportObject, int destinationVnode) {
        try {
            send_async(group, exportObject, destinationVnode);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ send: " + e);
        }
    }
    
    /**
     * Sends a double arry to a single specified vnode.  Wrapper to CCJ <code>send_async</code>
     * <b>Note</b>: this is NOT a blocking operation.
     *
     * @param exportObject the double arrray to be sent
     * @param destinationVnode the vnode ID of the receiving vnode
     */
    public void send(double[] exportObject, int destinationVnode) {
        try {
            send_async(group, exportObject, destinationVnode);
        }
        catch (CCJException e) {
            System.err.println("Error in CCJ send: " + e);
        }
    }
    
    /**
     * not used but required to extend CCJ <code>ColMember</code>
     */
    public void run() {
        // empty
    }
}
