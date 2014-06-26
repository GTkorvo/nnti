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
 * Comm is the interface for the all the Jpetra communication
 * classes.
 *
 * Terminology:  A vnode is considered to be a virtual node (a memory image) and not a physical machine.
 * A physical node is considered to be a physical machine, and may have many vnodes on it.
 * The term node has no meaning, and should be qualified by virtual or physical.
 *
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @author  Jason Cross
 */
public interface Comm {
    
    /**
     * Causes each vnode in the communicator to wait until all
     * vnodes have arrived.
     * No-op for a serial communicator.
     */
    public void barrier();
    
    /**
     * Broadcasts any Serializable object from the <code>root</code>
     * vnode to all other vnodes.  Requires the return object to
     * be cast into a more useful object.
     *
     * @param value object that is broadcast by the root vnode
     * @param root root vnode ID, most likely 0
     * @return all vnodes recieve this object
     */
    public Serializable broadcast(Serializable value, int root);
    
    /**
     * Broadcasts an int from the <code>root</code>
     * vnode to all other vnodes.
     *
     * @param value int that is broadcast by the root vnode
     * @param root root vnode ID, most likely 0
     *
     * @return all vnodes recieve this int
     */
    public int broadcast(int value, int root);
    
    /**
     * Broadcasts an double from the <code>root</code>
     * vnode to all other vnodes.
     *
     * @param value double that is broadcast by the root vnode
     * @param root root vnode ID, most likely 0
     *
     * @return all vnodes recieve this int
     */
    public double broadcast(double value, int root);
    
    /**
     * Takes a list of Serilizable objects from all vnodes in the communicator and
     * creates an ordered contiguous list of those objects in each vnode.
     * Requires the return object to be cast into a more useful object.
     *
     * @param myElements     in on entry; the list of values to be sent to all vnodes
     *
     * @return allElements   out on exit; the list of values from all vnodes
     */
    public Serializable[] gatherAll(Serializable [] myElements);
    
    /**
     * Takes a list of values from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myElements     in on entry; the list of values to be sent to all vnodes
     *
     * @return allElements   out on exit; the list of values from all vnodes
     */
    public double[] gatherAll(double [] myElements);
    
    /**
     * Takes a list of values from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myElements    in on entry; the list of values to be sent to all vnodes
     *
     * @return              out on exit; the list of values from all vnodes
     */
    public int[] gatherAll(int [] myElements);
    
    public int[][] gatherAll2dArray(int[] myElements);
    
    /**
     * Takes a value from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myInt     in on entry; the list of values to be sent to all vnodes
     *
     * @return          out on exit; the list of values from all vnodes
     */
    public int[] gatherAll(int myInt);
    
    /**
     * Takes a value from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myDouble     in on entry; the list of values to be sent to all vnodes
     *
     * @return          out on exit; the list of values from all vnodes
     */
    public double[] gatherAll(double myDouble);
    
    /**
     * Takes a list of values from all vnodes in the communicator, computes
     * the sum, and returns the sum to all vnodes.
     *
     * @param partialSums   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed across all vnodes
     */
    public double[] sumAll(double [] partialSums);
    
    /**
     * Takes a list of values from all vnodes in the communicator, computes
     * the sum, and returns the sum to all vnodes.
     *
     * @param partialSums   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed
     *                      across all vnodes
     */
    public int[] sumAll(int [] partialSums);
    
    /**
     * Takes a list of values form all vnodes in the communicator, computes
     * the max, and returns the max to all vnodes.
     *
     * @param partialMaxs   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed across all vnodes
     */
    public double[] maxAll(double [] partialMaxs);
    
    /**
     * Takes a list of values form all vnodes in the communicator, computes
     * the max, and returns the max to all vnodes.
     *
     * @param partialMaxs   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed across all vnodes
     */
    public int[] maxAll(int [] partialMaxs);
    
    public Serializable maxAll(Serializable partialMaxs);
    
    /**
     * Takes a list of values form all vnodes in the communicator, computes
     * the min, and returns the min to all vnodes.
     *
     * @param partialMins   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed across all vnodes
     */
    public double[] minAll(double [] partialMins);
    
    public Serializable minAll(Serializable partialMins);
    
    /**
     * Takes a list of values form all vnodes in the communicator, computes
     * the min, and returns the min to all vnodes.
     *
     * @param partialMins   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed across all vnodes
     */
    public int[] minAll(int [] partialMins);
    
    /**
     * Takes a list of values from all vnodes in the communicator, computes
     * the scan sum, and returns it to all vnodes such that vnode i contains
     * the sum of values from vnodes up to and including vnode i.
     *
     * @param myElements    in on entry; the values to be summed across all vnodes
     *
     * @return              out on exit; the list of values summed
     *                      across vnodes 0 through i
     */
    public double[] scanSums(double [] myElements);
    
    /**
     * Takes a list of values from all vnodes in the communicator, computes
     * the scan sum, and returns it to all vnodes such that vnode i contains
     * the sum of values from vnodes up to and including vnode i.
     *
     * @param myElements    in on entry; the values to be summed across all vnodes
     *
     * @return              out on exit; the list of values summed
     *                      across vnodes 0 through i
     */
    public int[] scanSums(int [] myElements);
    
    /**
     * Accessor for the number of vnodes in the commmunicator.
     *
     * @return number of vnodes
     */
    public int getNumVnodes();
    
    /**
     * Accessor for the rank of the calling vnode.
     *
     * @return the rank of the calling vnode in MPI (CCJ); 0 in serial mode
     */
    public int getVnodeId();
    
    /**
     * Sets the vnode ID for the calling vnode.
     * Can be used to facilitate threaded programming acroass an MPI
     * application by associating several MPI vnodes with a single
     * vnode. By default, each MPI vnode is associated with a single
     * vnode wit the same ID.
     *
     * @param newVnodeID new vnode ID
     */
    public void setMyVnodeID(int newVnodeID);
    
    /**
     * Sends an int arry to a single specified vnode.
     * <b>Note<b>: this is NOT a blocking operation.
     *
     * @param exportObject the int array to be sent
     * @param destinationVnode the vnode ID of the receiving vnode
     */
    public void send(int[] exportObject, int destinationVnode);
    
    /**
     * Sends a double arry to a single specified vnode.
     * <b>Note<b>: this is NOT a blocking operation.
     *
     * @param exportObject the double array to be sent
     * @param destinationVnode the vnode ID of the receiving vnode
     */
    public void send(double[] exportObject, int destinationVnode);
    
    /**
     * Sends a Serializable arry to a single specified vnode.
     * <b>Note<b>: this is NOT a blocking operation.
     *
     * @param exportObject the Serializable array to be sent
     * @param destinationVnode the vnode ID of the receiving vnode
     */
    public void send(Serializable exportObject, int destinationVnode);
    
    /**
     * <code>getReceives</code> does all the work receiving all expected messages at once and then returns them.
     * <b>Note<b>: this IS a blocking operation.
     *
     * return the objects received
     */
    public Serializable receive(int senderId);
    
    public int[] scatter2dArray(int[][] in);
    
    public int[] scatterIntArray(int[] in);
    
    public int[][] gather(int[] in);
    
    /**
     * Accessor to see if the communicator is on a serial or parallel machine.
     *
     * @return <code>true</code> in serial mode, otherwise <code>false</code>
     */
    public boolean isSerial();
    
    public Directory createDirectory(VectorSpace vectorSpace);
    
    public Distributor createDistributor();
}