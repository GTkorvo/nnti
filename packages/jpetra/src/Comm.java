package Jpetra;

import java.io.Serializable;
/*
 * Comm.java
 *
 * Created on May 29, 2001, 1:25 PM
 */


/**
 * Comm is the interface for the all the Jpetra communication
 * classes.
 *
 * NOTE:  A vnode is considered to be a memory image and not a physical machine.
 *
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @author  Jason Cross
 * @version 
 */
public interface Comm {
    
    /**
     * Causes each vnode in the communicator to wait until all
     * vnodes have arrived.
     * No-op for a serial communicator. 
     */
    public void barrier();

    /**
     * Causes each thread in a given vnode on each machine
     * to wait until all threads in that vnode have arrived.
     * No-op for a serial communicator.
     */        
    public void threadBarrier();
    
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
     * @return all vnodes recieve this int
     */
    public int broadcast(int value, int root);
   
   /**
     * Broadcasts an double from the <code>root</code>
     * vnode to all other vnodes.
     * 
     * @param value double that is broadcast by the root vnode
     * @param root root vnode ID, most likely 0
     * @return all vnodes recieve this int
     */
    public double broadcast(double value, int root);
    
    /**
     * Takes a list of Serilizable objects from all vnodes in the communicator and
     * creates an ordered contiguous list of those objects in each vnode.
     * Requires the return object to be cast into a more useful object.
     *
     * @param myElements     in on entry; the list of values to be sent to all vnodes
     * @return allElements   out on exit; the list of values from all vnodes
     */
    public Serializable[] gatherAll(Serializable [] myElements);

    /**
     * Takes a list of values from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myElements     in on entry; the list of values to be sent to all vnodes
     * @return allElements   out on exit; the list of values from all vnodes
     */
    public double[] gatherAll(double [] myElements);
    
    /**
     * Takes a list of values from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myElements    in on entry; the list of values to be sent to all vnodes
     * @return              out on exit; the list of values from all vnodes
     */
    public int[] gatherAll(int [] myElements);
    
    /**
     * Takes a value from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myInt     in on entry; the list of values to be sent to all vnodes
     * @return          out on exit; the list of values from all vnodes
     */
    public int[] gatherAll(int myInt);
    
    /**
     * Takes a value from all vnodes in the communicator and
     * creates an ordered contiguous list of those values in each vnode.
     *
     * @param myDouble     in on entry; the list of values to be sent to all vnodes
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
    
    /**
     * Takes a list of values form all vnodes in the communicator, computes
     * the min, and returns the min to all vnodes.
     *
     * @param partialMins   in on entry; the list of values, usually computed
     *                      locally, to be summed across all vnodes
     * @return              out on exit; the list of values summed across all vnodes
     */
    public double[] minAll(double [] partialMins);
    
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
     * the scan sum, and returns it to all vnodes such that vnodeor i contains
     * the sum of values from vnodes up to and including vnodeor i.
     *
     * @param myElements    in on entry; the values to be summed across all vnodes
     * @return              out on exit; the list of values summed
     *                      across vnodes 0 through i
     */
    public double[] scanSums(double [] myElements);
    
    /**
     * Takes a list of values from all vnodes in the communicator, computes
     * the scan sum, and returns it to all vnodes such that vnodeor i contains
     * the sum of values from vnodes up to and including vnodeor i.
     *
     * @param myElements    in on entry; the values to be summed across all vnodes
     * @return              out on exit; the list of values summed
     *                      across vnodes 0 through i
     */
    public int[] scanSums(int [] myElements);

    /**
     * Accessor for the number of threads in this vnode.
     *
     * @return number of threads in this vnode
     */
    public int getNumMyThreads();
    
    /**
     * Accessors for the number of threads in all proccesses.
     *
     * @return number of threads in all proccesses
     */
    public int getNumThreads();
    
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
    public int getVnodeID();
    
    /**
     * Ask about this one.
     */
    public int getThreadID();
    
    /**
     * Sets the thread ID for the calling vnode.
     * Can be used to facilitate threaded programming across an MPI
     * application by allowing multiple MPI vnodes to be considered
     * threads of a virtual shared memory vnode. Threads and vnodes
     * should be used together.
     *
     * @param newThreadID new thread ID
     */
    public int setThreadID(int newThreadID);
    
    /**
     * Sets the vnode ID for the calling vnode.
     * Can be used to facilitate threaded programming acroass an MPI
     * application by associating several MPI vnodes with a single
     * vnode. By default, each MPI vnode is associated with a single
     * vnode wit the same ID.
     * 
     * @param newVnodeID new vnode ID
     */
    public int setMyVnodeID(int newVnodeID);
    
    /**
     * Accessor to see if the communicator is on a serial or parallel machine.
     *
     * @return <code>true</code> in serial mode, otherwise <code>false</code>
     */
    public boolean getIsSerial();    
}

