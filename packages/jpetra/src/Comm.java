package Jpetra;
/*
 * Comm.java
 *
 * Created on May 29, 2001, 1:25 PM
 */


/**
 * Comm is the interface for the all the Jpetra communication
 * classes.
 *
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @version 
 */
public interface Comm {
    
    /**
     * Causes each processor in the communicator to wait until all
     * processors have arrived.
     * No-op for a serial communicator. 
     */
    public void barrier();
    
    /**
     * Causes each process on a given node in the communicator
     * to wait until all processes on that node have arrived.
     * No-op for a serial communicator.
     */
    public void nodeBarrier();
    
    /**
     * Takes a list of input values from the root processor
     * and sends them to all other processors.
     * No-op for a serial communicator.
     *
     * @param numElements   in on entry; the number of elements in the elements array
     * @param elements      in on entry; the array of elements to process
     * @param root          in on entry; the processor from which all processors 
     *                      will receive a copy of the elements
     */
    public int broadcast(int numElements, double [] elements, int root);
    
    /**
     * Takes a list of input values from the root processor
     * and sends them to all other processors.
     * No-op for a serial communicator.
     *
     * @param numElements   in on entry; the number of elements in the elements array
     * @param elements      in on entry; the array of elements to process
     * @param root          in on entry; the processor from which all 
     *                      processors will receive a copy of the elements
     */
    public int broadcast(int numElements, int [] elements, int root);
    
    /**
     * Takes a list of values from all porcessors in the communicator and
     * creates an ordered contiguous list of those values on each processor.
     *
     * @param numElements   in on entry; the number of elements in myElements
     * @param myElements    in on entry; the list of values to be sent to all processors
     * @param allElements   out on exit; the list of values from all processors
     */
    public int gatherAll(int numElements, double [] myElements, double [] allElements);
    
    /**
     * Takes a list of values from all porcessors in the communicator and
     * creates an ordered contiguous list of those values on each processor.
     *
     * @param numElements   in on entry; the number of elements in myElements
     * @param myElements    in on entry; the list of values to be sent to all processors
     * @param allElements   out on exit; the list of values from all processors
     */
    public int gatherAll(int numElements, int [] myElements, int [] allElements);
    
    /**
     * Takes a list of values from all processors in the communicator, computes
     * the sum, and returns the sum to all processors.
     *
     * @param numElements   in on entry; the length of partialSums
     * @param partialSums   in on entry; the list of values, usually computed
     *                      locally, to be summed across all processors
     * @param globalSums    out on exit; the list of values summed across all processors
     */
    public int sumAll(int numElements, double [] partialSums, double [] globalSums);
    
    /**
     * Takes a list of values from all processors in the communicator, computes
     * the sum, and returns the sum to all processors.
     *
     * @param numElements   in on entry; the length of partialSums
     * @param partialSums   in on entry; the list of values, usually computed
     *                      locally, to be summed across all processors
     * @param globalSums    out on exit; the list of values summed
     *                      across all processors
     */
    public int sumAll(int numElements, int [] partialSums, int [] globalSums);
    
    /**
     * Takes a list of values form all processors in the communicator, computes
     * the max, and returns the max to all processors.
     *
     * @param numElements   in on entry; the length of partialMaxs
     * @param partialMaxs   in on entry; the list of values, usually computed
     *                      locally, to be summed across all processors
     * @param globalMaxs    out on exit; the list of values summed across all processors
     */
    public int maxAll(int numElements, double [] partialMaxs, double [] globalMaxs);
    
    /**
     * Takes a list of values form all processors in the communicator, computes
     * the max, and returns the max to all processors.
     *
     * @param numElements   in on entry; the length of partialMaxs
     * @param partialMaxs   in on entry; the list of values, usually computed
     *                      locally, to be summed across all processors
     * @param globalMaxs    out on exit; the list of values summed across all processors
     */
    public int maxAll(int numElements, int [] partialMaxs, int [] globalMaxs);
    
    /**
     * Takes a list of values form all processors in the communicator, computes
     * the min, and returns the min to all processors.
     *
     * @param numElements   in on entry; the length of partialMins
     * @param partialMins   in on entry; the list of values, usually computed
     *                      locally, to be summed across all processors
     * @param globalMaxs    out on exit; the list of values summed across all processors
     */
    public int minAll(int numElements, double [] partialMins, double [] globalMins);
    
    /**
     * Takes a list of values form all processors in the communicator, computes
     * the min, and returns the min to all processors.
     *
     * @param numElements   in on entry; the length of partialMins
     * @param partialMins   in on entry; the list of values, usually computed
     *                      locally, to be summed across all processors
     * @param globalMaxs    out on exit; the list of values summed across all processors
     */
    public int minAll(int numElements, int [] partialMins, int [] globalMins);
    
    /**
     * Takes a list of values from all processors in the communicator, computes
     * the scan sum, and returns it to all processors such that processor i contains
     * the sum of values from processors up to and including processor i.
     *
     * @param numElements   in on entry; the length of myElememtns
     * @param myElements    in on entry; the values to be summed across all processors
     * @param sums          out on exit; the list of values summed
     *                      across processors 0 through i
     */
    public int scanSums(int numElements, double [] myElements, double [] sums);
    
    /**
     * Takes a list of values from all processors in the communicator, computes
     * the scan sum, and returns it to all processors such that processor i contains
     * the sum of values from processors up to and including processor i.
     *
     * @param numElements   in on entry; the length of myElememtns
     * @param myElements    in on entry; the values to be summed across all processors
     * @param sums          out on exit; the list of values summed
     *                      across processors 0 through i
     */
    public int scanSums(int numElements, int [] myElements, int [] sums);
    
    /**
     * Ask about this one.
     */
    public int getMyNode();

    /**
     * Accessor for the number of nodes in the commmunicator.
     */
    public int getNumNodes();
    
    /**
     * Accessor for the number of processors in the communicator.
     */
    public int getNumProc();
    
    /**
     * Accessor for the rank of the calling process.
     *
     * @return the rank of the calling process in MPI; 0 in serial mode
     */
    public int getPID();
    
    /**
     * Ask about this one.
     */
    public int getThreadID();
    
    /**
     * Sets the thread ID for the calling process.
     * Can be used to facilitate threaded programming across an MPI
     * application by allowing multiple MPI processes to be considered
     * threads of a virtual shared memory process. Threads and nodes
     * should be used together.
     */
    public int setThreadID(int newThreadID);
    
    /**
     * Sets the node ID for the calling process.
     * Can be used to facilitate threaded programming acroass an MPI
     * application by associating several MPI processes with a single
     * node. By default, each MPI process is associated with a single
     * node wit the same ID.
     */
    public int setMyNodeID(int newNodeID);
    
    /**
     * Accessor to see if the communicator is on a serial or parallel machine.
     *
     * @return <code>true</code> in serial mode, otherwise <code>false</code>
     */
    public boolean getIsSerial();    
}

