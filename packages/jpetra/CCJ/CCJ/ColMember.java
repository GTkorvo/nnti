/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

// import java.rmi.server.*;
// import java.rmi.registry.*;
import java.rmi.Naming;
import java.rmi.AlreadyBoundException;
// import java.net.*;
import java.net.InetAddress;
import java.io.Serializable;

public abstract class ColMember implements ColMemberInterface, Runnable {

    private static final boolean DEBUG = false;
    private static final boolean ASSERT = false;
    private static final boolean DEBUG2 = false; // Debug level 2
    private static final boolean DO_TIMINGS = false; // Also update ThreadPool.DO_TIMINGS
    private boolean USE_THREADPOOL = false;
    // DEFAULT_POOLSIZE is number of threads in the threadpool
    private static final int DEFAULT_POOLSIZE = 2;
    // DEFAULT_NUMBEROFHASHBUCKETS must be a power of 2
    private static final int DEFAULT_NUMBEROFHASHBUCKETS = 16;
    private static final int DEFAULT_NUMBEROFPREALLOCATEDMESSAGES = 10;

    // DEBUG timing variable
    private long receiveTime = 0;

    private Thread myThread;
    protected MyList mqueue;

    private String hostName;
    private String lookupName;
    private ThreadPool threadPool;

    private ColGroup cachedGroup;
    private int cachedRank;
    private int cachedSize;
    // If not using threadpool... a message is needed to send data
    private static int numberOfMembersDone = 0;
    private static int numberOfLocalMembers = 0;
    private String myId;
    private String myPort;
    
    boolean record = false;
    int[] records;
    int numRecords;
    
    public ColMember() throws CCJException {
	this(false, DEFAULT_POOLSIZE, DEFAULT_NUMBEROFHASHBUCKETS,
	     DEFAULT_NUMBEROFPREALLOCATEDMESSAGES);
    }


    int cachedRank() {
	return cachedRank;
    }


    int cachedSize() {
	return cachedSize;
    }


    public ColMember(boolean USE_THREADPOOL,
		     int threadPoolSize,
		     int numberOfHashBuckets,
		     int numberOfPreallocatedMessages)
	    throws CCJException {

	boolean ok;

	this.USE_THREADPOOL = USE_THREADPOOL;
	cachedGroup = null;
	cachedRank = -1;
	cachedSize = 0;

	if (DEBUG2) {
	    System.out.println("New ColMember is being created");
	}

	myThread = null;

	threadPool = new ThreadPool(threadPoolSize);

	try {
	    mqueue = new MyList(this, numberOfHashBuckets, numberOfPreallocatedMessages);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException("Fail to create my message queue " + e);
	}

	/*try {
	    hostName = InetAddress.getLocalHost().getHostName();
	} catch (java.net.UnknownHostException e) {
	    throw new CCJException("Couldn't get local hostname : " + e);
	}*/

	myId = ColGroupMaster.getId();
	
    
	lookupName = "//" + ColGroupMaster.hostPort() + "/"+myId;
	ok = false;

	while (!ok) {
	    try {
		ok = true;
		Naming.bind("//" + ColGroupMaster.hostPort() + "/" + myId, mqueue);
	    } catch (java.net.MalformedURLException e) {
		throw new CCJException("Not a well-formed host name: " + myId + " : " + e);
	    } catch (AlreadyBoundException e) {
		ok = false;
		myId = ColGroupMaster.getId();
		lookupName = "//" + ColGroupMaster.hostPort() + "/" + myId;
System.err.println("bind throws " + e);
	    } catch (java.rmi.RemoteException e) {
		ok = false;
		myId = ColGroupMaster.getId();
		lookupName = "//" + ColGroupMaster.hostPort() + "/" + myId;
System.err.println("bind throws " + e);
	    }
	}

	if (ColGroupMaster.groupCentral == null) {
	    throw new CCJException("Internal CCJ error: ColGroupMaster not initialized yet");
	}

	synchronized (this) {
	    numberOfLocalMembers++;
	}

	if (DEBUG) {
	    System.out.println("Colmember: " + lookupName + " created. "+
			       "DO_TIMINGS = " + DO_TIMINGS +
			       " USE_THREADPOOL = " + USE_THREADPOOL +
			       " numberOfThreads = " + threadPool.numberOfThreads);
	    System.out.flush();
	}
    }


    public void begin() throws CCJException {
	myThread = new Thread(this);
	myThread.start();
    }


    void done() throws CCJException {
	threadPool.waitForThreads();

//          StringBuffer s = new StringBuffer();
//          s.append("****************** MEMBER " + cachedRank);
//          s.append(" ************************************\n");
//          s.append("USE_THREADPOOL              : " + USE_THREADPOOL + "\t\t");
//          s.append("threadPoolSize              : " +
//                   threadPool.numberOfThreads + "\n");
//          s.append("numberOfHashBuckets         : " + mqueue.BUCKETS + "\t\t");
//          s.append("numberOfPreallocatedMessages: " +
//                   mqueue.STARTING_POOLSIZE + "\n");
//          s.append("max queue-size              : " +
//                   (mqueue.poolSize() - 1) + "\t\t");
//          s.append("NumberOfNext                : " +
//                   mqueue.numberOfNext + "\n");
//          s.append("pool.hits                   : " + threadPool.hits + "\t\t\t");
//          s.append("pool.misses                 : " + threadPool.misses +
//                   "\n");

//          if (DO_TIMINGS) {
//              s.append("SendTime                    : " +
//                       threadPool.sendTime + "\t\t");
//              s.append("ReceiveTime                 : " + receiveTime + "\n");
//              s.append("ThreadPool.waitTime         : " +
//                       threadPool.waitTime + "\n");
//          }

//          s.append("******************************************************"+
//                   "**********\n");
//          System.out.print(s);

	if (mqueue.size() > 0) {
	    System.err.println("WARNING: ColObject: "+ lookupName + " has "+
			       mqueue.size() + " messages left after done!");
	}

	threadPool.done();

	// Wait for all local colmembers to stop.
	int i;

	synchronized (this) {
	    i = ++numberOfMembersDone;

	    while (numberOfMembersDone < numberOfLocalMembers) {
		try {
		    wait();
		} catch (InterruptedException e) {
		    // Do nothing
		}

	    }

	    notifyAll();
	}

	if (i == numberOfLocalMembers) {
	    if (DEBUG) {
		System.out.println("ColGroupCentral: Waiting for members"+
				   " to end.");
	    }

	    try {
		// Wait for all threads to end.
		Thread.currentThread().sleep(1000);
	    } catch (InterruptedException e) {
		// Do nothing
	    }

	    Thread.yield();

	    // System.exit(0);
	}

	try {
	    Naming.unbind(myId);
	} catch (java.rmi.NotBoundException e) {
	    throw new CCJException("unbind(" + myId + ") fails; " + e);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException("unbind(" + myId + ") fails; " + e);
	} catch (java.net.MalformedURLException e) {
	    throw new CCJException("unbind(" + myId + ") fails; " + e);
	}
    }


    private Serializable receive(int tag, int source) throws CCJException {
	Serializable data;
	long begin = 0, end = 0;

	if (DO_TIMINGS) {
	    begin = System.currentTimeMillis();
	}

	data = mqueue.getFirst(tag, source);

	if (DO_TIMINGS) {
	    end = System.currentTimeMillis();
	    receiveTime += end - begin;
	}

	if (DEBUG) {
	    System.out.println(cachedRank + ": Got message from "+
			       source + ", tag= "+ tag + " data = "+ data);
	}
    
    if (this.record) {
        this.records[numRecords++] = source;
    }
    
	return data;
    }


    private void send(ColGroup group, int destination, Serializable objToSend, int sender, int tag) throws java.rmi.RemoteException {

	long beginTime = 0;
	MyListInterface dst = (MyListInterface)group.memberLists.elementAt(destination);

	if (DO_TIMINGS) {
	    beginTime = System.currentTimeMillis();
	}

//        if (USE_THREADPOOL) {
//            threadPool.send(tag, sender, objToSend, dst, beginTime);
//            Thread.yield();
//        } else {

	    if (DEBUG) {
		System.out.println("ColMember.send: "+ sender + ": Really sending NO_THREAD.");
	    }

	    dst.add(tag, sender, objToSend);

	    if (DO_TIMINGS) {
		long endTime = System.currentTimeMillis();
		threadPool.sendTime += endTime - beginTime;
	    }
//        }
    }


    private void confirmedSend(MyListInterface dst,
			       Serializable objToSend,
			       int sender,
			       int tag)
	    throws java.rmi.RemoteException, CCJException {
	// Send a blocking RMI to other member. Used for barrier.
	long beginTime = 0;

	if (DO_TIMINGS) {
	    beginTime = System.currentTimeMillis();
	}

	if (DEBUG) {
	    System.out.println(sender + ": Really sending CONFIRMED. tag =" + tag);
	}

	dst.confirmedAdd(tag, sender, objToSend);

	if (DO_TIMINGS) {
	    long endTime = System.currentTimeMillis();
	    threadPool.sendTime += endTime - beginTime;
	}

	if (DEBUG) {
	    System.out.println(sender + ": Confirmed send DONE.tag ="+ tag);
	}
    }


    public void send_sync(ColGroup group, Serializable objectToDistribute, int destination) throws CCJException, NoSuchMemberException {

	int tag = group.getNewUnicastSendTag(destination);

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("sender rank " + rank);
	}
	if (destination < 0 || destination >= size) {
	    throw new NoSuchMemberException("destination rank " + destination);
	}

	try {
	    confirmedSend(((MyListInterface) group.memberLists.elementAt(destination)), objectToDistribute, rank, tag);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public void send_async(ColGroup group,
			   Serializable objectToDistribute,
			   int destination)
	    throws CCJException, NoSuchMemberException {

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("sender rank " + rank);
	}
	if (destination < 0 || destination >= size) {
	    throw new NoSuchMemberException("destination rank " + destination);
	}

	int tag = group.getNewUnicastSendTag(destination);

	MyListInterface dest = (MyListInterface) group.memberLists.elementAt(destination);
	try {
	    send(group, destination, objectToDistribute, rank, tag);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public Serializable receive(ColGroup group, int source)
	    throws CCJException, NoSuchMemberException {

	int tag = group.getNextUnicastReceiveTag(source);
	return receive(tag, source);
    }


    public Serializable send_receive(ColGroup send_group,
				     Serializable objectToDistribute,
				     int destination,
				     ColGroup receive_group,
				     int source)
	    throws CCJException, NoSuchMemberException {

	send_async(send_group, objectToDistribute, destination);
	return receive(receive_group, source);
    }


    public Serializable rendezVous(ColGroup group, Serializable objectToDistribute, int destination) throws CCJException, NoSuchMemberException {

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank < 0 || rank >= cachedSize) {
	    throw new NoSuchMemberException("sender rank " + rank);
	}
	if (destination < 0 || destination >= size) {
	    throw new NoSuchMemberException("destination rank " + destination);
	}

	if (rank < destination) {
	    try {
		int tag = group.getNewUnicastSendTag(destination);
		MyListInterface dest = (MyListInterface) group.memberLists.elementAt(destination);
		return dest.rendezVousAdd(tag, rank, objectToDistribute);
	    } catch (java.rmi.RemoteException e) {
		throw new CCJException(e.toString());
	    }
	} else {
	    int tag = group.getNextUnicastReceiveTag(destination);
	    return mqueue.rendezVousGet(tag, destination, objectToDistribute);
	}
    }


    public Serializable broadcast(ColGroup group, Serializable objectToDistribute, int root) throws CCJException, NoSuchMemberException {
	return broadcast(group, objectToDistribute, root, ColGroup.BROADCASTTAG);
    }


    private Serializable broadcast(ColGroup group, Serializable objectToDistribute, int root, int tag) throws CCJException, NoSuchMemberException {
    
    if (DEBUG) System.out.println("Doing broadcast");
    
	int rank_rel;
	int mask, sum;
	Serializable result = null;

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("local rank " + rank);
	}
	if (root < 0 || root >= size) {
	    throw new NoSuchMemberException("root rank " + root);
	}

	if (size == 1) {
	    return objectToDistribute;
	}

	int local_tag = group.getNewGroupTag(tag);

	if (DEBUG) {
	    System.out.println(rank + ": broadcast started, root = " +
			       root + " groupsize="+ size + " counter << = " +
			       local_tag);
	}

	if (rank == root) {
	    result = objectToDistribute;
	}

	rank_rel = rel_rank(rank, root, size);

	for (mask = 1; mask < size; mask *= 2);

	mask /= 2;
	sum = 0;

	try {
	    while (mask > 0) {
		if (sum + (mask & rank_rel) == rank_rel) {
		    /* do something in this step */

		    if ((mask & rank_rel) != 0) {
			if (DEBUG2) {
			    System.out.println(rank + ": ColMember.broadcast: receive from "+ abs_rank(rank_rel - mask, root, size));
			}

			result = receive(local_tag, abs_rank(rank_rel - mask, root, size));
		    } else if ((rank_rel + mask) < size) {
			if (DEBUG2) {
			    System.out.println(rank + ": ColMember.broadcast: send to " + abs_rank(rank_rel + mask, root, size));
			}
			send(group, abs_rank(rank_rel + mask, root, size), result, rank, local_tag);
		    }
		}

		sum += (mask & rank_rel);
		mask /= 2;
	    }
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}

	if (DEBUG) {
	    System.out.println(rank + ": broadcast done.");
	}
	return result;
    }


    public Serializable reduce(ColGroup group,
			       Serializable dataObject,
			       Reducible reductionObject,
			       int root)
	    throws CCJException, NoSuchMemberException {
	return reduce(group, dataObject, reductionObject, root, ColGroup.REDUCETAG);
    }


    private Serializable reduce(ColGroup group,
				Serializable dataObject,
				Reducible reductionObject,
				int root,
				int tag)
	    throws CCJException, NoSuchMemberException {

    if (DEBUG) System.out.println("Doing reduce");
	        
	int mask, relrank, peer;
	Serializable result = dataObject;
	Serializable data;

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("local rank " + rank);
	}
	if (root < 0 || root >= size) {
	    throw new NoSuchMemberException("root rank " + root);
	}

	int local_tag = group.getNewGroupTag(tag);

	if (DEBUG) {
	    System.out.println(rank + ": Reduce started. Root = " + root + " groupsize = " + size + " counter<< = " + local_tag);
	}

	mask = 1;

	relrank = (rank - root + size) % size;

	try {
	    while (mask < size) {
		if ((mask & relrank) == 0) { /* receive and reduce */
		    peer = (relrank | mask);

		    if (peer < size) {
			peer = (peer + root) % size;

			if (DEBUG2) {
			    System.out.println(rank + ": Reduce receive from "+ peer);
			    System.out.flush();
			}

			data = receive(local_tag, peer);

			if (DEBUG) {
			    System.out.println(rank + ": reduce call "+ reductionObject + " arg1 = "+ result + " arg2 = "+ data);
			    System.out.flush();
			}

			result = reductionObject.reduce(result, data);
			// DebugInterface(result);
		    }
		} else { /* send and terminate */
		    peer = ((relrank & (~mask)) + root) % size;

		    if (DEBUG2) {
			System.out.println(rank + ": Reduce send to "+ peer);
			System.out.flush();
		    }

		    // DebugInterface(result);
		    send(group, peer, result, rank, local_tag);

		    break;
		}

		mask <<= 1;
	    }
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}

	if (DEBUG) {
	    System.out.println(rank + ": reduce done.");
	}

	if (rank == root) {
	    return result;
	} else {
	    return null;
	}
    }


//    native void DebugInterface(Object s);


    public Serializable allReduce(ColGroup group,
				  Serializable dataObject,
				  Reducible reductionObject)
	    throws CCJException, NoSuchMemberException {

	int size = group.size();

    if (DEBUG) System.out.println("Doing allReduce");
    
	// If the number of processors, we prefer to do a butterfly allReduce.
	// That is most efficient.
	/*MODIFIED
	The allReduce_Butterfly method seems to be flawed and needs to be further tested and fixed.
	*/
	for (int i = 0; i < 31; i++) {
	    if (size == (1 << i)) {
	//    System.out.println("Doing butterfly reduce");
		return allReduce_butterfly(group, dataObject, reductionObject);
	    }
	}
    
    //System.out.println("Doing bcast reduce");
	return allReduce_reduce_bcast(group, dataObject, reductionObject);
    }


    private Serializable allReduce_reduce_bcast(
				ColGroup group,
				Serializable dataObject,
				Reducible reductionObject)
	    throws CCJException, NoSuchMemberException {
	Serializable result;
    if (DEBUG) System.out.println("Doing allReduce_reduce_bcast");
	
	result = reduce(group, dataObject, reductionObject, 0, ColGroup.ALLREDUCETAG1);
    
//        DebugInterface(result);
	return broadcast(group, result, 0, ColGroup.ALLREDUCETAG2);

    }


    private Serializable allReduce_butterfly(
				ColGroup group,
				Serializable dataObject,
				Reducible reductionObject)
	    throws CCJException, NoSuchMemberException {

    if (DEBUG) System.out.println("Doing allReduce_butterfly");
    
	int local_tag = group.getNewGroupTag(ColGroup.ALLREDUCETAG1);

	Serializable local_result, remote_result;

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException();
	}

	if (size == 1) {
	    return dataObject;
	}

	boolean even = ((rank & 1) == 0);
	int next_rank, prev_rank;
	int send_size = 1;
	int offset = 0;

	local_result = dataObject;

	if (even) {

	    try {
		while (send_size < size) {

		    offset += send_size;
		    next_rank = (rank + offset) % size;

		    // exchange data with the machine which is 2^phase further away
		    MyListInterface list = (MyListInterface) group.memberLists.elementAt(next_rank);
		    remote_result = list.rendezVousAdd(local_tag, rank, local_result);
		    // System.out.println(rank + ": Reduce objects: local = " + local_result + " remote Add from " + next_rank + " = " + remote_result);
		    local_result = reductionObject.reduce(local_result, remote_result);
		    send_size *= 2;
		}
	    } catch (java.rmi.RemoteException e) {
		throw new CCJException(e.toString());
	    }

	} else {

	    while (send_size < size) {

		offset += send_size;
		prev_rank = (rank + size - offset) % size;

		// exchange data with the machine which is 2^phase 'in front' away
		remote_result = mqueue.rendezVousGet(local_tag, prev_rank, local_result);
// System.out.println(rank + ": Reduce objects: local = " + local_result + " remote Get from " + prev_rank + " = " + remote_result);
		local_result = reductionObject.reduce(local_result, remote_result);
		send_size *= 2;
	    }
	}

	return local_result;
    }


    public Partitionable flatGather(ColGroup group,
				    Partitionable rootObject,
				    Serializable dataObject,
				    int root)
	    throws CCJException, NoSuchMemberException {

	// For some reason he starts with a FLATGATHERTAG but then switches to GATHERTAG when doing the sending and receiving ???

    if (DEBUG) System.out.println("Doing flatGather");
	// Just do a flat gather. So every member sends its data to root
	int local_tag = group.getNewGroupTag(ColGroup.FLATGATHERTAG);

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("local rank " + rank);
	}
	if (root < 0 || root >= size) {
	    throw new NoSuchMemberException("root rank " + root);
	}

	if (rank != root) {
	    try {
		send(group, root, dataObject, rank, local_tag);
	    } catch (java.rmi.RemoteException e) {
		throw new CCJException(e.toString());
	    }
	    return null;
	}

	for (int i = 0 ; i < size ; i++) {
	    Serializable data;

	    if (i == root) {
		rootObject.setElementAt(i, size, dataObject);
	    } else {
		data = receive(local_tag, i);
		rootObject.setElementAt(i, size, data);
	    }
	}

	return rootObject;
    }


//    native void DebugInterface(Object s);
//    native void DebugObj(Object o);

    public Partitionable gather(ColGroup group,
				Partitionable rootObject,
				Serializable dataObject,
				int root)
	    throws CCJException, NoSuchMemberException {
	// Don't do a flat gather. But an altered reduce.

	//System.err.println("1");

    if (DEBUG) System.out.println("Doing gather");
	int mask, relrank, peer;
	Serializable result = dataObject;

	Serializable [] resultArray = null;
	Serializable data;

	int local_tag = group.getNewGroupTag(ColGroup.GATHERTAG);

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("local rank " + rank);
	}
	if (root < 0 || root >= size) {
	    throw new NoSuchMemberException("root rank " + root);
	}

	if (DEBUG) {
	    System.out.println(rank + ": gather started. Root = "+ root + " groupsize = " + size + " counter<< = " + local_tag);
	}

	if (size == 1) {
	    rootObject.setElementAt(0, 1, dataObject);
	    return rootObject;
	}

	mask = 1;

	relrank = (rank - root + size) % size;

	try {
	    while (mask < size) {
		if ((mask & relrank) == 0) {

		    //System.err.println("4A");

		    peer = (relrank | mask);

		    if (peer < size) {

			peer = (peer + root) % size;

			if (DEBUG2) {
			    System.out.println("gather: "+rank + " waiting for data from: " + peer);
			}

			data = receive(local_tag, peer);

			if (DEBUG) {
			    System.out.println(rank + ": ColMember.gather: got "+ "resultarray from " + peer);
			}

			if (resultArray == null) {
			    resultArray = (Serializable [])data;
			} else {
			    // Merge the two arrays
			    Serializable[] tmp = (Serializable[]) data;

			    for (int i = 0 ; i < size ; i++) {
				if (tmp[i] != null) {
				    resultArray[i] = tmp[i];
				}
			    }
			}
		    }

		} else { /* send and terminate */
		    //System.err.println("4B1");

		    peer = ((relrank & (~mask)) + root) % size;

		    if (DEBUG2) {
			System.out.println("Gather: "+rank + " sending data to "+ peer);
		    }

		    if (resultArray == null) {
			// This member is a leaf node.
			resultArray = new Serializable[size];
		    }

		    resultArray[rank] = dataObject; // Add data to array
		    send(group, peer, resultArray, rank, local_tag);
		    resultArray = null;

		    break;
		}

		mask <<= 1;
	    }
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}

	//System.err.println("5");

	if (DEBUG) {
	    System.out.println(rank + ": gather done.");
	}

	//System.err.println("6");

	if (rank != root) {
	    //System.err.println("7B");
	    return null;
	}

	//System.err.println("7A1");

	// Add object of root to array.
	resultArray[root] = dataObject;

	//System.err.println("7A2");

	// Let the root node create the new object.
	for (int i = 0 ; i < size ; i++) {
	    //System.err.println("7A21");

	    //DebugObj( resultArray[i] );
	    rootObject.setElementAt(i, size, resultArray[i]);
	}

	//System.err.println("7A3");
	return rootObject;
    }


    public Partitionable allGather(ColGroup group,
				   Partitionable rootObject,
				   Serializable dataObject)
	    throws CCJException, NoSuchMemberException {
	int size = group.size();

    if (DEBUG) System.out.println("Doing allGather");
	// If the number of processors, we prefer to do a butterfly allGather.
	// That is most efficient.
	for (int i = 0; i < 31; i++) {
	    if (size == (1 << i)) {
		return allGather_butterfly(group, rootObject, dataObject);
	    }
	}

	return allGather_double_ring(group, rootObject, dataObject);
    }


    private Partitionable allGather_butterfly(ColGroup group,
		Partitionable rootObject,
		Serializable dataObject)
	    throws CCJException, NoSuchMemberException {

    if (DEBUG) System.out.println("Doing allGather_butterfly");
	// A butterfly allGather. Works only if the number of members is
	// a power of two.
	//
	// Do a RMI-aware gather. Instead of assuming that we have an asynchronous send (which we don't have),
	// we make use of the fact that RMI always sends a reply. As a result, we can do an exchange allgather,
	// where two machines exchange the date they have at the moment (in one RMI).

	int local_tag = group.getNewGroupTag(ColGroup.ALLGATHERTAG);

	Serializable[] local_results, remote_results;

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException();
	}

	if (size == 1) {
	    rootObject.setElementAt(0, 1, dataObject);
	    return rootObject;
	}

	boolean even = ((rank & 1) == 0);

	local_results = new Serializable[size];
	local_results[rank] = dataObject;

	int next_rank, prev_rank;
	int send_size = 1;
	int offset = 0;

	if (even) {
	    try {
		while (send_size < size) {

		    offset += send_size;
		    next_rank = (rank + offset) % size;

		    // exchange data with the machine which is 2^phase further away
		    MyListInterface list = (MyListInterface) group.memberLists.elementAt(next_rank);

		    if (ASSERT) {
			for (int i = 0; i < send_size; i++) {
			    int x = (rank + i) % size;
			    if (local_results[x] == null) {
				throw new CCJException(rank + ": Internal error: In phase -pow " + send_size + " I haven't got data[" + x + "] to send to " + next_rank);
			    }
			}
		    }

		    remote_results = (Serializable[])list.rendezVousAdd(local_tag, rank, local_results);

		    for (int j=0;j<size;j++) {
			if (remote_results[j] != null) {
			    local_results[j] = remote_results[j];
			}
		    }

		    send_size *= 2;

		    if (ASSERT) {
			for (int i = 0; i < send_size; i++) {
			    int x = (rank + i) % size;
			    if (local_results[x] == null) {
				throw new CCJException(rank + ": Internal error: In phase -pow " + (send_size / 2) + " I haven't received data[" + x + "] from " + next_rank);
			    }
			}
		    }
		}
	    } catch (java.rmi.RemoteException e) {
		throw new CCJException(e.toString());
	    }

	} else {
	    while (send_size < size) {

		offset += send_size;
		prev_rank = (rank + size - offset) % size;

		if (ASSERT) {
		    for (int i = 0; i < send_size; i++) {
			int x = (rank + size - i) % size;
			if (local_results[x] == null) {
			    throw new CCJException(rank + ": Internal error: In phase -pow " + send_size + " I haven't got data[" + x + "] to send to " + prev_rank);
			}
		    }
		}

		// exchange data with the machine which is 2^phase 'in front' away
		remote_results = (Serializable[])mqueue.rendezVousGet(local_tag, prev_rank, local_results);

		for (int j=0;j<size;j++) {
		    if (remote_results[j] != null) {
			local_results[j] = remote_results[j];
		    }
		}

		send_size *= 2;

		if (ASSERT) {
		    for (int i = 0; i < send_size; i++) {
			int x = (rank + size - i) % size;
			if (local_results[x] == null) {
			    throw new CCJException(rank + ": Internal error: In phase -pow " + (send_size / 2) + " I haven't received data[" + x + "] from " + prev_rank);
			}
		    }
		}
	    }
	}

	for (int i=0;i<size;i++) {
	    rootObject.setElementAt(i, size, local_results[i]);
	}

	int log_size = 0;
	for (log_size = 0; (1 << log_size) < size; log_size++);

	return rootObject;
    }


    private Partitionable allGather_double_ring(
		ColGroup group,
		Partitionable rootObject,
		Serializable dataObject)
	    throws CCJException, NoSuchMemberException {

    if (DEBUG) System.out.println("Doing allGather_double_ring");
	// Do a double ring allGather. Less efficient than a butterfly
	// allGather, but works for any number of hosts.
	// Algorithm:
	//   loop invariant:
	//      in iteration i, I possess data[me - i, ..., me + i].
	//      I receive data[me - i - 1] from my lower neighbour,
	//      data[me + i + 1] from my upper neighbour.
	//      I post data[me + i] for my lower neighbour,
	//      data[me - i] for my upper neighbour.
	//
	// Do a RMI-aware gather. Instead of assuming that we have an asynchronous send (which we don't have),
	// we make use of the fact that RMI always sends a reply. As a result, we can do an exchange allgather,
	// where two machines exchange the date they have at the moment (in one RMI).

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException();
	}

	rootObject.setElementAt(rank, size, dataObject);

	if (size == 1) {
	    return rootObject;
	}

	int local_tag = group.getNewGroupTag(ColGroup.ALLGATHERTAG);

	int next_rank = (rank + 1) % size;
	int prev_rank = (rank + size - 1) % size;

	int rounds     = (size - 1) / 2;
	boolean add_half_round = (((size - 1) % 2) != 0);

	MyListInterface next_list = (MyListInterface) group.memberLists.elementAt(next_rank);

	Serializable sending_back    = dataObject;
	Serializable sending_forward = dataObject;

	try {
	    for (int i = 0; i < rounds; i++) {
		// Invariant:
		//     sending_back = data[me + i]
		//     sending_forward = data[me - i]
		// then from the exchange:
		//     sending_back := data[me + i + 1]
		//     sending_forward := data[me - i - 1]

		// System.out.println("Doing futureSet");

		mqueue.futureSet(sending_back);

		// System.out.println("Doing futureAdd");

		sending_back = next_list.futureAdd(sending_forward);
		rootObject.setElementAt((rank + (i + 1)) % size, size,
					sending_back);

		// System.out.println("Doing futureGet");
		sending_forward = mqueue.futureGet();
		rootObject.setElementAt((rank + size - (i + 1)) % size, size,
					sending_forward);

		// System.out.println("Doing futureDone");
	    }

	    if (add_half_round) {
		if ((rank & 1) == 0) {
		    // we must send
		    sending_back = next_list.futureAdd(sending_forward);
		    rootObject.setElementAt((rank + size - rounds - 1) % size,
					    size, sending_back);
		} else {
		    // we must receive
		    mqueue.futureSet(sending_back);
		    sending_forward = mqueue.futureGet();
		    rootObject.setElementAt((rank + rounds + 1) % size, size,
					    sending_forward);
		}
	    }
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}

	return rootObject;
    }


    public Serializable flatScatter(ColGroup group,
				    Partitionable rootObject,
				    int root)
	    throws CCJException, NoSuchMemberException {
	// Do a flat tree like magpie
	// For some reason he starts with a FLATSCATTERTAG but then switches to SCATTERTAG when doing the sending and receiving ???

    if (DEBUG) System.out.println("Doing flatScatter");
    
	int local_tag = group.getNewGroupTag(ColGroup.FLATSCATTERTAG);

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("local rank " + rank);
	}
	if (root < 0 || root >= size) {
	    throw new NoSuchMemberException("root rank " + root);
	}

	if (rank != root) {
	    return receive(local_tag, root);
	}

	try {
	    for (int i = 0 ; i < size ; i++) {
		if (i != root) {
		    Serializable dataObject = rootObject.elementAt(i, size);
		    send(group, i, dataObject, root, local_tag);
		}
	    }
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
	return rootObject.elementAt(root, size);
    }


    public Serializable scatter(ColGroup group,
				Partitionable rootObject,
				int root)
	    throws CCJException, NoSuchMemberException {

	// Don't do a flat tree like magpie, but bcast everything and filter.

    if (DEBUG) System.out.println("Doing scatter");

	// NOTE: Used to be a Serializable [] -> manta dies !
	Serializable[] objectArray = null;
	int rank_rel;
	int mask, sum;
	Serializable result = null;
	Serializable data;

	// NOTE: Used to be a Serializable [] -> manta dies !
	Serializable[] objectsToSend;
	int local_tag = group.getNewGroupTag(ColGroup.SCATTERTAG);

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (rank == -1) {
	    throw new NoSuchMemberException("local rank " + rank);
	}
	if (root < 0 || root >= size) {
	    throw new NoSuchMemberException("root rank " + root);
	}

	if (size == 1) {
	    return rootObject.elementAt(0, 1);
	}

	if (DEBUG) {
	    System.out.println(rank + ": scatter started, root = " + root + " groupsize = "+ size + " counter << = " + local_tag);
	}

	if (rank == root) {
	    // Split the object in all the small objects
	    objectArray = new Serializable[size];

	    for (int i = 0 ; i < size ; i++) {
		if (i != root) {
		    objectArray[i] = rootObject.elementAt(i, size);
		}
	    }
	}

	rank_rel = rel_rank(rank, root, size);

	for (mask = 1; mask < size; mask *= 2);

	mask /= 2;

	sum = 0;

	try {
	    while (mask > 0) {
		if (sum + (mask & rank_rel) == rank_rel) {
		    /* do something in this step */

		    if ((mask & rank_rel) != 0) {
			if (DEBUG2) {
			    System.out.println("ColMember.broadcast: Member "+
					       rank + " calls receive from "+
					       abs_rank(rank_rel - mask, root, size));
			}

			data = receive(local_tag, abs_rank(rank_rel - mask, root, size));

			objectArray = (Serializable[]) data;

			result = objectArray[rank];
		    } else {
			if ((rank_rel + mask) < size) {
			    if (DEBUG2) {
				System.out.println("ColMember.broadcast: Member "+ rank + " sends to "+ abs_rank(rank_rel + mask, root , size));
			    }

			    int dst = abs_rank(rank_rel + mask, root, size);
			    int from = (dst < size ? dst : size - 1);
			    int till = ((dst + mask - 1) < size ? (dst + mask - 1) : size - 1);

			    objectsToSend = new Serializable[size];

			    System.arraycopy(objectArray, from, objectsToSend, from, till - from + 1);

			    send(group, abs_rank(rank_rel + mask, root, size), objectsToSend, rank, local_tag);
			}
		    }
		}

		sum += (mask & rank_rel);
		mask /= 2;
	    }
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}

	if (DEBUG) {
	    System.out.println(rank + ": scatter done.");
	}

	if (rank == root) {
	    result = rootObject.elementAt(root, size);
	}

	return result;
    }


    private int barrier_count;

    public void barrier(ColGroup group)
	    throws CCJException, NoSuchMemberException {

    if (DEBUG) System.out.println("Doing barrier");
    
	MyEntry e_left = null;
	MyEntry e_right = null;

	int local_tag = group.getNewGroupTag(ColGroup.BARRIERTAG);

	if (cachedGroup != group) {
	    cachedRank = group.getRank(this);
	    cachedGroup = group;
	    cachedSize  = group.size();
	}
	int rank = cachedRank;
	int size = cachedSize;

	if (size == 1) {
	    return;
	}

	if (DEBUG) {
	    System.out.println(cachedRank + ": Barrier started.");
	}

	try {
	    // Wait for left child
	    if (2 * rank + 1 < size) {
		e_left = mqueue.rendezVousGetStart(local_tag, 2 * rank + 1);
		if (ASSERT && e_left == null) {
		    System.err.println(cachedRank + ": e_left is null");
		}
	    }
	    // Wait for right child
	    if (2 * rank + 2 < size) {
		e_right = mqueue.rendezVousGetStart(local_tag, 2 * rank + 2);
		if (ASSERT && e_right == null) {
		    System.err.println(cachedRank + ": e_right is null");
		}
	    }

	    if (rank != 0) {
		int parent = (rank - 1) / 2;

		// Kick parent
		// Wait for parent
		MyListInterface dest = (MyListInterface)group.memberLists.elementAt(parent);
		dest.rendezVousAdd(local_tag, rank, null);
	    }

	    // Kick left child
	    if (e_left != null) {
		mqueue.rendezVousGetFinish(e_left, null);
	    }
	    // Kick right child
	    if (e_right != null) {
		mqueue.rendezVousGetFinish(e_right, null);
	    }

	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public void waitForCompletion() {
	if (!USE_THREADPOOL) {
	    return;
	}

	threadPool.waitForThreads();
    }

    
    public void setupRecords(int nrecvs) {
        this.record = true;
        this.numRecords = 0;
        this.records = new int[nrecvs];
    }
    
    public void endRecords() {
        this.record = false;
    }
    
    public int[] getRecords() {
        return this.records;
    }
    
    public int getMyRank(ColGroup group) throws NoSuchMemberException {
	return group.getRank(this);
    }


    private int rel_rank(int ABS, int ROOT, int SIZE) {
	return ((SIZE + ABS - ROOT) % SIZE);
    }


    private int abs_rank(int REL, int ROOT, int SIZE) {
	return ((REL + ROOT) % SIZE);
    }

}
