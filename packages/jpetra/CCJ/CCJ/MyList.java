/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.rmi.server.UnicastRemoteObject;
import java.io.Serializable;

final class MyList extends UnicastRemoteObject implements MyListInterface {
    // Message list with a node queue with at least 1 Entry ready to be used

    private static final boolean DEBUG = true;
    private static final boolean ASSERT = false;
    protected static int BUCKETS = 2; // Must be a power of 2
    protected static int STARTING_POOLSIZE = 1;

    private MyEntry head[];
    private MyEntry tail[];
    private MyEntry pool;
    private ColMember member;
    // Debug stats
    protected int numberOfNext;

    protected MyEntry future;

    private static final int WAIT_FOR_SET = 0;
    private static final int WAIT_FOR_ADD = 1;
    private static final int WAIT_FOR_GET = 2;

    protected int future_state = WAIT_FOR_SET;

    MyList(ColMember member, int BUCKETS, int STARTING_POOLSIZE)
	    throws java.rmi.RemoteException {
	super();

	this.BUCKETS = BUCKETS;
	this.STARTING_POOLSIZE = STARTING_POOLSIZE;
	head = new MyEntry[BUCKETS];
	tail = new MyEntry[BUCKETS];

	for (int i = 0 ; i < BUCKETS ; i++) {
	    head[i] = tail[i] = null;
	}

	pool = createPool();
	this.member = member;
	// Init debug variable
	numberOfNext = 0;

	future = new MyEntry();
    }


    private int hash(int tag, int source) {
	return (tag + source) & (BUCKETS - 1);
    }


    private MyEntry createPool() {
	MyEntry e, head;
	head = new MyEntry();

	for (int i = 0 ; i < STARTING_POOLSIZE ; i++) {
	    e = new MyEntry();
	    e.next = head;
	    head = e;
	}

	return head;
    }


    private int size(int b) {
	int size = 0;

	for (MyEntry e = head[b]; e != null ; size++, e = e.next);

	return size;
    }


    protected int size() {
	int size = 0;

	for (int b = 0 ; b < BUCKETS ; b++) {
	    size += size(b);
	}

	return size;
    }


    protected int poolSize() {
	MyEntry e = pool;
	int i;

	for (i = 0 ; e != null ; i++, e = e.next);
	return i;
    }


    private MyEntry poolGet() {
	// Must be synchronized (this)
	MyEntry e = pool;

	if (e == null) {
	    e = new MyEntry();
	} else {
	    pool = pool.next;
	}
	if (ASSERT) {
	    e.in_use = true;
	}

	return e;
    }


    private void poolPut(MyEntry e) throws CCJException {
	// Must be synchronized (this)
	e.next = pool;
	if (ASSERT) {
	    if (! e.in_use) {
		throw new CCJException("Internal error: Clear an entry that's already cleared by " + e.clearer);
	    }
	    e.in_use = false;
	    if (e.received) {
		throw new CCJException("Internal error: Clear an entry with received = true");
	    }
	}
	pool = e;
    }


    private MyEntry getFirstEntry(int tag, int source) throws CCJException {
	int bucket = hash(tag, source);

	if ((ASSERT && member.cachedRank() == 0 && bucket == 6) || DEBUG) {
	    System.out.println(member.cachedRank() + ": MyList looking for message, tag = "+
			       tag + " source = " + source + " size =" + size() +
			       " in bucket[" + bucket + "] = " + head[bucket] +
			       " length " + size(bucket));
	}

	for (;;) {
	    MyEntry prev = null;
	    for (MyEntry e = head[bucket]; e != null; e = e.next) {
		if (ASSERT) {
		    if (! e.in_use) {
			throw new CCJException("Internal error: In queue: entry that's already cleared by " + e.clearer);
		    }
		}

		if (e.source == source && e.tag == tag) {
		    if (prev == null) {
			head[bucket] = e.next;
		    } else {
			prev.next = e.next;
		    }

		    if (tail[bucket] == e) {
			tail[bucket] = prev;
		    }

		    if ((ASSERT && member.cachedRank() == 0 && bucket == 6) || DEBUG) {
			System.out.println(member.cachedRank() + ": Found a message " + e +
					    " bucket " + bucket +
					    " head[] = " + head[bucket] +
					    "; length := " + size(bucket));
		    }

		    if (e.confirm) {
			synchronized (e) {
			    e.received = true;
			    e.notifyAll();
			}
		    }

		    return e;
		}

		// Debug var update
		numberOfNext++;
		prev = e;
	    }

	    try {
		wait();
	    } catch (InterruptedException someException) {
		/* poll and maybe go on waiting */
	    }
	}
    }


    synchronized Serializable getFirst(int tag, int source)
	    throws CCJException {
	MyEntry e = getFirstEntry(tag, source);
	Serializable data = e.data;

	if (! e.confirm) {
	    if (ASSERT) {
		e.clearer = "getFirst e.confirm = " + e.confirm;
	    }
	    poolPut(e);
	}

	return data;
    }


    private synchronized MyEntry add(int tag, int source, Serializable data,
				     boolean confirm) {
	int bucket = hash(tag, source);

	MyEntry e;

	e = poolGet();
	e.confirm = confirm;
	e.data    = data;
	e.tag     = tag;
	e.source  = source;
	e.next = null;

	if (head[bucket] == null) {
	    head[bucket] = e;
	} else {
	    tail[bucket].next = e;
	}
	tail[bucket] = e;

	if ((ASSERT && member.cachedRank() == 0 && bucket == 6) || DEBUG) {
	    System.out.println(member.cachedRank() + ": MyList.add called, " +
			       e + " size = " + size() +
			       " bucket[" + bucket + "] " + head[bucket] +
			       " length " + size(bucket) + " entry " + e);
	}

	notifyAll();

	return e;
    }


    public void add(int tag, int source, Serializable data)
	    throws java.rmi.RemoteException {
	add(tag, source, data, false /* no need to confirm */);
    }


    public void confirmedAdd(int tag, int source, Serializable data)
	    throws java.rmi.RemoteException, CCJException {
	if (DEBUG) {
	    System.out.println(member.cachedRank() + ": MyList.confirmedAdd called, tag = " +
			       tag + " source = " + source + " size = " +
			       size() + " Data = " + data);
	}

	MyEntry e = add(tag, source, data, true /* need to confirm */);

	synchronized (e) {
	    while (! e.received) {
		try {
		    e.wait();
		} catch (InterruptedException exc) {
		    if (DEBUG) {
			System.out.println("conFirmedAdd: Entry woken up, tag = " +
					   tag + " source = " + source + " size = " +
					   size() + " Data = " + data +
					   " Exception = " + exc.toString());
		    }
		}
	    }
	}

	synchronized (this) {
	    if (ASSERT) {
		e.clearer = "confirmedAdd";
	    }
	    e.received = false;
	    e.confirm = false;
	    poolPut(e);
	}

	if (DEBUG) {
	    System.out.println("conFirmedAdd: Entry woken up, tag = " +
			       tag + " source = " + source + " size = " +
			       size() + " Data = " + data);
	}
    }


    private int filled_dims(Serializable data) {
	int n = 0;
	if (data == null) {
	    return -1;
	}

	try {
	    Object[] b = (Object[])data;
	    for (int i = 0; i < b.length; i++) {
		if (b[i] != null) {
		    n++;
		}
	    }
	    return n;
	} catch (ClassCastException e) {
	    return 0;
	}
    }


    public Serializable rendezVousAdd(int tag, int source, Serializable data)
	    throws java.rmi.RemoteException, CCJException {

	MyEntry e;

	synchronized (this) {

	    if (DEBUG) {
		System.out.println(member.cachedRank() + ": MyList.rendezVousAdd called, tag = " +
				   tag + " source = " + source + " size = " +
				   size() + " Data = " + data);
	    }

	    e = add(tag, source, data, false /* no need to confirm */);
	    if (ASSERT) {
		if (e.received) {
		    System.err.println("Seem to heave synced/received before I'm here");
		}
		if (e.reply != null) {
		    System.err.println("Seem to heave synced/reply before I'm here");
		}
	    }
	}

	Serializable reply = null;

	synchronized (e) {
	    if (ASSERT) {
		if (! e.in_use) {
		    throw new CCJException("Internal error: add gave me an entry that's already cleared by " + e.clearer);
		}
	    }

	    while (! e.received) {
		try {
		    e.wait();
		} catch (InterruptedException exc) {
		    System.out.println("rendezVousAdd: Entry woken up, tag = " +
				       tag + " source = " + source + " size = " +
				       size() + " Data = " + data +
				       " Exception = " + exc.toString());
		}
	    }

	    if (ASSERT) {
		if (! e.in_use) {
		    throw new CCJException("Internal error: after wait in rendezVousAdd: an entry that's already cleared by " + e.clearer);
		}
	    }

	    reply = e.reply;
	}

	synchronized (this) {
	    // Clear the ListEntry and put it into the pool
	    if (ASSERT) {
		e.clearer = "rendezVousAdd";
	    }
	    e.received = false;
	    poolPut(e);
	}

	if (DEBUG) {
	    System.out.println("rendezVousAdd: Entry woken up, tag = " +
			       tag + " source = " + source + " size = " +
			       size() + " Data = " + data);
	}

	return reply;
    }


    Serializable rendezVousGet(int tag, int source, Serializable reply)
	    throws CCJException {
	MyEntry e = rendezVousGetStart(tag, source);
	return rendezVousGetFinish(e, reply);
    }


    MyEntry rendezVousGetStart(int tag, int source) throws CCJException {

	MyEntry e;

	synchronized (this) {
	    if (DEBUG) {
		System.out.println(member.cachedRank() + ": MyList looking for message, tag = " +
				   tag + " source = " + source + " size =" + size());
	    }

	    e = getFirstEntry(tag, source);
	}

	return e;
    }


    Serializable rendezVousGetFinish(MyEntry e, Serializable reply) {

	Serializable data;

	synchronized (e) {
	    e.reply = reply;
	    e.received = true;
	    e.notifyAll();
	    data = e.data;
	}

	return data;
    }


    void futureSet(Serializable reply) {

	// The idea of this is that we can post the reply to
	// an RMI ahead of time, so the CPU can do something
	// useful (like doing an RMI) while waiting.

//        System.out.println("futureSet saving " + reply);

	synchronized (future) {

	    while (future_state != WAIT_FOR_SET) {
		try {
		    future.wait();
		} catch (InterruptedException e) {
		    System.out.println("futureGet got exception " + e);
		}
	    }

	    future.reply = reply;
//            future.received = false;
	    future_state = WAIT_FOR_ADD;
	    future.notifyAll();
	}
    }


    protected Serializable futureGet() {

	// Use this to retrieve the data sent by the RMI that ate your reply.
	Serializable data = null;

	synchronized (future) {

	    while (future_state != WAIT_FOR_GET) {
		try {
		    future.wait();
		} catch (InterruptedException e) {
		    System.out.println("futureGet got exception " + e);
		}
	    }

	    data = future.data;
//            future_used = false;
	    future_state = WAIT_FOR_SET;
	    future.notifyAll();
	}
//        System.out.println("futureGet returning " + data);

	return data;
    }


    public Serializable futureAdd(Serializable data)
	    throws java.rmi.RemoteException {

	synchronized (future) {

	    while (future_state != WAIT_FOR_ADD) {
		// We're early -> wait for the result to appear.
		try {
		    future.wait();
		} catch (InterruptedException e) {
		    System.out.println("futureGet got exception " + e);
		}
	    }

//            System.out.println("futureAdd saving " + data);
	    future.data = data;
//            future.received = true;
	    future_state = WAIT_FOR_GET;
	    future.notifyAll();

//            System.out.println("futureAdd returning " + future.reply);
	    return future.reply;
	}
    }

}
