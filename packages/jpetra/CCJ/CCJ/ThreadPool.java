/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.rmi.RemoteException;
import java.io.Serializable;

final class ThreadPool {

    private static final boolean DEBUG = false;
    private static final boolean DO_TIMINGS = false; // Also update ColMember.DO_TIMINGS
    private RemoteException lastException;
    protected static int numberOfThreads;
    private int nextThread;
    private int numberOfActiveThreads;

    private SendThread[] pool;
    private boolean[] active;

    // Debug timing variables and statistics
    protected long sendTime;
    protected int hits, misses;
    protected int numberOfWaits;
    protected long waitTime;
    protected long searchTime;

    ThreadPool(int numberOfThreads) {
	// Initialize debug timing variables
	hits = 0;
	misses = 0;
	numberOfActiveThreads = 0;
	sendTime = 0;
	waitTime = 0;
	numberOfWaits = 0;

	// Init cache variables
	nextThread = 0;
	this.numberOfThreads = numberOfThreads;
	// Init sendthreadpool with numberOfThread threads
	pool = new SendThread[numberOfThreads];
	active = new boolean[numberOfThreads];

	for (int i = 0 ; i < numberOfThreads ; i++) {
	    pool[i] = new SendThread(this, i);
	    active[i] = false;
	    pool[i].start();
	}
    }


    protected synchronized void send(int tag,
				     int sender,
				     Serializable objectToSend,
				     MyListInterface dst,
				     long startTime)
	    throws CCJException {
	// Debug timing variables
	long st = 0, et;

	if (DEBUG) {
	    System.out.println("ThreadPool.send called");
	}

	for (;;) {
	    if (numberOfActiveThreads < numberOfThreads) {
		if (DO_TIMINGS) {
		    st = System.currentTimeMillis();
		}

		// Find a sendThread that is not busy yet
		// Try the next thread it may be idle
		int start;

		if (!active[nextThread]) {
		    start = nextThread;
		    hits++;
		} else {
		    start = 0;
		    misses++;
		}

		// Find an idle thread, starting with the most likely candidate.
		int t = start;

		for (int i = 0 ; i < numberOfThreads ; i++) {
		    if (!active[t]) {
			active[t] = true;
			numberOfActiveThreads++;
			nextThread = (t + 1) % numberOfThreads;
			pool[t].fill(tag, sender, objectToSend, dst,
				 startTime);
			// Update debug info
			if (DO_TIMINGS) {
			    et = System.currentTimeMillis();
			    searchTime += et - st;
			}

			return;
		    }

		    t = (t+1) % numberOfThreads;
		}

		throw new CCJException("Internal error: no pool thread available");
	    }

	    try {
		// No thread available waiting for a thread.
		if (DO_TIMINGS) {
		    st = System.currentTimeMillis();
		}

		wait();

		if (DO_TIMINGS) {
		    et = System.currentTimeMillis();
		    waitTime += et - st;
		}
	    } catch (InterruptedException e) {}

	}
    }


    synchronized void threadDone(int i, long startTime) {

	// The thread has completed the RMI, and is now available.

	// Update stats, thread is available
	active[i] = false;

	numberOfActiveThreads--;

	if (DO_TIMINGS) {
	    long endTime = System.currentTimeMillis();
	    sendTime += endTime - startTime;
	}

	// Wake up, a thread is available
	notify();
    }


    protected synchronized RemoteException waitForThreads() {
	while (numberOfActiveThreads > 0) {
	    try {
		wait();
	    } catch (InterruptedException e) {
		// Wake up!
	    }

	}

	RemoteException e = lastException;

	lastException = null;
	return e;
    }


    protected void done() {
	// The system should shut down.
	// Check whether it is ok to shutdown, no operations pending?
	if (numberOfActiveThreads > 0) {
	    System.err.println("WARNING Threads are still active, during threadpool.done!");
	}

	// Tell all sendthreads to shutdown.
	for (int i = 0 ; i < numberOfThreads ; i++) {
	    pool[i].done();

	    try {
		pool[i].join(500); // Wait at most 500 msec for thread to end.
	    }
	    catch (InterruptedException e) {
		System.err.println("SendThread number "+ i + " of threadpool, did not die, "+
				   "during shutdown.");
	    }
	}
    }

}
