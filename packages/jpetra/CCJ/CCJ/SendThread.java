/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.rmi.RemoteException;
import java.io.Serializable;

final class SendThread extends java.lang.Thread {
    private static final boolean DEBUG = false;

    private Message m; // Cached message, to be filled with data to send
    private MyListInterface dst;// Reference to the list to place the data that is send
    private boolean filled; // Busy flag
    private ThreadPool pool; // Variables to let the threadpool know when this thread is done
    private int myId;
    private boolean quit; // Flag to indicate when to shutdown.
    private long startTime; // Debug variable

    SendThread(ThreadPool pool, int myId) {
	super();
	this.pool = pool;
	this.myId = myId;
	filled = false;
	quit = false;
	m = new Message();
    }


    protected synchronized void fill(int tag, int sender,
				     Serializable objectToSend,
				     MyListInterface dst, long startTime) {
	this.filled = true;
	m.tag = tag;
	m.source = sender;
	m.data = objectToSend;
	this.dst = dst;
	this.startTime = startTime;

	notify();
    }


    protected synchronized void done() {
	this.filled = true;
	this.quit = true;
	notify();
    }


    protected synchronized boolean send() throws CCJException {

	while (!filled) {
	    try {
		wait();
	    } catch (InterruptedException e) {
	    }
	}

	try {
	    if (quit) {
		return false;
	    }

	    if (DEBUG) {
		System.out.println("sendThread actually sending to: "+ dst);
	    }

	    dst.add(m.tag, m.source, m.data);

	    if (DEBUG) {
		System.out.println("sendThread has sent");
	    }
	} catch (RemoteException e) {
	    throw new CCJException("Send error, aborting..." + e);
	}

	filled = false;
	pool.threadDone(myId, startTime);
	return true;
    }


    public void run() {
	try {
	    while (send()) {
		// Just keep sending until done has been called.
	    }
	} catch (CCJException e) {
	    System.err.println("Catch CCJException " + e + "; call it quits...");
	}
    }

}
