/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.io.Serializable;

interface MyListInterface extends java.rmi.Remote {
    // These methods must be public due to RMI
    public void add(int tag, int source, Serializable data)
			    throws java.rmi.RemoteException;
    public void confirmedAdd(int tag, int source, Serializable data)
			    throws java.rmi.RemoteException, CCJException;
    public Serializable rendezVousAdd(int tag, int source, Serializable data)
			    throws java.rmi.RemoteException, CCJException;
    public Serializable futureAdd(Serializable m)
			    throws java.rmi.RemoteException;
}
