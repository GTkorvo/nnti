/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

public interface ColGroupCentralInterface extends java.rmi.Remote {

    public void addMember(String groupName, MyListInterface member)
		    throws  java.rmi.RemoteException,
			    CCJException,
			    GroupAlreadyActiveException;

    public void addMember(String groupName, int rank, MyListInterface member)
		    throws  java.rmi.RemoteException,
			    CCJException,
			    GroupAlreadyActiveException,
			    RankAlreadyInUseException;

    public void deleteMember(String groupName, MyListInterface member)
		    throws  java.rmi.RemoteException,
			    CCJException;

    public ColGroup getGroup(String groupName, int numberOfMembers)
		    throws  java.rmi.RemoteException,
			    NoSuchGroupException;

    public void copyGroup(String destName, String sourceName)
		    throws  java.rmi.RemoteException,
			    CCJException,
			    NoSuchGroupException,
			    GroupAlreadyActiveException;

    public void deleteGroup(String groupName)
		    throws  java.rmi.RemoteException,
			    CCJException,
			    NoSuchGroupException;
}
