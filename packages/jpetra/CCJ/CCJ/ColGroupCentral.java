/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.rmi.server.UnicastRemoteObject;
import java.rmi.Naming;
import java.util.Vector;
import java.rmi.RemoteException;

final public class ColGroupCentral
	extends UnicastRemoteObject
	implements ColGroupCentralInterface {

    private static final boolean DEBUG = false;

    // If these numbers are changed: see ColMember.java
    private static final int GROUPSHIFT = 4; // Number of bits groups will be shifted into tags.
    // This is 2log (ColMember.MAXTAG) rounded up.
    private static final int MAXNUMBEROFGROUPS = 4095; // 2^^(16-GROUPSHIFT) - 1

    private Vector groups;
    private static int currentId = 0;

    ColGroupCentral(String host)
	    throws  java.net.MalformedURLException,
		    java.rmi.AlreadyBoundException,
		    java.rmi.RemoteException {
	groups = new Vector();
	Naming.bind("//" + host +"/GroupCentral", this);
	if (DEBUG) {
	    System.out.println("GroupCentral bound");
	}
    }


    public synchronized void addMember(String groupName, MyListInterface member)
	    throws java.rmi.RemoteException, CCJException, GroupAlreadyActiveException {

	if (DEBUG) {
	    System.out.println("central.addmember " + member + " to group " + groupName);
	}

	for (int i = 0 ; i < groups.size(); i++) {
	    ColGroup group = (ColGroup) groups.elementAt(i);

	    if (group.name.equals(groupName)) {
		if (DEBUG) {
		    System.out.println("Found group "+groupName);
		}

		if (group.isActive()) {
		    throw new GroupAlreadyActiveException("Group name = " + groupName);
		}

		group.addMember(member);
		notifyAll();
		return;
	    }
	}

	// Group does not yet exist
	if (DEBUG) {
	    System.out.println("GroupCentral: group "+groupName + " does not yet exist,"+
			       "creating");
	}

	if (groups.size() >= MAXNUMBEROFGROUPS) {
	    throw new CCJException("Too many groups created. Only " +
				   MAXNUMBEROFGROUPS + " permitted.");
	}

	ColGroup group =
	    new ColGroup(groupName, groups.size() << GROUPSHIFT);
	group.addMember(member);
	groups.addElement(group);
    }


    public void addMember(String groupName, int rank, MyListInterface member)
	    throws  java.rmi.RemoteException,
		    CCJException,
		    GroupAlreadyActiveException,
		    RankAlreadyInUseException {

	if (DEBUG) {
	    System.out.println("central.addmember to group named: "+
			       groupName + " with specified rank: " + rank);
	}

	for (int i = 0 ; i < groups.size(); i++) {
	    ColGroup group = (ColGroup) groups.elementAt(i);

	    if (group.name.equals(groupName)) {
		if (DEBUG) {
		    System.out.println("Found group "+groupName);
		}

		if (group.isActive()) {
		    throw new GroupAlreadyActiveException("Group name = " + groupName);
		}

		group.addMember(member, rank);
		notifyAll();
		return;
	    }
	}

	// Group does not yet exist
	if (DEBUG) {
	    System.out.println("GroupCentral: group "+groupName +
			       " does not yet exist, creating");
	}

	if (groups.size() >= MAXNUMBEROFGROUPS) {
	    throw new CCJException("Too many groups created. Only " +
				   MAXNUMBEROFGROUPS + " permitted.");
	}

	ColGroup group = new ColGroup(groupName, groups.size() << GROUPSHIFT);
	group.addMember(member);
	groups.addElement(group);
    }


    public synchronized void deleteMember(String groupName, MyListInterface member)
	    throws  java.rmi.RemoteException,
		    CCJException {
	ColGroup group = locateGroup(groupName);
	group.deleteMember(member);
	if (group.size() == 0) {
	    deleteGroup(groupName);
	}
    }


    private ColGroup locateGroup(String groupName) throws NoSuchGroupException {
	for (int i = 0 ; i < groups.size(); i++) {
	    ColGroup group = (ColGroup) groups.elementAt(i);

	    if (group.name.equals(groupName)) {
		return group;
	    }
	}

	throw new NoSuchGroupException("Cannot find group " + groupName);
    }


    public synchronized ColGroup getGroup(String groupName, int numberOfMembers)
	    throws java.rmi.RemoteException, NoSuchGroupException {

	if (DEBUG) {
	    System.out.println("GroupCentral: getGroup ("+ groupName + " , "+ numberOfMembers + ")");
	}

	for (int i = 0 ; i < groups.size(); i++) {
	    ColGroup group = (ColGroup) groups.elementAt(i);

	    if (group.name.equals(groupName)) {
		// have all members been added to the group yet?
		while (group.size() < numberOfMembers) {
		    try {
			wait();
		    } catch (InterruptedException e) {
			// Do nothing
		    }

		}

		// Enough members have been added

		if (group.size() > numberOfMembers) {
		    System.err.println("Warning group "+groupName + " has more members "+"than it should have, during ColGroupCentral.getGroup");
		}

		group.setActive(); // The group will now be active.
		return group;
	    }
	}

	throw new NoSuchGroupException("In ColGroupCentral.getGroup");
    }


    public synchronized void copyGroup(String destName, String sourceName)
	    throws  java.rmi.RemoteException,
		    CCJException,
		    NoSuchGroupException,
		    GroupAlreadyActiveException {
	if (groups.size() >= MAXNUMBEROFGROUPS) {
	    throw new CCJException("Too many groups created. Only " +
				   MAXNUMBEROFGROUPS + " permitted.");
	}

	// Check to see if the destination name is already in use.
	for (int i = 0 ; i < groups.size(); i++) {
	    ColGroup group = (ColGroup) groups.elementAt(i);

	    if (group.name.equals(destName)) {
		throw new GroupAlreadyActiveException("copyGroup: Destination " + destName + " already in use");
	    }
	}

	for (int i = 0 ; i < groups.size(); i++) {
	    ColGroup newGroup;
	    ColGroup group = (ColGroup) groups.elementAt(i);

	    if (group.name.equals(sourceName)) {
		// The source group is found: copy it.
		newGroup = new ColGroup(destName, groups.size() << GROUPSHIFT);
		newGroup.addAll(group);
		groups.addElement(newGroup);
		return;
	    }
	}

	throw new NoSuchGroupException("ColGroupCentral.copyGroup: " +
				       sourceName);
    }


    synchronized public void deleteGroup(String groupName)
	    throws  java.rmi.RemoteException,
		    CCJException,
		    NoSuchGroupException {
	ColGroup group = locateGroup(groupName);
	groups.removeElement(group);
	if (groups.size() == 0) {
	    try {
		Naming.unbind("GroupCentral");
	    } catch (java.net.MalformedURLException e) {
		throw new CCJException("Seems group master name is malformed?? " + e);
	    } catch (java.rmi.NotBoundException e) {
		throw new CCJException("Group Master already unbound " + e);
	    }
	}
    }

}
