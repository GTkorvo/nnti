/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.util.Vector;
import java.rmi.RemoteException;
import java.io.Serializable;

final public class ColGroup implements Serializable {

    // NOTE: Due to the counter[] in this object, it can only be used by a single GroupMember -- Jason
    private static final boolean DEBUG = false;

    // The tags to identify messages from different operations.
    // AllReduce and allGather have two phases, so need two tags.
    protected static final int BROADCASTTAG = 0;
    protected static final int REDUCETAG = 1;
    protected static final int GATHERTAG = 2;
    protected static final int SCATTERTAG = 3;
    protected static final int BARRIERTAG = 4;
    protected static final int ALLREDUCETAG1 = 5;
    protected static final int ALLREDUCETAG2 = 6;
    protected static final int ALLGATHERTAG = 7;
    protected static final int ALLGATHERPART2TAG = 8;
    protected static final int FLATGATHERTAG = 9;
    protected static final int FLATSCATTERTAG = 10;
    protected static final int UCASTTAG = 11;
    protected static final int MAXTAG = 11;

    // The next 2 constants are for creating tags. A tag is creating using:
    // <MESSAGECOUNTER>  |  <GROUPID> |  <OPERATIONTAG>
    //     16 bits           12 bits         4 bits // Totaling an int of 32 bits.

    // If these numbers are changed:
    // also change ColGroupCentral.MAXNUMBEROFGROUPS and
    // ColGroupCentral.GROUPSHIFT.

    // Message counter will be shifted COUNTERSHIFT  number of bits.
    private static final int COUNTERSHIFT = 16;
    // COUNTERMASK is 2^^15-1 Mask to circulate messagecounter.
    // Without the highest bit to prevent negative values
    private static final int COUNTERMASK = 0x7FFF;

    protected String name;
    protected boolean active;
    // protected
	public Vector memberLists;
    protected int id;

    private int group_counter[], ucast_send_counter[], ucast_receive_counter[];

    ColGroup(String name, int id) {
	this.name = name;
	this.id = id;
	memberLists = new Vector();
	group_counter = new int[MAXTAG + 1];
    }


    public boolean equals(ColGroup otherGroup) {
	return name.equals(otherGroup.name);
    }


    public int getRank(ColMember member) throws NoSuchMemberException {
	int i = memberLists.indexOf((MyListInterface) member.mqueue);

	if (DEBUG) {
	    System.out.println("ColGroup getRank looking for: " +
			       (MyListInterface)member.mqueue + " found: "+ i);

	    for (int j = 0 ; j < memberLists.size(); j++) {
		System.out.println("ColGroup member " + j + " = " +
				    (MyListInterface)memberLists.elementAt(j));
	    }
	}

	if (i == -1) {
	    throw new NoSuchMemberException();
	}

	return i;
    }


    public MyListInterface getMember(int rank) throws NoSuchMemberException {
	Object memberList;

	try {
	    memberList = memberLists.elementAt(rank);
	} catch (ArrayIndexOutOfBoundsException e) {
	    memberList = null;
	}

	if (DEBUG) {
	    System.out.println("ColGroup.getMember(" + rank + ") list: " + memberList);
	}

	if (memberList == null) {
	    throw new NoSuchMemberException();
	}

	return (MyListInterface) memberList;
    }


    public int size() {
	int maxSize = memberLists.size();
	int size = 0;

	for (int i = 0 ; i < maxSize ; i++) {
	    if (memberLists.elementAt(i) != null) {
		size++;
	    }
	}

	return size;
    }


    protected void addMember(MyListInterface memberList) {
	memberLists.addElement(memberList);
    }


    void deleteMember(MyListInterface memberList) {
	String hostName = null;
	try {
	    hostName = java.net.InetAddress.getLocalHost().getHostName();
	} catch (java.net.UnknownHostException e) {
	    System.err.println("Cannot get local Hostname " + e);
	}
	memberLists.remove(memberList);
    }


    protected void addMember(MyListInterface memberList, int rank)
	    throws RankAlreadyInUseException {
	memberLists.ensureCapacity(rank);

	if (memberLists.elementAt(rank) != null) {
	    throw new RankAlreadyInUseException(name + ": rank = " + rank);
	}

	memberLists.add(rank, memberList);
    }


    protected boolean isActive() {
	return active;
    }


    protected void setActive() {
	active = true;
	ucast_send_counter = new int[memberLists.size()];
	ucast_receive_counter = new int[memberLists.size()];
    }


    protected void addAll(ColGroup other) {
	for (int i = 0 ; i < other.size(); i++) {
System.err.println("from addAll: add this queue " + other.memberLists.elementAt(i) + " to my queue " + memberLists);
	    memberLists.addElement(other.memberLists.elementAt(i));
	}
    }


    protected int getNewGroupTag(int tag) throws CCJException {

	if (tag == UCASTTAG) {
	    throw new CCJException("Internal error: Can not create a tag for a non-group operation!");
	}

	int temp = group_counter[tag] = (group_counter[tag] + 1) & COUNTERMASK;
	temp = temp << COUNTERSHIFT;

	return (id + temp + tag);
    }


    protected int getNewUnicastSendTag(int dest) {

	int temp = ucast_send_counter[dest] = (ucast_send_counter[dest] + 1) & COUNTERMASK;
	temp = temp << COUNTERSHIFT;

	return (id + temp + UCASTTAG);
    }


    protected int getNextUnicastReceiveTag(int source) {

	int temp = ucast_receive_counter[source] = (ucast_receive_counter[source] + 1) & COUNTERMASK;
	temp = temp << COUNTERSHIFT;

	return (id + temp + UCASTTAG);
    }

}
