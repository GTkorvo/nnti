/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.rmi.registry.Registry;
import java.rmi.registry.LocateRegistry;
import java.rmi.RemoteException;
import java.rmi.Naming;
import java.net.*;
import java.io.DataInputStream;
// import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.DataOutputStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;

final public class ColGroupMaster {

    private static final boolean DEBUG = true;

    private static final int MAX_RETRIES = 10;

    private static int currentId = 0;
    private static Vector host_names = new Vector();
    private static Vector host_ports = new Vector();    
    protected static ColGroupCentralInterface groupCentral = null;
    private static int numberOfLocalMembers = 0;
    private static int numberOfMembersDone = 0;
    private int host;
    private static int port;
    private int numberOfCPUs;

    private static String hostName;
    private String masterName;
    private String masterPort;
    private CcjHosts myHosts;
    
    private Registry local;

    public static String hostPort () {
        return hostName + ":" + port;
    }
    
    public ColGroupMaster(String filePath)
	    throws  CCJException,
		    GroupAlreadyActiveException {

    /*
     * The Vnode network reconition code has essentially been completly rewritten from the orginal.
     * It is now handled by the CcjHosts class.
     */
    
    myHosts = new CcjHosts(filePath, host_names, host_ports);
    
    hostName = myHosts.getMyHost();
    port = myHosts.getMyPort();
    host = myHosts.getMyHostNum();
    numberOfCPUs = myHosts.getNumCpus();
    
    System.out.println("hostName: "+hostName+" port: "+port);
    
    /*END ALL MODIFICATIONS*/
	if (host == -1) {
	    throw new CCJException("Cannot find my host name " + hostName + " in list");
	}

	if (DEBUG) {
	    System.out.println("ColGroupMaster created on CPU number "+ host + " named: "+hostName);
	}

	try {
	    startRMIenv();
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException("Failed to create registry : " + e);
	}

	masterName = (String)host_names.get(0);
    masterPort = (String)host_ports.get(0);
	try {
	    getCentral();
	} catch (java.net.MalformedURLException e) {
	    throw new CCJException("Not a well-formed URL: \"GroupCentral\": " + e);
	} catch (java.rmi.AlreadyBoundException e) {
	    throw new CCJException("AlreadyBoundException " + e);
	} catch (java.rmi.AccessException e) {
	    throw new CCJException("AccessException " + e);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException("RemoteException " + e);
	}
    }


    public int getNumberOfCpus() {
	    return numberOfCPUs;
    }

    public static String getHostPort (int id) {
        return (String)host_ports.get(id);
    }

    public static String getHostName (int id) {
        return (String)host_names.get(id);
    }
    
    public void addMember(String groupName, ColMember member)
	    throws CCJException, GroupAlreadyActiveException {
	try {
	    groupCentral.addMember(groupName, (MyListInterface) member.mqueue);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public void addMember(String groupName, int rank, ColMember member)
	    throws  CCJException,
		    GroupAlreadyActiveException,
		    RankAlreadyInUseException {
	try {
	    groupCentral.addMember(groupName, rank, (MyListInterface) member.mqueue);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public void deleteMember(String groupName, ColMember member)
	    throws  CCJException,
		    NoSuchMemberException {
	try {
	    groupCentral.deleteMember(groupName, member.mqueue);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
	member.done();
    }


    public ColGroup getGroup(String groupName, int numberOfMembers)
	    throws CCJException, NoSuchGroupException {
	try {
	    return groupCentral.getGroup(groupName, numberOfMembers);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public void copyGroup(String destName, String sourceName)
	    throws  CCJException,
		    NoSuchGroupException,
		    GroupAlreadyActiveException {

	try {
	    groupCentral.copyGroup(destName, sourceName);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    public void deleteGroup(String groupName)
	    throws CCJException, NoSuchGroupException {
	try {
System.err.println("delete group " + groupName);
	    groupCentral.deleteGroup(groupName);
	} catch (java.rmi.RemoteException e) {
	    throw new CCJException(e.toString());
	}
    }


    private void getCentral()
	    throws  CCJException,
		    java.net.MalformedURLException,
		    java.rmi.AccessException,
		    java.rmi.AlreadyBoundException,
		    java.rmi.RemoteException {
	// Try to read the hostname of the master from the shared file
	int	sleep = 1000;
	int	i;

	if (host == 0) {
	    groupCentral = new ColGroupCentral(hostName+":"+port);
	}

	for (i = 0; i < MAX_RETRIES; i++) {

	    if (DEBUG) {
		System.out.println("Naming.lookup(" + masterName + ":" + masterPort + "/GroupCentral)");
	    }

	    try {
		groupCentral = (ColGroupCentralInterface)Naming.lookup("//" + masterName  + ":" + masterPort + "/GroupCentral");
		break;
	    } catch (java.rmi.NotBoundException e) {
		System.err.println("Retrying lookup of GroupCentral...");
	    } catch (java.rmi.RemoteException e) {
		System.err.println("Retrying lookup of GroupCentral...");
	    }

	    try {
		Thread.sleep(sleep);
		sleep += 1000;
	    } catch (InterruptedException e) {
		// Ignore
	    }
	}

	if (i == MAX_RETRIES) {
	    throw new CCJException(host + " lookup timed out");
	}

	if (DEBUG) {
	    System.out.println("Host " + host + " done lookup");
	}
    }


    private void startRMIenv() throws java.rmi.RemoteException {
	// Start the registry.
	try {
	    // First try to create a registry object
	    local = LocateRegistry.createRegistry(port);
	} catch (java.rmi.RemoteException e) {
	    //added to support multiply Vnodes on the same IP by using ports
	    port = myHosts.findNextHostPort();
	    if (port != 0) {startRMIenv();}
	    //from original
	    // Didn't work; registry allready present ?
	    //local = LocateRegistry.getRegistry();
	}

	// Start a security manager.
	// We should try do figure out when to install this thing .....
	//        System.setSecurityManager(new RMISecurityManager());
    }


    public static String getId() {
	return new String(Integer.toString(currentId++));
    }

}
