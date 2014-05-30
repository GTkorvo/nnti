// @HEADER
// ***********************************************************************
//
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

package CCJ;

import java.util.Vector;
import java.io.*;
import java.net.*;
import java.rmi.registry.Registry;

/**
 * Reads a file from the file name path passed to it, and then reads in a host per line in the format:<br>
 * <code>host:port</code><br>
 * The port number will be added to the default port value specified by
 * <code>import java.rmi.registry.Registry;</code>.
 *<p>
 * To deal with the fact that InetAddress.getLocalHost() has some very odd behavior
 * and the fact that machines may have more than one valid IP address
 * this modification uses NetworkInterface.getByInetAddress to determine if
 * the IP address specified is bound to any local network interface.
 * Effectivly this means that any local interface or IP can be used.
 *<p>
 * The terminology here is a mix of that from Jpetra and CCJ.
 *<p>
 * NOTE/BUGS: Additional error checking should be included to eliminate the ability
 * to specify loopback interfaces as valid IPs.  Also, the first IP found to be bound
 * to an interface will be used.  This will prevent multiple processes on the same
 * physical node from using multiple IP addresses.
 *
 * @author  Jason Cross
 */

public class CcjHosts {
    
    /**
     * all the IPs of the vnodes
     */
    private String[] hosts;
    
    /**
     * all the ports for the IPs of the vnodes
     */
    private String[] ports;
    
    /**
     * this vnode's IP; may not initially be the correct IP
     */
    private String hostName;
    
    /**
     * this vnode's port; may not initially be the correct port
     */
    private int port;
    
    /**
     * this vnode's rank
     */
    private int host;
    
    /**
     * number of tatal vnodes
     */
    private int numberOfCPUs;
    
    /**
     * all the IPs of the vnodes; is a reference to the static host_names in ColGroupMaster
     */
    private Vector host_names;
    
    /**
     * all the ports of the vnodes; is a reference to the static host_ports in ColGroupMaster
     */
    private Vector host_ports;
    
    
    /**
     * Determines the IP addresses of all vnodes and identifies this vnodes' IP.
     *
     * @param filePath path and file name to the ccjhosts file
     * @param host_names static reference passed from ColGroupMaster
     * @param host_ports static reference passed from ColGroupMaster
     */
    public CcjHosts(String filePath, Vector host_names, Vector host_ports) throws CCJException {
        this.host_names = host_names;
        this.host_ports = host_ports;
        
        //try to open and read one line from the given ccjhosts file
        BufferedReader in = null;
        String line = null;
        try {
            //creates a BufferedReader to read the file with for added effeciency and for
            //the <code>readLine()</code> method
            in = new BufferedReader(new FileReader(filePath));
            
            line = in.readLine();
        }
        catch (Exception e) {
            throw new CCJException("An error occurred while reading from the CcjHosts file: " + e);
        }
        
        //keeps track of number of while loops excuted; used to determine host
        int i = 0;
        NetworkInterface activeInterface = null;
        host=-1;
        numberOfCPUs = 0;
        
        //goes through each line of the ccjhosts file, parses it, and then addes the vnode
        //IP and port to host_names and host_ports
        while (line != null) {
            
            //checks for comments begining with / or #
            if (('/' != line.charAt(0)) && ('#' != line.charAt(0))) {
                //seperate the IP and the port
                String[] lineSplit = line.split(":");
                //set the initial port
                int hostPort=Registry.REGISTRY_PORT;
                InetAddress address;
                
                //make sure the line isn't blank
                if (lineSplit != null) {
                    try {
                        //find the a vnode's IP or valid the IP
                        address = InetAddress.getByName(lineSplit[0]);
                        String name = address.getHostAddress();
                        
                        host_names.add(name);
                        
                        //make sure that the port was specified
                        if (lineSplit.length == 2) {
                            //add the port number to the default port
                            hostPort += new Integer(lineSplit[1]).intValue();
                        }
                        
                        host_ports.add(""+hostPort);
                        
                        
                        //sees if this vnode's IP has been found
                        if (activeInterface == null) {
                            try {
                                activeInterface = NetworkInterface.getByInetAddress(address);
                            }
                            //exception is useless information
                            catch (Exception e) {
                                //empty
                            }
                            
                            //sets the hostName, host, and port if this IP is for this vnode
                            if (activeInterface != null) {
                                hostName = name;
                                host=i;
                                port=hostPort;
                            }
                        }
                        numberOfCPUs++;
                        i++;
                    }
                    catch (java.net.UnknownHostException e) {
                        throw new CCJException("Could not find host name " + lineSplit[0] + " " + e);
                    }
                }
            }
            
            try {
                line = in.readLine();
            }
            catch (Exception e) {
                throw new CCJException("An error occurred reading from the CcjHosts file: " + e);
                
            }
        }
        
        //close the ReadBuffer/file
        try {
            in.close();
        }
        catch (Exception e) {
            throw new CCJException("An error occurred while trying to close the CcjHosts file: " + e);
        }
    }
    
    /**
     * accessor for hostName
     *
     * @return this vnode's IP
     */
    public String getMyHost() {
        return hostName;
    }
    
    /**
     * accessor for port
     *
     * @return this vnode's port
     */
    public int getMyPort() {
        return port;
    }
    
    /**
     * accessor for host
     *
     * @return this vnode's rank
     */
    public int getMyHostNum() {
        return host;
    }
    
    /**
     * accessor for numberOfCPUs
     *
     * @return the number of vnode's in the group
     */
    public int getNumCpus() {
        return numberOfCPUs;
    }
    
    /**
     * Used by ColGroupMaster when trying to identify which port to bind to.
     * Steps through each IP and port combination, and finds the current IP and port
     * then finds the next IP and port combination, and trys to bind to that port.
     *
     * @return port to try to bind to or 0 if is there is no new port to try to bind to
     */
    public int findNextHostPort() {
        boolean isNext = false;
        for(int i=0; i < host_names.size(); i++) {
            //sets isNext to true if the next port found to match this vnode's ip address should be used to bind to
            if (hostName.equals((String)host_names.get(i)) && new Integer((String) host_ports.get(i)).intValue() == port) {
                isNext=true;
            }
            //this port is the port to try to bind to
            else if (hostName.equals((String)host_names.get(i)) && isNext == true) {
                String sPort = (String) host_ports.get(i);
                port = new Integer(sPort).intValue();
                return port;
            }
        }
        return 0;
    }
}