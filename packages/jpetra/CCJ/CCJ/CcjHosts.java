package CCJ;

import java.util.Vector;
import java.io.*;
import java.net.*;
import java.rmi.registry.Registry;
/*
 * CcjHosts.java
 *
 * Created on Mon June 23 22:05:00 CDT 2003
 */

/**
 * Reads a file from the file name path passed to it, and then reads in a host per line in the format:<br>
 * <code>host:port</code><br>
 * The port number will be added to the default port value specified by <code>import java.rmi.registry.Registry;</code>.
 *
 * To deal with the fact that InetAddress.getLocalHost() has some very odd behavior
 * and the fact that machines may have more than one valid IP address
 * this modification uses NetworkInterface.getByInetAddress to determine if
 * the IP address specified is bound to any local network interface.
 * Effectivly this means that any local interface or IP could be used.
 *
 * NOTE/BUGS: Additional error checking should be included to eliminate the ability
 * to specify loopback interfaces as valid IPs.  Also, the first IP found to be bound
 * to an interface will be used.  This will prevent multiple processes on the same
 * node from using multiple IP addresses.
 *
 * @author  Jason Cross 
 */

public class CcjHosts {
    private String[] hosts;
    private String[] ports;
    private String hostName;
    private int port;
    private int host;
    private int numberOfCPUs;
    private Vector host_names;
    private Vector host_ports;
    
    public CcjHosts (String filePath, Vector host_names, Vector host_ports) throws CCJException {
        this.host_names = host_names;
        this.host_ports = host_ports;
                
        BufferedReader in = null;
        String line = null;
        try {
            in = new BufferedReader(new FileReader(filePath));
            
            line = in.readLine();
        }
        catch (Exception e) {
            throw new CCJException("An error occurred while reading from the CcjHosts file: " + e);
        }
       
        int i = 0;
        NetworkInterface activeInterface = null;
        host=-1;
        numberOfCPUs = 0;
        while (line != null) {
            String[] lineSplit = line.split(":");
            int hostPort=Registry.REGISTRY_PORT;
            InetAddress address;
            
            if (lineSplit != null) {
                try {
                    address = InetAddress.getByName(lineSplit[0]);
                    String name = address.getHostAddress();
            
                    host_names.add(name);                
                
                if (lineSplit.length == 2) {
                    hostPort += new Integer(lineSplit[1]).intValue();
                }

                host_ports.add(""+hostPort);
                
                if (activeInterface == null) {
                    try {
                        activeInterface = NetworkInterface.getByInetAddress(address);
                    }
                    //exception is useless information
                    catch (Exception e) {}
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

            try {
                line = in.readLine();
            }
            catch (Exception e) {
                throw new CCJException("An error occurred reading from the CcjHosts file: " + e);
                
            }
        }
        try {
            in.close();
        }
        catch (Exception e) {
            throw new CCJException("An error occurred while trying to close the CcjHosts file: " + e);
        }    
    }
    
    public String getMyHost () {
        return hostName;
    }
    
    public int getMyPort () {
        return port;
    }
    
    public int getMyHostNum () {
        return host;
    }
    
    public int getNumCpus () {
        return numberOfCPUs;
    }
    
    public int findNextHostPort () {
        boolean isNext = false;
        for(int i=0; i < host_names.size(); i++) {
            if (hostName.equals((String)host_names.get(i)) && new Integer((String) host_ports.get(i)).intValue() == port) {
                isNext=true;
            }
            else if (hostName.equals((String)host_names.get(i)) && isNext == true) {
                String sPort = (String) host_ports.get(i);
                port = new Integer(sPort).intValue();
                return port;
            }
        }
        return 0;
    }
}