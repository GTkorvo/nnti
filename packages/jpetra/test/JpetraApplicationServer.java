package test;

import Jpetra.*;

import java.net.ServerSocket;
import java.net.Socket;
import java.net.InetAddress;
import java.io.PrintWriter;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.security.SecureClassLoader;
/**
 *
 * @author  Jason Cross
 */
public class JpetraApplicationServer extends JpetraObject {
    Comm comm;
    ServerSocket serverSocket;
    
    /** Creates a new instance of JpetraApplicationServer */
    public JpetraApplicationServer(String[] args) {
        initializeOutput();
        this.comm = new CcjComm(args[0]);
        
        this.createServerSocket(args[1]);
        try {
            JpetraClassLoader classLoader = new JpetraClassLoader();
            while(true) {
                Socket connection = this.serverSocket.accept();
                ObjectInputStream ois = new ObjectInputStream(connection.getInputStream());
                String classToLoad = (String) ois.readObject();
                String[] appArgs = (String[]) ois.readObject();
                classLoader.setObjectInputStream(ois);
                Class applicationClass = classLoader.findClass(classToLoad);
                ois.close();
                this.runApplication((JpetraServerApplication) applicationClass.newInstance(), appArgs);
                applicationClass = null;
                classToLoad = null;
                connection.close();
                connection = null;
                this.cleanupAfterApplicationRun();
            }
        }
        catch (java.io.IOException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
        catch (java.lang.ClassNotFoundException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
        catch (java.lang.InstantiationException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
        catch (java.lang.IllegalAccessException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
    }
    
    public void runApplication(JpetraServerApplication application, String[] appArgs) {
        application.setComm(comm);
        application.runApplication(appArgs);
    }
    
    public void cleanupAfterApplicationRun() {
        System.gc();
    }
    
    public void createServerSocket(String fileName) {
        try {
            this.serverSocket = new ServerSocket(0);
        }
        catch (java.io.IOException e) {
            System.err.println(e.toString());
        }
        int[] myIpAndPort = new int[5];
        //byte[] myRawIp = this.serverSocket.getInetAddress().getAddress();
        byte[] myRawIp = null;
        try {
            myRawIp = InetAddress.getLocalHost().getAddress();
            System.out.println(InetAddress.getLocalHost().getHostAddress());
        }
        catch (java.net.UnknownHostException e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        
        System.out.println(myRawIp.length);
        myIpAndPort[0] = myRawIp[0];
        System.out.println(myRawIp[0]);
        myIpAndPort[1] = myRawIp[1];
        System.out.println(myRawIp[1]);
        myIpAndPort[2] = myRawIp[2];
        System.out.println(myRawIp[2]);
        myIpAndPort[3] = myRawIp[3];
        System.out.println(myRawIp[3]);
        myIpAndPort[4] = this.serverSocket.getLocalPort();
        int[][] nodes = this.comm.gather(myIpAndPort);
        
        if (this.comm.getVnodeId() == 0) {
            try {
                PrintWriter outFile = new PrintWriter(new FileOutputStream(fileName));
                outFile.println(this.comm.getNumVnodes());
                for(int i=0; i < nodes.length; i++) {
                    for(int j=0; j < nodes[i].length; j++) {
                        System.out.println(nodes[i][j] + "");
                        outFile.println(nodes[i][j] + "");
                    }
                }
                outFile.close();
            }
            catch (java.io.IOException e) {
                System.err.println(e.toString());
            }
        }
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length < 2) {
            System.out.println("Usage: java test/JpetraApplicationServer [ccjhosts.txt] [nodeIpsAndPorts.txt]");
            System.exit(0);
        }
        
        new JpetraApplicationServer(args);
    }
    
}
