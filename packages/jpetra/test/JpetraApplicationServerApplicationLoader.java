package test;

import Jpetra.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.FileReader;
import java.net.InetAddress;
import java.net.Socket;
import java.io.OutputStream;
import java.io.ObjectOutputStream;

/**
 *
 * @author  Jason Cross
 */
public class JpetraApplicationServerApplicationLoader extends JpetraObject {
    
    /** Creates a new instance of JpetraApplicationServer */
    public JpetraApplicationServerApplicationLoader(String[] args) {
        String className = args[2];
        String[] appArgs;
        if (args.length > 3) {
            appArgs = new String[args.length-3];
            System.arraycopy(args, 3, appArgs, 0, appArgs.length);
        }
        else {
            appArgs = new String[0];
        }
        try {
            File inFile = new File(args[1] + "/" + className + ".class");
            long inFileLength = inFile.length();
        /*
        int num2dBytes = (int) (inFileLength/Integer.MAX_VALUE);
        num2dBytes++;
        long bytesLeft = inFileLength;
        byte[][] rawFileData = new byte[num2dBytes][];
        FileInputStream fis = new FileInputStream(inFile);
        for(int i=0; i < num2dBytes; i++) {
            if ((bytesLeft - Integer.MAX_VALUE) > 0) {
                rawFileData[i] = new byte[Integer.MAX_VALUE];
                bytesLeft -= Integer.MAX_VALUE;
            } else {
                if (bytesLeft == 0) {
                    System.err.println("Fatal Error: bytesLeft > 0 and rawFileData is not filled yet.");
                    System.exit(1);
                }
                rawFileData[i] = new byte[(int) bytesLeft];
                bytesLeft = 0;
            }
            fis.read(rawFileData[i]);
        }*/
            
            if (inFileLength > Integer.MAX_VALUE) {
                System.err.println("FATAL ERROR: The class file is too big to be sent.");
                System.exit(1);
            }
            FileInputStream fis = new FileInputStream(inFile);
            byte[] rawInFileData = new byte[(int) inFileLength];
            fis.read(rawInFileData);
            fis.close();
            BufferedReader nodesFile = new BufferedReader(new FileReader(args[0]));
            String line = nodesFile.readLine();
            int numNodes = 0;
            byte[] temp = new byte[4];
            if (line != null) {
                numNodes = Integer.parseInt(line);
            }
            InetAddress[] nodes = new InetAddress[numNodes];
            int[] ports = new int[numNodes];
            for(int i=0; i < nodes.length; i++) {
                temp[0] = Byte.parseByte(nodesFile.readLine());
                temp[1] = Byte.parseByte(nodesFile.readLine());
                temp[2] = Byte.parseByte(nodesFile.readLine());
                temp[3] = Byte.parseByte(nodesFile.readLine());
                ports[i] = Integer.parseInt(nodesFile.readLine());
                nodes[i] = InetAddress.getByAddress(temp);
            }
            nodesFile.close();
            
            for(int i=0; i < nodes.length; i++) {
                Socket connection = new Socket(nodes[i], ports[i]);
                OutputStream os = connection.getOutputStream();
                ObjectOutputStream oos = new ObjectOutputStream(os);
                oos.writeObject(className);
                oos.writeObject(appArgs);
                oos.writeObject(rawInFileData);
                
                oos.close();
                connection.close();
            }
        }
        catch (java.io.FileNotFoundException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
        catch (java.io.IOException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length < 3) {
            System.out.println("Usage: java test/JpetraApplicationServerClassLoader [ipsAndPorts.txt] [pathToClass] [JpetraServerApplicationClassName] (application args)");
            System.exit(0);
        }
        
        new JpetraApplicationServerApplicationLoader(args);
    }
    
}