package JpetraVis;

/*
 * MatView.java
 *
 * Created on February 28, 2002
 * Last modified on March 19, 2002
 *
 */

/**
 * @author  kahansen
 * @version 
 */

import javax.swing.*;
import Jpetra.*;

public class MatView {

    private String matrix;
    
    /** Creates new MatView 
    
	 @param matrixMarket the matrix data in Matrix Market format 
    */
    public MatView(String matrixMarket) {
    	matrix = matrixMarket;
    }

    /**  Displays the MatView object graphically to the screen  */
   public void display() {
       int port = 1099;
     //new Thread(new matvis.SocketListener(null, null,
                  //new matvis.matvisMain(), port)).start();
     //System.out.println("main, launched the SocketListener thread...");

     java.net.InetAddress inetadd = getLocalInetAddress();
     if (inetadd == null) System.exit(1);
     System.out.println("Inet address = "+inetadd+" Port = "+port);

     java.net.Socket socket = null;
     try {
         socket = new java.net.Socket(inetadd, port);
     }
     catch(Exception e) {
         System.out.println("main, failed to create socket...");
         System.exit(1);
     }

     System.out.println("main, got a client socket. now sending data...");

     java.io.OutputStream outStream = null;
     try {
         outStream = socket.getOutputStream();
     }
     catch(Exception e2) {
         System.out.println("main, exception getting socket's outStream.");
         System.exit(1);
     }

     java.io.OutputStreamWriter strmWriter =
         new java.io.OutputStreamWriter(outStream);

     try {
         strmWriter.write(matrix, 0, matrix.length());
     }
     catch(Exception e3) {
         System.out.println("Exception writing to socket's stream.");
         System.exit(1);
     }

     try {
         strmWriter.close();
     }
     catch(Exception e4) {
         System.out.println("Exception closing socket's stream.");
     }

     try {
         socket.close();
     }
     catch(Exception e4) {
         System.out.println("Exception closing socket.");
     }

     System.out.println("main, closed socket, all done.");
     int n  = JOptionPane.showConfirmDialog(
                      null, "Did a new MatrixFrame appear?",
                            "Confirm matrix socket test",
                            JOptionPane.YES_NO_OPTION);
     if (n != JOptionPane.YES_OPTION) {
         System.exit(1);
     }
   }
   
   public java.net.InetAddress getLocalInetAddress() {
        java.net.InetAddress iadd = null;
        try {
            iadd = java.net.InetAddress.getLocalHost();
        }
        catch(java.net.UnknownHostException uhe) {
            System.out.println("!! unknown-host-exception !!");
         return(null);
        }

     return(iadd);
   }

}
