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

package Jpetra.MatrixMarketIO;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.net.Socket;

import Jpetra.*;

/**
 * Writes a <code>CisMatrix</code> to a plain text MatrixMarket formatted
 * file or to a TCP/IP socket for display in Matvis.
 *
 * @author Jason Cross
 */
public class CisMatrixWriter extends JpetraObject {
    /**
     * Opens a TCP/IP connection to Matvis and sends the <code>CisMatrix</code> to Matvis for dispaly.
     *
     * @param server The server to send to.  The server can be an IP address or a name.  Ex: 127.0.0.1 or localhost
     * @param port The TCP/IP port on which Matvis is listening for connections.
     * @param cisMatrix The <code>CisMatri</code> to send to Matvis.
     *
     * @throws IOException Any errors that occur while trying to connect to or send the <code>CisMatrix></code> to Matvis are thrown.
     */
    public static void sendToMatvis(String server, int port, CisMatrix cisMatrix, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() == 0) {
            Socket socket = new Socket(server, port);
            write(new PrintWriter(socket.getOutputStream()), cisMatrix);
            socket.close();
        }
    }
    
    /**
     * Writes a <code>CisMatrix</code> to a plan text file in MatrixMarket format.
     *
     * @param fileName the absolute or relative path of the file to write to
     * @param cisMatrix the <code>CisMatrix<code> to write to the file
     *
     * @throws IOException Any errors that occur while trying to open or read form the file are thrown.
     */
    public static void write(String fileName, CisMatrix cisMatrix, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() == 0) {
            write(new PrintWriter(new FileOutputStream(fileName)), cisMatrix);
        }
    }
    
    /**
     * Writes a <code>CisMatrix</code> in MatrixMarket format in plain text to the specified <code>PrintWriter</code>.
     *
     * @param out the <code>PrintStream</code> that the <code>CisMatrix</code> is printed to
     * @param cisMatrix the <code>CisMatrix</code> to print out
     *
     * @throws IOException Any IO errors that occur while tring to print to the <code>PrintStream</code> are thrown.
     */
    private static void write(PrintWriter out, CisMatrix cisMatrix) throws java.io.IOException {
        out.println("%%MatrixMarket matrix coordinate real general");
        out.println(cisMatrix.getNumMyRows() + " " + cisMatrix.getNumMyColumns() + " " + cisMatrix.getNumMyNonZeros());
        out.println();
        int[] numEntries = cisMatrix.getNumEntriesArray();
        double[] entries = cisMatrix.getEntriesArray();
        int[] nonZeros = cisMatrix.getGraph().getNonZeroEntriesArray();
        
        if (cisMatrix.isRowOriented()) {
            JpetraObject.println("STD", "Matrix is Row Oriented.");
            int row = 0;
            int index = 0;
            for(int i=0; i < numEntries.length; i++) {
                for(int j=0; j < numEntries[i]; j++) {
                    out.println((i + 1) + " " + (nonZeros[index] + 1) + " " + entries[index++]);
                }
            }
        }
        else {
            JpetraObject.println("STD", "Matrix is Col Oriented.");
            int col = 0;
            int index = 0;
            for(int i=0; i < numEntries.length; i++) {
                for(int j=0; j < numEntries[i]; j++) {
                    out.println((nonZeros[index] + 1) + " " + (i + 1) + " " + entries[index++]);
                }
            }
            
        }
        out.close();
    }
}
