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

import Jpetra.*;

/**
 *
 * @author  Jason Cross
 */
public class GraphWriter extends JpetraObject {
    
    public GraphWriter() {
    }
    
    public static void write(String fileName, Graph graph, CisMatrix cisMatrix, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            return;
        }
        
        PrintWriter out = new PrintWriter(new FileOutputStream(fileName));
        
        out.println("%%MatrixMarket matrix coordinate pattern general");
        out.println(cisMatrix.getNumMyRows() + " " + cisMatrix.getNumMyColumns() + " " + graph.getNumNonZeros());
        out.println();
        int[] numEntries = graph.getNumEntriesArray();
        int[] nonZeros = graph.getNonZeroEntriesArray();
        
        if (cisMatrix.isRowOriented()) {
            JpetraObject.println("STD", "Graph is Row Oriented.");
            int row = 0;
            int index = 0;
            for(int i=0; i < numEntries.length; i++) {
                for(int j=0; j < numEntries[i]; j++) {
                    out.println((i + 1) + " " + (nonZeros[index++] + 1));
                }
            }
        }
        else {
            JpetraObject.println("STD", "Graph is Col Oriented.");
            int col = 0;
            int index = 0;
            for(int i=0; i < numEntries.length; i++) {
                for(int j=0; j < numEntries[i]; j++) {
                    out.println((nonZeros[index++] + 1) + " " + (i + 1));
                }
            }
            
        }
        out.close();
    }
}
