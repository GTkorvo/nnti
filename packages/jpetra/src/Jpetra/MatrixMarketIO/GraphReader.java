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

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.net.URL;

import java.util.Set;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import Jpetra.*;

public class GraphReader extends JpetraObject {
    
    public GraphReader() {
    }
    
        public static Graph read(String fileName, boolean rowOriented, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            return null;
        }
        
        FileInputStream fis = new FileInputStream(fileName);
        return doRead(fileName, fis, rowOriented, comm);
    }
    
    public static Graph readUrl(String urlString, boolean rowOriented, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            return null;
        }
        
        URL url = new URL(urlString);
        return doRead(urlString, url.openStream(), rowOriented, comm);
    }
    
    // need to make it generate its own vector space
    private static Graph doRead(String fileName, InputStream is, boolean rowOriented, Comm comm) throws java.io.IOException {
        // Open file for reading. If it's compressed, use on the fly
        // decompression
        InputStreamReader isr = null;
        if (fileName.endsWith("gz"))
            isr = new InputStreamReader(new GZIPInputStream(is));
        else
            isr = new InputStreamReader(is);
        MatrixVectorReader mvr = new MatrixVectorReader(isr);
        
        // Read header
        MatrixInfo info = mvr.readMatrixInfo();
        String[] comments = mvr.readComments();
        MatrixSize size = mvr.readMatrixSize(info);
        
        
        Graph result = null;
        
        int[] row = null, col = null;
        Object data = null, dataR = null, dataI = null;
        
        // Call appropriate parser
        if (info.isDense()) {
            if (info.isInteger()) {
                //data = new int[size.numEntries()];
                //mvr.readArray((int[]) data);
                // !! need to change type ?
            } else if (info.isReal()) {
                //data = new double[size.numEntries()];
                //mvr.readArray((double[]) data);
            } else if (info.isComplex()) {
                //dataR = new double[size.numEntries()];
                //dataI = new double[size.numEntries()];
                //mvr.readArray((double[]) dataR, (double[]) dataI);
                JpetraObject.println("ERR", "Complex data type not supported by Graph.");
            } else
                throw new IOException("Parser error");
        } else {
            row = new int[size.numEntries()];
            col = new int[size.numEntries()];
            if (info.isInteger()) {
                //data = new int[size.numEntries()];
                //mvr.readCoordinate(row, col, (int[]) data);
                // !! need to do something here
            } else if (info.isReal()) {
                //data = new double[size.numEntries()];
                //mvr.readCoordinate(row, col, (double[]) data);
            } else if (info.isComplex()) {
                /*dataR = new double[size.numEntries()];
                dataI = new double[size.numEntries()];
                mvr.readCoordinate(
                row,
                col,
                (double[]) dataR,
                (double[]) dataI);
                 shouldn't end up here.*/
            } else if (info.isPattern()) {
                mvr.readPattern(row, col);
                ElementSpace myElementSpace;
                if (rowOriented) {
                    myElementSpace = new ElementSpace(size.numRows(), 0, comm);
                }
                else {
                    myElementSpace = new ElementSpace(size.numColumns(), 0, comm);
                }
                VectorSpace myVectorSpace = new VectorSpace(myElementSpace);
                result = buildGraphFromSparse(row, col, myVectorSpace, rowOriented);
            }
            // !! need to add some support for this one!
            else
                throw new IOException("Parser error");
        }
        
        mvr.close();
        is.close();
        
        return result;
    }
    
    public static Graph buildGraphFromSparse(int[] row, int[] col, VectorSpace vectorSpace, boolean rowOriented) {
        TreeMap OuterTree = new TreeMap();
        int[] rowColMajor;
        int[] rowColMinor;
        if (rowOriented) {
            rowColMajor = row;
            rowColMinor = col;
        }
        else {
            rowColMinor = row;
            rowColMajor = col;
        }
        
        int numElements = 0;
        int numTotalEntries = 0;
        // need to see if the TreeMap exists for the specified row/col
        for(int i=0; i < rowColMajor.length; i++) {
            TreeMap rowColTreeMap;
            if (!OuterTree.containsKey(new Integer(rowColMajor[i]))) {
                rowColTreeMap = new TreeMap();
                //this.println("STD", "globalRowCol does not exist, creating...");
                OuterTree.put(new Integer(rowColMajor[i]), rowColTreeMap);
                numElements++;
            }
            else {
                //this.println("STD", "globalRowCol exists, setting rowColTreeMap to existing TreeMap...");
                rowColTreeMap = (TreeMap) OuterTree.get(new Integer(rowColMajor[i]));
            }
            
            // now that we know the row/col exists, insert entries into the row/col
            // and sum them with any pre-existing entries
            //this.println("STD", "CisMatrix is internally inserting an entry...");
            Object temp;
            int toInsert;
            //this.println("STD", "Inserting index " + index);
            temp = rowColTreeMap.get(new Integer(rowColMinor[i]));
            rowColTreeMap.put(new Integer(rowColMinor[i]), new Integer(1));
            numTotalEntries++;
            /*if (temp == null) {
                //this.println("STD", "Index " + index + " does not exist, creating...");
                rowColTreeMap.put(new Integer(rowColMinor[i]), new Integer(1));
            }
            else {
                //this.println("STD", "Index " + index + " does exists, adding to prevous value...");
                //rowColTreeMap.put(new Integer(index), new Double(entry + ((Double) temp).doubleValue()));
            }*/
        }
        
        //this.doubleValues = new double[this.numTotalEntries];
        int[] numEntries = new int[numElements];
        
        Set outerKeysValues = OuterTree.entrySet();
        Iterator outerIterator = outerKeysValues.iterator();
        Map.Entry outerMapEntry;
        
        TreeMap innerTree;
        Set innerKeysValues;
        Iterator innerIterator;
        Map.Entry innerMapEntry;
        
        int nextEntryIndex=0;
        int numEntriesColRow;
        int[] tempGraph = new int[numTotalEntries];
        int i=0;
        //startIndex = 0;
        while(outerIterator.hasNext()) {
            //this.println("STD", "Doing an outer loop...");
            outerMapEntry = (Map.Entry) outerIterator.next();
            innerTree = (TreeMap) outerMapEntry.getValue();
            innerKeysValues = innerTree.entrySet();
            innerIterator = innerKeysValues.iterator();
            numEntriesColRow = 0;
            while (innerIterator.hasNext()) {
                //this.println("STD", "Doing an inner loop...");
                innerMapEntry = (Map.Entry) innerIterator.next();
                //this.doubleValues[nextEntryIndex] = ((Double) innerMapEntry.getValue()).doubleValue();
                tempGraph[nextEntryIndex++] = ((Integer) innerMapEntry.getKey()).intValue() - 1;
                numEntriesColRow++;
            }
            //if (this.maxSecondaryId < tempGraph[nextEntryIndex-1]) {
            //    this.maxSecondaryId = tempGraph[nextEntryIndex-1];
            //}
            //this.startIndex[i] = startIndex;
            //startIndex += numEntriesColRow;
            numEntries[i++] = numEntriesColRow;
        }
        
        return new Graph(vectorSpace, tempGraph, numEntries);
    }
}
