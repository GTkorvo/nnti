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

import Jpetra.*;

public class MultiVectorReader extends JpetraObject {
    
    public MultiVectorReader() {
        // empty
    }
    
    public static MultiVector read(String fileName, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            VectorSpace vectorSpace = new VectorSpace(new ElementSpace(-1, 0, 0, comm));
            MultiVector multiVector = new MultiVector(vectorSpace);
            return multiVector;
        }
        
        FileInputStream fis = new FileInputStream(fileName);
        return doRead(fileName, fis, comm);
    }
    
    public static MultiVector readUrl(String urlString, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            VectorSpace vectorSpace = new VectorSpace(new ElementSpace(-1, 0, 0, comm));
            MultiVector multiVector = new MultiVector(vectorSpace);
            return multiVector;
        }
        
        URL url = new URL(urlString);
        return doRead(urlString, url.openStream(), comm);
    }
    
    private static MultiVector doRead(String fileName, InputStream is, Comm comm) throws java.io.IOException {
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
        
        ElementSpace myElementSpace = new ElementSpace(-1, size.numRows(), 0, comm);
        VectorSpace myVectorSpace = new VectorSpace(myElementSpace);
        MultiVector result = null;
        
        int[] row = null, col = null;
        Object data = null, dataR = null, dataI = null;
        
        // Call appropriate parser
        if (info.isDense()) {
            if (info.isInteger()) {
                data = new int[size.numEntries()];
                mvr.readArray((int[]) data);
                // strip out zeros and add values to the CisMatrix
                result = buildMultiVectorFromIntDense(size, (int[]) data, myVectorSpace);
            } else if (info.isReal()) {
                data = new double[size.numEntries()];
                mvr.readArray((double[]) data);
                // strip out zeros and add values to the CisMatrix
                result = buildMultiVectorFromDense(size, (double[]) data, myVectorSpace);
            } else if (info.isComplex()) {
                JpetraObject.println("FATALERR", "In MultiVectorReader: Complex data type not supported by MutliVector.");
                System.exit(1);
            } else
                throw new IOException("MultiVectorReader does not support the vector format of the dense vector being read in.  Supported types are integer, pattern, and real.");
        } else {
            row = new int[size.numEntries()];
            col = new int[size.numEntries()];
            if (info.isInteger()) {
                data = new int[size.numEntries()];
                mvr.readCoordinate(row, col, (int[]) data);
                result = buildMultiVectorFromIntSparse(row, col, size, (int[]) data, myVectorSpace);
            } else if (info.isReal()) {
                data = new double[size.numEntries()];
                mvr.readCoordinate(row, col, (double[]) data);
                result = buildMultiVectorFromSparse(row, col, size, (double[]) data, myVectorSpace);
            } else if (info.isComplex()) {
                JpetraObject.println("FATALERR", "In MultiVectorReader: Complex data type not supported by MutliVector.");
                System.exit(1);
            } else if (info.isPattern()) {
                mvr.readPattern(row, col);
                result = buildMultiVectorFromPatternSparse(row, col, size, myVectorSpace);
            }
            else
                throw new IOException("MultiVectorReader does not support the vector format of the sparse vector being read in.  Supported types are integer, pattern, and real.");
        }
        
        mvr.close();
        is.close();
        
        return result;
    }
    
    private static MultiVector buildMultiVectorFromPatternSparse(int[] row, int[] col, MatrixSize size, VectorSpace myVectorSpace) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        double[][] values = new double[numCols][numRows];
        for(int i=0; i < size.numEntries(); i++) {
            values[col[i]-1][row[i]-1] = 1;
        }
        
        return new MultiVector(myVectorSpace, values);
    }
    
    private static MultiVector buildMultiVectorFromSparse(int[] row, int[] col, MatrixSize size, double[] data, VectorSpace myVectorSpace) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        double[][] values = new double[numCols][numRows];
        for(int i=0; i < size.numEntries(); i++) {
            values[col[i]-1][row[i]-1] = data[i];
        }
        
        return new MultiVector(myVectorSpace, values);
    }
    
    private static MultiVector buildMultiVectorFromIntSparse(int[] row, int[] col, MatrixSize size, int[] data, VectorSpace myVectorSpace) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        double[][] values = new double[numCols][numRows];
        for(int i=0; i < size.numEntries(); i++) {
            values[col[i]-1][row[i]-1] = data[i];
        }
        
        return new MultiVector(myVectorSpace, values);
    }
    
    private static MultiVector buildMultiVectorFromDense(MatrixSize size, double[] data, VectorSpace myVectorSpace) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        int entryIndex = 0;
        double[][] values = new double[numCols][numRows];
        for(int col=0; col < numCols; col++) {
            for(int row=0; row < numRows; row++) {
                values[col][row] = data[entryIndex++];
            }
        }
        
        return new MultiVector(myVectorSpace, values);
    }
    
    private static MultiVector buildMultiVectorFromIntDense(MatrixSize size, int[] data, VectorSpace myVectorSpace) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        int entryIndex = 0;
        double[][] values = new double[numCols][numRows];
        for(int col=0; col < numCols; col++) {
            for(int row=0; row < numRows; row++) {
                values[col][row] = data[entryIndex++];
            }
        }
        
        return new MultiVector(myVectorSpace, values);
    }
}
