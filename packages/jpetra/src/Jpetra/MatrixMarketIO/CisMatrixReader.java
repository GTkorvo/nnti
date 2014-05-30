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

/**
 * <code>CisMatrixReader</code> can read a MatrixMarket formatted dense or sparse matrix
 * from a plain text or gzipped (.gz) file or url.  The url protocol has to be supported
 * by Java such as HTTP and FTP.  In order for the gzip decompression to work properly
 * the file or url name must end in gz.
 */
public class CisMatrixReader extends JpetraObject {
    /**
     * Reads a MatrixMarket formatted plain text or gzipped file from a url.
     *
     * @param urlString The url in for the format protocol://location.  Ex: http://test.com/matrix.gz
     * @param rowOriented <code>true</code> if the matrix should be read into row major form
     * @param comm the <code>comm</code> object to be associated with the new CisMatrix
     *
     * @return A CisMatrix that has had fillComplete called on it.
     * @throws IOException Any IO errors that occur while reading from the url will be thrown.
     */
    public static CisMatrix readUrl(String urlString, boolean rowOriented, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            VectorSpace vectorSpace = new VectorSpace(new ElementSpace(-1, 0, 0, comm));
            CisMatrix cisMatrix = new CisMatrix(vectorSpace, rowOriented);
            cisMatrix.fillComplete();
            return cisMatrix;
        }
        
        URL url = new URL(urlString);
        return doRead(urlString, url.openStream(), rowOriented, comm);
    }
    
    /**
     * Reads from a MatrixMarket formatted plain text or gzipped file.
     *
     * @param fileName The relative or absolute file name of the file to be read in.
     * @param rowOriented <code>true</code> if the matrix should be read into row major form
     * @param comm the <code>comm</code> object to be associated with the new CisMatrix
     *
     * @return A CisMatrix that has had fillComplete called on it.
     * @throws IOException Any IO errors that occur while reading from the file will be thrown.
     */
    public static CisMatrix read(String fileName, boolean rowOriented, Comm comm) throws java.io.IOException {
        if (comm.getVnodeId() != 0) {
            VectorSpace vectorSpace = new VectorSpace(new ElementSpace(-1, 0, 0, comm));
            CisMatrix cisMatrix = new CisMatrix(vectorSpace, rowOriented);
            cisMatrix.fillComplete();
            return cisMatrix;
        }
        
        FileInputStream fis = new FileInputStream(fileName);
        return doRead(fileName, fis, rowOriented, comm);
    }
    
    /**
     * Reads a MatrixMarket formatted from an InputStream.
     *
     * @param fileName The url string or file name of the file or url being read.  If it ends in gz then the <code>InputStream is</code> treated as a compressed gzip stream.
     * @param is The <code>InputStream</code> to read from.
     * @param rowOriented <code>true</code> if the matrix should be read into row major form
     * @param comm the <code>comm</code> object to be associated with the new CisMatrix
     *
     * @return A CisMatrix that has had fillComplete called on it.
     * @throws IOException Any IO errors that occur while reading from the <code>InputStream</code> will be thrown.
     */
    private static CisMatrix doRead(String fileName, InputStream is, boolean rowOriented, Comm comm) throws java.io.IOException {
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
        
        ElementSpace myElementSpace;
        if (rowOriented) {
            myElementSpace = new ElementSpace(-1, size.numRows(), 0, comm);
        }
        else {
            myElementSpace = new ElementSpace(-1, size.numColumns(), 0, comm);
        }
        VectorSpace myVectorSpace = new VectorSpace(myElementSpace);
        CisMatrix result = new CisMatrix(myVectorSpace, rowOriented);
        
        int[] row = null, col = null;
        Object data = null, dataR = null, dataI = null;
        
        // Call appropriate parser
        if (info.isDense()) {
            if (info.isInteger()) {
                data = new int[size.numEntries()];
                mvr.readArray((int[]) data);
                // strip out zeros and add values to the CisMatrix
                buildCisMatrixFromIntDense((int[]) data, size, result);
            } else if (info.isReal()) {
                data = new double[size.numEntries()];
                mvr.readArray((double[]) data);
                // strip out zeros and add values to the CisMatrix
                buildCisMatrixFromDense((double[]) data, size, result);
            } else if (info.isComplex()) {
                JpetraObject.println("FATALERR", "In CisMatrixReader: Complex data type not supported by CisMatrix.");
                System.exit(1);
            } else
                throw new IOException("CisMatrixReader does not understand the type of dense matrix being read.  Supported types are interger and real.");
        } else {
            row = new int[size.numEntries()];
            col = new int[size.numEntries()];
            if (info.isInteger()) {
                data = new int[size.numEntries()];
                mvr.readCoordinate(row, col, (int[]) data);
                buildCisMatrixFromIntSparse(row, col, (int[]) data, result);
            } else if (info.isReal()) {
                data = new double[size.numEntries()];
                mvr.readCoordinate(row, col, (double[]) data);
                buildCisMatrixFromSparse(row, col, (double[]) data, result);
            } else if (info.isComplex()) {
                JpetraObject.println("FATALERR", "In CisMatrixReader: Complex data type not supported by CisMatrix.");
                System.exit(1);
            } else if (info.isPattern()) {
                mvr.readPattern(row, col);
                buildCisMatrixFromPatternSparse(row, col, result);
            }
            else
                throw new IOException("CisMatrixReader does not understand the type of sparse matrix being read.  Supported types are interger and real.");
        }
        
        mvr.close();
        is.close();
        
        return result;
    }
    
    /**
     * Inserts values from a double array into the <code>CisMatrix</code>.
     *
     * @param row for int i, row[i] is the row index for entry data[i]
     * @param col for int i, col[i] is the row index for entry data[i]
     * @param data the entry values of the matrix that was read in
     * @param cisMatrix the <code>CisMatrix</code> to insert values into
     */
    private static void buildCisMatrixFromSparse(int[] row, int[] col, double[] data, CisMatrix cisMatrix) {
        if (cisMatrix.isRowOriented()) {
            for(int i=0; i < row.length; i++) {
                cisMatrix.insertEntry(row[i], col[i]-1, data[i], DistObject.REPLACE);
            }
        }
        else {
            for(int i=0; i < col.length; i++) {
                cisMatrix.insertEntry(col[i], row[i]-1, data[i], DistObject.REPLACE);
            }
        }
        cisMatrix.fillComplete();
    }
    
    /**
     * Inserts 1's into the <code>CisMatrix</code>.
     *
     * @param row for int i, row[i] is the row index for entry data[i]
     * @param col for int i, col[i] is the row index for entry data[i]
     * @param cisMatrix the <code>CisMatrix</code> to insert values into
     */
    private static void buildCisMatrixFromPatternSparse(int[] row, int[] col, CisMatrix cisMatrix) {
        if (cisMatrix.isRowOriented()) {
            for(int i=0; i < row.length; i++) {
                cisMatrix.insertEntry(row[i], col[i]-1, 1, DistObject.REPLACE);
            }
        }
        else {
            for(int i=0; i < col.length; i++) {
                cisMatrix.insertEntry(col[i], row[i]-1, 1, DistObject.REPLACE);
            }
        }
        cisMatrix.fillComplete();
    }
    
    /**
     * Inserts values from an int array into the <code>CisMatrix</code>.
     *
     * @param row for int i, row[i] is the row index for entry data[i]
     * @param col for int i, col[i] is the row index for entry data[i]
     * @param data the entry values of the matrix that was read in
     * @param cisMatrix the <code>CisMatrix</code> to insert values into
     */
    private static void buildCisMatrixFromIntSparse(int[] row, int[] col, int[] data, CisMatrix cisMatrix) {
        if (cisMatrix.isRowOriented()) {
            for(int i=0; i < row.length; i++) {
                cisMatrix.insertEntry(row[i], col[i]-1, data[i], DistObject.REPLACE);
            }
        }
        else {
            for(int i=0; i < col.length; i++) {
                cisMatrix.insertEntry(col[i], row[i]-1, data[i], DistObject.REPLACE);
            }
        }
        cisMatrix.fillComplete();
    }
    
    /**
     * Inserts values from a double array into the <code>CisMatrix</code>.
     *
     * @param data the entry values of the matrix that was read in
     * @param size contains the number of rows and columns associated with the matrix that was read in
     * @param cisMatrix the <code>CisMatrix</code> to insert values into
     */
    private static void buildCisMatrixFromDense(double[] data, MatrixSize size, CisMatrix cisMatrix) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        int entryIndex = 0;
        if (cisMatrix.isRowOriented()) {
            for(int col=0; col < numCols; col++) {
                for (int row=0; row < numRows; row++) {
                    if (data[entryIndex] != 0) {
                        cisMatrix.insertEntry(row, col, data[entryIndex++], DistObject.REPLACE);
                    }
                }
            }
        }
        else {
            for(int col=0; col < numCols; col++) {
                for (int row=0; row < numRows; row++) {
                    if (data[entryIndex] != 0) {
                        cisMatrix.insertEntry(col, row, data[entryIndex++], DistObject.REPLACE);
                    }
                }
            }
        }
        
        cisMatrix.fillComplete();
    }
    
    /**
     * Inserts values from an int array into the <code>CisMatrix</code>.
     *
     * @param data the entry values of the matrix that was read in
     * @param size contains the number of rows and columns associated with the matrix that was read in
     * @param cisMatrix the <code>CisMatrix</code> to insert values into
     */
    private static void buildCisMatrixFromIntDense(int[] data, MatrixSize size, CisMatrix cisMatrix) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        int entryIndex = 0;
        if (cisMatrix.isRowOriented()) {
            for(int col=0; col < numCols; col++) {
                for (int row=0; row < numRows; row++) {
                    if (data[entryIndex] != 0) {
                        cisMatrix.insertEntry(row, col, data[entryIndex++], DistObject.REPLACE);
                    }
                }
            }
        }
        else {
            for(int col=0; col < numCols; col++) {
                for (int row=0; row < numRows; row++) {
                    if (data[entryIndex] != 0) {
                        cisMatrix.insertEntry(col, row, data[entryIndex++], DistObject.REPLACE);
                    }
                }
            }
        }
        
        cisMatrix.fillComplete();
    }
    
    /* maybe change the following code to work for column major form and use it... but not for now
    public static void buildCisMatrixFromDense(double[] data, MatrixSize size, CisMatrix cisMatrix) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        int numNonZeros;
        int colIndex;
        int entryIndex = 0;
        for (int i=0; i < numRows; i++) {
            numNonZeros = 0;
            // first find the number of non-zero entries
            for(int j=0; j < numCols; j++) {
                if (data[j] != 0) {
                    numNonZeros++;
                }
            }
            double[] entries = new double[numNonZeros];
            int[] indicies = new int[numNonZeros];
            colIndex=0;
            // now build a new entries and indicies array with only the non-zero entries
            for(int j=0; j < numCols; j++) {
                if (data[entryIndex] != 0) {
                    entries[colIndex] = data[entryIndex];
                    indicies[colIndex] = j;
                }
                entryIndex++;
            }
            cisMatrix.insertEntries(i, indicies, entries);
        }
     
        cisMatrix.fillComplete();
    }*/
}
