package Jpetra.MatrixMarketIO;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;

import Jpetra.*;

public class CisMatrixReader extends JpetraObject {
    
    public CisMatrixReader() {
    }
    
    public static CisMatrix read(String fileName, boolean rowOriented, Comm comm) throws java.io.IOException {
        // Open file for reading. If it's compressed, use on the fly
        // decompression
        FileInputStream fis = new FileInputStream(fileName);
        InputStreamReader isr = null;
        if (fileName.endsWith("gz"))
            isr = new InputStreamReader(new GZIPInputStream(fis));
        else
            isr = new InputStreamReader(fis);
        MatrixVectorReader mvr = new MatrixVectorReader(isr);
        
        // Read header
        MatrixInfo info = mvr.readMatrixInfo();
        String[] comments = mvr.readComments();
        MatrixSize size = mvr.readMatrixSize(info);
        
        ElementSpace myElementSpace = new ElementSpace(size.numRows(), comm);
        VectorSpace myVectorSpace = new VectorSpace(myElementSpace);
        CisMatrix result = new CisMatrix(myVectorSpace, rowOriented);
        
        int[] row = null, col = null;
        Object data = null, dataR = null, dataI = null;
        
        // Call appropriate parser
        if (info.isDense()) {
            if (info.isInteger()) {
                //data = new int[size.numEntries()];
                //mvr.readArray((int[]) data);
                // !! need to change type ?
            } else if (info.isReal()) {
                data = new double[size.numEntries()];
                mvr.readArray((double[]) data);
                // strip out zeros and add values to the CisMatrix
                buildCisMatrixFromDense((double[]) data, size, result);
            } else if (info.isComplex()) {
                //dataR = new double[size.numEntries()];
                //dataI = new double[size.numEntries()];
                //mvr.readArray((double[]) dataR, (double[]) dataI);
                JpetraObject.println("ERR", "Complex data type not supported by CisMatrix.");
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
                data = new double[size.numEntries()];
                mvr.readCoordinate(row, col, (double[]) data);
                buildCisMatrixFromSparse(row, col, (double[]) data, result);
            } else if (info.isComplex()) {
                /*dataR = new double[size.numEntries()];
                dataI = new double[size.numEntries()];
                mvr.readCoordinate(
                row,
                col,
                (double[]) dataR,
                (double[]) dataI);
                 shouldn't end up here.*/
            } else if (info.isPattern())
                mvr.readPattern(row, col);
            // !! need to add some support for this one!
            else
                throw new IOException("Parser error");
        }
        
        mvr.close();
        fis.close();
        
        return result;
    }
    
    public static void buildCisMatrixFromSparse(int[] row, int[] col, double[] data, CisMatrix cisMatrix) {
        if (cisMatrix.isRowOriented()) {
            for(int i=0; i < row.length; i++) {
                cisMatrix.insertEntry(row[i], col[i]-1, data[i]);
            }
        }
        else {
            for(int i=0; i < col.length; i++) {
                cisMatrix.insertEntry(col[i], row[i]-1, data[i]);
            }
        }
        cisMatrix.fillComplete();
    }
    
    public static void buildCisMatrixFromDense(double[] data, MatrixSize size, CisMatrix cisMatrix) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        int entryIndex = 0;
        if (cisMatrix.isRowOriented()) {
            for (int i=0; i < numCols; i++) {
                for(int j=0; j < numRows; j++) {
                    if (data[entryIndex] != 0) {
                        cisMatrix.insertEntry(j, i, data[entryIndex++]);
                    }
                }
            }
        }
        else {
            for (int i=0; i < numCols; i++) {
                for(int j=0; j < numRows; j++) {
                    if (data[entryIndex] != 0) {
                        cisMatrix.insertEntry(i, j, data[entryIndex++]);
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
