package Jpetra.MatrixMarketIO;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;

import Jpetra.*;

public class MultiVectorReader extends JpetraObject {
    
    public MultiVectorReader() {
    }
    
    public static MultiVector read(String fileName, Comm comm) throws java.io.IOException {
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
        MultiVector result = null;
        
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
                result = buildMultiVectorFromDense(size, (double[]) data, myVectorSpace);
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
                result = buildMultiVectorFromSparse(row, col, size, (double[]) data, myVectorSpace);
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
    
    public static MultiVector buildMultiVectorFromSparse(int[] row, int[] col, MatrixSize size, double[] data, VectorSpace myVectorSpace) {
        int numRows = size.numRows();
        int numCols = size.numColumns();
        double[][] values = new double[numCols][numRows];
        for(int i=0; i < size.numEntries(); i++) {
            values[col[i]-1][row[i]-1] = data[i];
        }
        
        return new MultiVector(myVectorSpace, values);
    }
    
    public static MultiVector buildMultiVectorFromDense(MatrixSize size, double[] data, VectorSpace myVectorSpace) {
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
