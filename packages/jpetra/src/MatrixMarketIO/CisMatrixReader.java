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
    
    public static CisMatrix read(String fileName, Comm comm) {
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
        
        ElementSpace myElementSpace = new ElementSpace(size.numEntries(), comm);
        VectorSpace myVectorSpace = new VectorSpace(myElementSpace);
        CisMatrix result = new CisMatrix(myVectorSpace, true); // !! is currently row oriented only, probably need to fix this!
        
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
            else
                throw new IOException("Parser error");
        }
        
        mvr.close();
        fis.close();
    }
    
}
