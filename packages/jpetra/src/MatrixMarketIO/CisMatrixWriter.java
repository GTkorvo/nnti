

package Jpetra.MatrixMarketIO;

import java.io.FileOutputStream;
import java.io.PrintWriter;

import Jpetra.*;

/**
 *
 * @author  Jason Cross
 */
public class CisMatrixWriter {
    
    public CisMatrixWriter() {
    }
    
    public static void write(String fileName, CisMatrix cisMatrix) throws java.io.IOException {
        PrintWriter out = new PrintWriter(new FileOutputStream(fileName));
        
        out.println("%%MatrixMarket matrix coordinate real general");
        out.println(cisMatrix.getNumRows() + " " + cisMatrix.getNumColumns() + " " + cisMatrix.getNumNonZeros());
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
