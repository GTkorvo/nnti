

package Jpetra.MatrixMarketIO;

import java.io.FileOutputStream;
import java.io.PrintWriter;

import Jpetra.*;

/**
 *
 * @author  Jason Cross
 */
public class GraphWriter {
    
    public GraphWriter() {
    }
    
    public static void write(String fileName, Graph graph, CisMatrix cisMatrix) throws java.io.IOException {
        PrintWriter out = new PrintWriter(new FileOutputStream(fileName));
        
        out.println("%%MatrixMarket matrix coordinate pattern general");
        out.println(cisMatrix.getNumRows() + " " + cisMatrix.getNumColumns() + " " + graph.getNumNonZeros());
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
