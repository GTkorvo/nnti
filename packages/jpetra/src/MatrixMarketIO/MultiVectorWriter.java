

package Jpetra.MatrixMarketIO;

import java.io.FileOutputStream;
import java.io.PrintWriter;

import Jpetra.*;

/**
 *
 * @author  Jason Cross
 */
public class MultiVectorWriter {
    
    public MultiVectorWriter() {
    }
    
    public static void write(String fileName, MultiVector multiVector) throws java.io.IOException {
        PrintWriter out = new PrintWriter(new FileOutputStream(fileName));
        
        out.println("%%MatrixMarket matrix array real general");
        out.println(multiVector.getNumRows() + " " + multiVector.getNumCols());
        out.println();
        
        double[][] values = multiVector.getValues();
        for(int col=0; col < values.length; col++) {
            for(int row=0; row < values[0].length; row++) {
                out.println(values[col][row]);
            }
        }
        out.close();
    }
}
