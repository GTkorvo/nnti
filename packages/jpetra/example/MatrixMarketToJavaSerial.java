package example;

import Jpetra.*;
import Jpetra.MatrixMarketIO.*;

/**
 *
 * @author  Jason Cross
 */
public class MatrixMarketToJavaSerial extends JpetraObject {
    
    /** Creates a new instance of MatrixMarketToJavaSerial */
    public MatrixMarketToJavaSerial(String inFileName, String outFileName) {
        initializeOutput();
        Comm comm = new SerialComm();
        
        Time time = new Time(comm);
        CisMatrix inMatrix = null;
        try {
            inMatrix = CisMatrixReader.read(inFileName, CisMatrix.ROW_ORIENTED, comm);
        }
        catch (java.io.IOException e) {
            this.println("FATALERR", e.toString());
            System.exit(1);
        }
        CisMatrix.writeToFile(outFileName, inMatrix);
        this.println("STD", "total time to read/write the matrix: " + time.getElapsedTime());
    }
    
    public static void main(String args[]) {
        if(args.length != 2) {
            System.out.println("Usage:");
            System.out.println("java example/MatrixMarketToJavaSerial [inMatrixFileName.mtx] [outMatrixFileName]");
            System.exit(0);
        }
        new MatrixMarketToJavaSerial(args[0], args[1]);
    }
}
