import Jpetra.*;

public class IntSDMTest {
    public static void main(String[] args) {
        new IntSDMTest();
    }

    public IntSDMTest() {
        int numRows = 5;
        int numCols = 5;
        IntSerialDenseMatrix isdm = new IntSerialDenseMatrix(numRows, numCols);
        
        // fill up isdm
        int i, j;
        int k = 0;
        for(i=0; i < numRows; i++) {
            for(j=0; j < numCols; j++) {
                isdm.setElement(i, j, k++);
            }
        }
        
        isdm.print();
        
        IntSerialDenseMatrix isdm1 = new IntSerialDenseMatrix(isdm);
        
        isdm1.print();
        
        System.out.println("isdm OneNorm: " + isdm.getOneNorm());
        System.out.println("isdm1 OneNorm: " + isdm1.getOneNorm());
        
        System.out.println("isdm InfoNorm: " + isdm.getInfNorm());
        System.out.println("isdm1 InfoNorm: " + isdm1.getInfNorm());
        
        isdm1.reshape(numRows+2, numCols+2);
        
        isdm.print();
        isdm1.print();
        
        System.out.println("isdm OneNorm: " + isdm.getOneNorm());
        System.out.println("isdm1 OneNorm: " + isdm1.getOneNorm());
        
        System.out.println("isdm InfoNorm: " + isdm.getInfNorm());
        System.out.println("isdm1 InfoNorm: " + isdm1.getInfNorm());
        
        isdm1.reshape(numRows-2, numCols-2);
        isdm.shape(numRows-3, numCols-3);
        
        isdm.print();
        isdm1.print();
        
        System.out.println("isdm OneNorm: " + isdm.getOneNorm());
        System.out.println("isdm1 OneNorm: " + isdm1.getOneNorm());
        
        System.out.println("isdm InfoNorm: " + isdm.getInfNorm());
        System.out.println("isdm1 InfoNorm: " + isdm1.getInfNorm());
    }
}