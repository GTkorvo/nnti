import Jpetra.*;
import java.util.Arrays;

public class CcjCommScanSumTest {
    boolean verbose = false;
    
    public static void main (String[] args) {
        new CcjCommScanSumTest(args);
     }
     
     public CcjCommScanSumTest(String[] args) {
         Jpetra.CcjComm comm = new Jpetra.CcjComm("ccjhosts.txt");
         
         if(args.length > 0 && args[0].equals("-v")) verbose = true;
      
         if(verbose) System.out.println("Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
             
         int[] myIntArray;
         double[] myDoubleArray;
         if (comm.getVnodeID()==0) {
            myIntArray = new int[]{5,7,2,1000};
            myDoubleArray = new double[]{5,7,2,1000};
         }
         else if (comm.getVnodeID()==1) {
            myIntArray = new int[]{2,1,98,3};
            myDoubleArray = new double[]{2,1,98,3};
         }
         else {
            myIntArray = new int[]{-2,4,102,9};
            myDoubleArray = new double[]{-2,4,102,9};
         }
         
         int[] myIntScanSum = comm.scanSums(myIntArray);
         double[] myDoubleScanSum = comm.scanSums(myDoubleArray);
                 
         
         if (verbose) {
             for(int i=0; i < myIntScanSum.length; i++) {
                 System.out.println(myIntScanSum[i] + " " + myDoubleScanSum[i]);   
             }       
         }
         
         System.exit(1);        
     }  
}