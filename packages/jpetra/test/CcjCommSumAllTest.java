import Jpetra.*;
import java.util.Arrays;

public class CcjCommSumAllTest {
    boolean verbose = false;
    
    public static void main (String[] args) {
        new CcjCommSumAllTest(args);
     }
     
     public CcjCommSumAllTest(String[] args) {
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

         int[] realIntSum;
         double[] realDoubleSum;         
         if (comm.getNumVnodes() == 1) {
             realIntSum = new int[]{5,7,2,1000};
             realDoubleSum = new double[]{5,7,2,1000};
         }
         else if (comm.getNumVnodes() == 2) {
             realIntSum = new int[]{7,8,100,1003};
             realDoubleSum = new double[]{7,8,100,1003};                 
         }
         else { 
             realIntSum = new int[]{5,12,202,1012};
             realDoubleSum = new double[]{5,12,202,1012};             
         }
         
         int[] myIntSum = comm.sumAll(myIntArray);
         double[] myDoubleSum = comm.sumAll(myDoubleArray);
                 
         
         if (!java.util.Arrays.equals(myIntSum, realIntSum)) {
             System.err.println("intSum Failed!");
         }
         if (!java.util.Arrays.equals(myDoubleSum, realDoubleSum)) {
             System.err.println("doubleSum Failed!");
         }
         
         if (verbose) {
             for(int i=0; i < myIntSum.length; i++) {
                 System.out.println(myIntSum[i] + " " + myDoubleSum[i]);   
             }       
         }
         
         System.exit(1);        
     }  
}