import Jpetra.*;
import java.util.Arrays;

public class CcjCommMinMaxAllTest {
    boolean verbose = false;
    
    public static void main (String[] args) {
        new CcjCommMinMaxAllTest(args);
     }
     
     public CcjCommMinMaxAllTest(String[] args) {
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
         
         int[] realIntMax;
         int[] realIntMin;
         double[] realDoubleMax;
         double[] realDoubleMin;
         
         if (comm.getNumVnodes() == 1) {
             realIntMax = new int[]{5,7,2,1000};
             realDoubleMax = new double[]{5,7,2,1000};
             
             realIntMin = new int[]{5,7,2,1000};
             realDoubleMin = new double[]{5,7,2,1000};
         }
         else if (comm.getNumVnodes() == 2) {
             realIntMax = new int[]{5,7,98,1000};
             realDoubleMax = new double[]{5,7,98,1000};
             
             realIntMin = new int[]{2,1,2,3};
             realDoubleMin = new double[]{2,1,2,3};                  
         }
         else { 
             realIntMax = new int[]{5,7,102,1000};
             realDoubleMax = new double[]{5,7,102,1000};
             
             realIntMin = new int[]{-2,1,2,3};
             realDoubleMin = new double[]{-2,1,2,3};               
         }
         
         int[] myIntMax = comm.maxAll(myIntArray);
         double[] myDoubleMax = comm.maxAll(myDoubleArray);
         
         int[] myIntMin = comm.minAll(myIntArray);
         double[] myDoubleMin = comm.minAll(myDoubleArray);        
         
         if (!java.util.Arrays.equals(myIntMax, realIntMax)) {
             System.err.println("intMax Failed!");
         }
         if (!java.util.Arrays.equals(myIntMin, realIntMin)) {
             System.err.println("intMin Failed!");
         }
         if (!java.util.Arrays.equals(myDoubleMax, realDoubleMax)) {
             System.err.println("doubleMax Failed!");
         }
         if (!java.util.Arrays.equals(myDoubleMin, realDoubleMin)) {
             System.err.println("doubleMin Failed!");
         }
         if (verbose) {
             for(int i=0; i < myIntMax.length; i++) {
                 System.out.println(myIntMax[i] + " " + myDoubleMax[i] + " " + myIntMin[i] + " " + myDoubleMin[i]);   
             }
         }        
         
         System.exit(1);        
     }  
}