import Jpetra.*;
import java.util.Arrays;

public class CcjCommSumAllTest2 {
    boolean verbose = false;
    
    public static void main (String[] args) {
        new CcjCommSumAllTest2(args);
     }
     
     public CcjCommSumAllTest2(String[] args) {
         Jpetra.CcjComm comm = new Jpetra.CcjComm("ccjhosts.txt");
         
         if(args.length > 0 && args[0].equals("-v")) verbose = true;
      
         if(verbose) System.out.println("Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
         
         int offSet = comm.getVnodeID()*4;
         int[] myIntArray = new int[]{offSet, offSet+1, offSet+2, offSet+3};
         double[] myDoubleArray = new double[]{offSet, offSet+1, offSet+2, offSet+3};  
         
         int[] finalIntArray = new int[4];
         double[] finalDoubleArray = new double[4];
         int j=0;
         for (int i = 0; i < comm.getNumVnodes(); i++) {
             j=i*4;
             
             finalIntArray[0] += j;
             finalDoubleArray[0] += j;
             finalIntArray[1] += j+1;
             finalDoubleArray[1] += j+1;
             finalIntArray[2] += j+2;
             finalDoubleArray[2] += j+2;            
             finalIntArray[3] += j+3;
             finalDoubleArray[3] += j+3;            
         }
         
         int[] globalIntSums = comm.sumAll(myIntArray);
         double[] globalDoubleSums = comm.sumAll(myDoubleArray);
         
         if (!java.util.Arrays.equals(globalIntSums, finalIntArray)) {
             System.err.println("intSumAll Failed!");
         }
         if (!java.util.Arrays.equals(globalDoubleSums, finalDoubleArray)) {
             System.err.println("doubleSumAll Failed!");
         }
         
         if (verbose) {
             for(int i=0; i < globalIntSums.length; i++) {
                 System.out.println(globalIntSums[i] + " " + globalDoubleSums[i]);   
             }
             System.out.println("------");
             for(int i=0; i < finalIntArray.length; i++) {
                System.out.println(finalIntArray[i] + " " + finalDoubleArray[i]);
             }        
         }
         
         System.exit(1);        
     }  
}