import Jpetra.*;
import java.util.Arrays;

public class CcjCommGatherAllTest {
    boolean verbose = false;
    
    public static void main (String[] args) {
        new CcjCommGatherAllTest(args);
     }
     
     public CcjCommGatherAllTest(String[] args) {
         Jpetra.CcjComm comm = new Jpetra.CcjComm("ccjhosts.txt");
         
         if(args.length > 0 && args[0].equals("-v")) verbose = true;
      
         if(verbose) System.out.println("Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
         
         int offSet = comm.getVnodeID()*4;
         int[] myIntArray = new int[]{offSet, offSet+1, offSet+2, offSet+3};
         double[] myDoubleArray = new double[]{offSet, offSet+1, offSet+2, offSet+3};        
         
         int myInt = comm.getVnodeID();
         double myDouble = myInt;
         int[] finalInts = new int[comm.getNumVnodes()];
         double[] finalDoubles = new double[comm.getNumVnodes()];
         for (int i=0; i < comm.getNumVnodes(); i++) {
             finalInts[i]=i;
             finalDoubles[i]=i;
         }
         
         int[] finalIntArray = new int[comm.getNumVnodes()*4];
         double[] finalDoubleArray = new double[comm.getNumVnodes()*4];
         for (int i = 0; i < comm.getNumVnodes()*4; i++) {
             finalIntArray[i] = i;
             finalDoubleArray[i] = i;
         }
         
         
         int[] allIntElements = comm.gatherAll(myIntArray);
         double[] allDoubleElements = comm.gatherAll(myDoubleArray);
         int[] allInts = comm.gatherAll(myInt);
         double[] allDoubles = comm.gatherAll(myDouble);
         
         if (!java.util.Arrays.equals(allIntElements, finalIntArray)) {
             System.err.println("intArrayGatherAllTest Failed!");
         }
 
         if (!java.util.Arrays.equals(allDoubleElements, finalDoubleArray)) {
             System.err.println("doubleArrayGatherAllTest Failed!");
         }
         
         if (!java.util.Arrays.equals(allInts, finalInts)) {
             System.err.println("intGatherAllTest Failed!");
         }
         
         if (!java.util.Arrays.equals(allDoubles, finalDoubles)) {
             System.err.println("doubleGatherAllTest Failed!");
         }
         
         
         if (verbose) {
             for(int i=0; i < allIntElements.length; i++) {
                 System.out.println(allIntElements[i] + " " + allDoubleElements[i]);    
             }
             
             System.out.println("-------");
             
             for(int i=0; i < allInts.length; i++) {
                 System.out.println(allInts[i] + " " + allDoubles[i]);
             }
         }
         
         System.exit(1);
    }
}