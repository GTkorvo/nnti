import Jpetra.*;

public class CcjCommTest {
    public static void main (String[] args) {
    
        int size = 1;
        int rank = 0;
        boolean verbose = false;
        
        if(args.length > 0 && args[0].equals("-v")) verbose = true;
    
        Jpetra.CcjComm comm = new Jpetra.CcjComm("ccjhosts.txt");
        
        int myPid = comm.getVnodeID();
        int numProc = comm.getNumVnodes();
        if(verbose) System.out.println("Processor "+myPid+" of "+numProc);
    
        System.exit(0);   
    }
}