import Jpetra.*;

public class CcjCommTest extends JpetraObject {
    public CcjCommTest(String[] args) {
        if(args.length > 0 && args[0].equals("-v")) setRootPrint("VERBOSE", true);
    
        int size = 1;
        int rank = 0;
    
        Jpetra.CcjComm comm = new Jpetra.CcjComm("ccjhosts.txt");
        
        int myPid = comm.getVnodeID();
        int numProc = comm.getNumVnodes();
        println("VERBOSE", "Processor "+myPid+" of "+numProc);
    
        System.exit(0); 
    }
    
    public static void main (String[] args) {
        new CcjCommTest(args);  
    }
}