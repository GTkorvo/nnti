import Jpetra.*;

class SerialCommTest {

    public static void main (String [] args) {
	
	int size = 1;
	int rank = 0;
	boolean verbose = false;
	
	if(args.length > 0 && args[0].equals("-v")) verbose = true;

	Jpetra.SerialComm comm = new Jpetra.SerialComm();
	int myPid = comm.getVnodeID();
	int numProce = comm.getNumVnodes();
	if(verbose) System.out.println("Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());

	System.exit(0);
    }
}
