import Jpetra.*;

class SerialCommTest {

    public static void main (String [] args) {
	
	int size = 1;
	int rank = 0;
	boolean verbose = false;
	
	if(args.length > 0 && args[0].equals("-v")) verbose = true;

	Jpetra.SerialComm comm = new Jpetra.SerialComm();
	int myPid = comm.getPID();
	int numProce = comm.getNumProc();
	if(verbose) System.out.println("Processor "+comm.getPID()+" of "+comm.getNumProc());

	System.exit(0);
    }
}
