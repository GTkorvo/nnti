class MapTest {

    public static void main(String [] args) {
	int size = 1;
	int rank = 0;
	boolean verbose = false;
	
	if(args.length > 0 && args[0].equals("-v")) verbose = true;
	
	Jpetra.SerialComm comm = new Jpetra.SerialComm();
	int pid = comm.getVnodeID();
	int numProc = comm.getNumVnodes();
	if(verbose) System.out.println("Processor "+pid+" of "+numProc+" is alive.");
	boolean verbose1 = verbose;
	verbose = (pid == 0);
	int numNodeElements = 10000;
	int numNodeElements1 = numNodeElements;
	int numGlobalElements = numNodeElements * numProc + Math.min(numProc, 3);
	if(pid < 3) numNodeElements++;
	int indexBase = 0;
	boolean isDistributedGlobal = (numGlobalElements > numNodeElements);

	// Test Jpetra defined uniform linear distribution constructor
	Jpetra.Map map = new Jpetra.Map(numGlobalElements, indexBase, comm);
	if(verbose) System.out.println("1) Checking Map(numGlobalElements, indexBase, comm)");
	int ierr = checkmap(map, numGlobalElements, numNodeElements, null,
			    indexBase, comm, isDistributedGlobal);
	if(verbose) {
	    if(ierr == 0) System.out.println("Checked OK");
	    else System.out.println("Error code: "+ierr);
	    System.out.println();
	}
	if(ierr != 0) System.exit(0);
	
	// Test user defined linear distribution constructor
	map = new Jpetra.Map(numGlobalElements, numNodeElements, indexBase, comm);
	if(verbose) System.out.println("2) Checking Map(numGlobalElements, numNodeElements, indexBae, comm)");
	ierr = checkmap(map, numGlobalElements, numNodeElements, null,
			indexBase, comm, isDistributedGlobal);
	if(verbose) {
	    if(ierr == 0) System.out.println("Checked OK");
	    else System.out.println("Error code: "+ierr);
	    System.out.println();
	}
	if(ierr != 0) System.exit(0);

	// Test user defined arbitrary distribution constructor
	// Generate global element list.
	int [] myGlobalElements = new int [numNodeElements];
	int maxGid = (comm.getVnodeID()+1)*numNodeElements-1+indexBase;
	if(comm.getVnodeID() > 2) maxGid += 3;
	for(int i=0; i<numNodeElements; i++) myGlobalElements[i] = maxGid-i;

	map = new Jpetra.Map(numGlobalElements, numNodeElements, myGlobalElements,
			     indexBase, comm);
	if(verbose) System.out.println("3) Cheking Map(numGlobalElements, numNodeElements, myGlobalElements, indexBase, comm)");
	ierr = checkmap(map, numGlobalElements, numNodeElements, myGlobalElements,
			indexBase, comm, isDistributedGlobal);
	if(verbose) {
	    if(ierr == 0) System.out.println("Checked OK");
	    else System.out.println("Error code: "+ierr);
	    System.out.println();
	}
	if(ierr != 0) System.exit(0);

	// Test LocalMap constructor
	Jpetra.LocalMap localMap = new Jpetra.LocalMap(numNodeElements1, indexBase, comm);
	if(verbose) System.out.println("4) Checking LocalMap(numNodeElements1, indexBase, comm)");
	ierr = checkmap(localMap, numNodeElements1, numNodeElements1, null,
			indexBase, comm, false);
	if(verbose) {
	    if(ierr == 0) System.out.println("Checked OK");
	    else System.out.println("Error code: "+ierr);
	    System.out.println();
	}
	if(ierr != 0) System.exit(0);

	System.exit(0);
    }

    public static int checkmap(Jpetra.Map map, int numGlobalElements, int numNodeElements,
			int [] globalElements, int indexBase, Jpetra.SerialComm comm,
			boolean isDistributedGlobal) {
	int i;
	
	if(!map.hasConstantElementSize()) return -1;
	if(isDistributedGlobal != map.isDistributedGlobal()) return -3;
	if(map.getElementSize() != 1) return -4;

	int [] myElementSizeList = new int [numNodeElements];
	if(map.getElementSizeList(myElementSizeList) != 0) return -5;
	for(i=0; i<numNodeElements; i++) {
	    if(myElementSizeList[i] != 1) return -5;
	}

	Jpetra.SerialComm comm1 = (Jpetra.SerialComm)map.getComm();

	if(comm1.getNumVnodes() != comm.getNumVnodes()) return -6;
	if(comm1.getVnodeID() != comm.getVnodeID()) return -7;
	if(map.getIndexBase() != indexBase) return -8;
	if(!map.isLinearMap() && globalElements == null) return -9;
	if(map.isLinearMap() && globalElements != null) return -9;
	if(map.getMaxAllGID() != numGlobalElements-1+indexBase) return -10;
	if(map.getMaxElementSize() != 1) return -11;
	if(map.getMaxLID() != numNodeElements-1+indexBase) return -12;

	int maxGid = (comm.getVnodeID()+1)*numNodeElements-1+indexBase;
	if(comm.getVnodeID() > 2) maxGid += 3;
	if(!isDistributedGlobal) maxGid = numNodeElements-1+indexBase;
	if(map.getMaxNodeGID() != maxGid) return -13;
	if(map.getMinAllGID() != indexBase) return -14;
	if(map.getMinElementSize() != 1) return -15;
	if(map.getMinLID() != indexBase) return -16;

	int minGid = comm.getVnodeID()*numNodeElements+indexBase;
	if(comm.getVnodeID() > 2) minGid += 3;
	if(!isDistributedGlobal) minGid = 0;
	if(map.getMinNodeGID() != minGid) return -17;

	int [] myGlobalElements1 = new int [numNodeElements];
	if(map.getGlobalElements(myGlobalElements1) != 0) return -18;

	if(globalElements == null) {
	    for(i=0; i<numNodeElements; i++) {
		if(myGlobalElements1[i] != minGid+i) return -18;
	    }
	}
	else {
	    for(i=0; i<numNodeElements; i++) {
		if(globalElements[i] != myGlobalElements1[i]) return -18;
	    }
	}

	if(map.getNumGlobalElements() != numGlobalElements) return -19;
	if(map.getNumGlobalEquations() != numGlobalElements) return -20;
	if(map.getNumNodeElements() != numNodeElements) return -21;
	if(map.getNumNodeEquations() != numNodeElements) return -22;

	return 0;
    }
}
