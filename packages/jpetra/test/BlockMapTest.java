class BlockMapTest {

    public static void main(String [] args) {
        
        int i;
    
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
    
        int numVnodeElements = 10000;
        int numGlobalElements = numVnodeElements*numProc+Math.min(numProc, 3);
        if(pid < 3) numVnodeElements++;
        int indexBase = 0;
        int elementSize = 7;
        boolean isDistributedGlobal = numGlobalElements > numVnodeElements;
    
        // Test uniform linear distribution constructor
        Jpetra.BlockMap map = new Jpetra.BlockMap(numGlobalElements, elementSize, indexBase, comm);
        
        if(verbose) {
            System.out.println("1) Checking BlockMap(numGlobalElements, elementSize, indexBase, comm)");
        }
        
        int ierr = checkmap(map, numGlobalElements, numVnodeElements, null, elementSize, null,
                            numGlobalElements*elementSize, numVnodeElements*elementSize,
                            indexBase, comm, isDistributedGlobal);
        if(verbose) {
            if(ierr == 0) {
                System.out.println("Checked OK");
            }
            else {
                System.out.println("Error code: "+ierr);
            }
            System.out.println();
        }
        
        if(ierr != 0) {
            System.exit(0);
        }
        
        // Test user defined linear distribution constructor
        map = new Jpetra.BlockMap(numGlobalElements, numVnodeElements, elementSize, indexBase, comm);
        
        if(verbose) {
            System.out.println("2) Checking BlockMap(numGlobalElements, numVnodeElements, elementSize, indexBase, comm)");
        }
        
        ierr = checkmap(map, numGlobalElements, numVnodeElements, null, elementSize, null,
                        numGlobalElements*elementSize, numVnodeElements*elementSize,
                        indexBase, comm, isDistributedGlobal);
    
        if(verbose) {
            if(ierr == 0) {
                System.out.println("Checked OK");
            }
            else {
                System.out.println("Error code: "+ierr);
            }
            System.out.println();
        }
    
        // Test user defined arbitrary distribution constructor
        // Generate global element list.
        int [] globalElements = new int [numVnodeElements];
        int maxGid = (comm.getVnodeID()+1)*numVnodeElements-1+indexBase;
        
        for(i=0; i<numVnodeElements; i++) {
            globalElements[i] = maxGid-i;
        }
    
        map = new Jpetra.BlockMap(numGlobalElements, numVnodeElements, globalElements, elementSize, indexBase, comm);
        
        if(verbose) {
            System.out.println("3) Checking BlockMap(numGlobalElements, numVnodeElements, globalElements, elementSize, indexBase, comm)");
        }
        
        ierr = checkmap(map, numGlobalElements, numVnodeElements, globalElements, elementSize, null,
                        numGlobalElements*elementSize, numVnodeElements*elementSize,
                        indexBase, comm, isDistributedGlobal);
        
        if(verbose) {
            if(ierr == 0) {
                System.out.println("Checked OK");
            }
            else {
                System.out.println("Error code: "+ierr);
            }
            System.out.println();
        }
        
        if(ierr != 0) {
            System.exit(0);
        }
    
        int [] elementSizeList = new int [numVnodeElements];
        int numVnodeEquations = 0;
        int numGlobalEquations = 0;
        
        for(i=0; i<numVnodeElements; i++) {
            elementSizeList[i] = i%6+2; // blocksizes go from 2 to 7
            numVnodeEquations += elementSizeList[i];
        }
        
        elementSize = 7; // Set to maximum for use in checkmap
        numGlobalEquations = comm.getNumVnodes()*numVnodeEquations;
    
        // Adjust numGlobalEquations based on processor ID
        if(comm.getNumVnodes() > 3) {
            if(comm.getVnodeID() > 2) {
                numGlobalEquations += 3*((numVnodeElements)%6+2);
            }
            else {
                numGlobalEquations -= (comm.getNumVnodes()-3)*((numVnodeElements-1)%6+2);
            }
        }
    
        for(int k=0; k<elementSizeList.length; k++) if(elementSizeList[k] <= 0) System.out.println(k+" "+elementSizeList[k]);
        map = new Jpetra.BlockMap(numGlobalElements, numVnodeElements, globalElements, elementSizeList,
                      indexBase, comm);
        if(verbose) System.out.println("4) Checking BlockMap(numGlobalElements, numVnodeElements, globalElements, elementSizeList, indexBase, comm)");
        ierr = checkmap(map, numGlobalElements, numVnodeElements, globalElements, elementSize, elementSizeList,
                numGlobalEquations, numVnodeEquations,
                indexBase, comm, isDistributedGlobal);
        if(verbose) {
            if(ierr == 0) System.out.println("Checked OK");
            else System.out.println("Error code: "+ierr);
            System.out.println();
        }
        if(ierr != 0) System.exit(0);
    
        // Test copy constructor
        Jpetra.BlockMap map1 = new Jpetra.BlockMap(map);
    
        if(verbose) System.out.println("5) Checking BlockMap(map)");
        ierr = checkmap(map1, numGlobalElements, numVnodeElements, globalElements, elementSize, elementSizeList,
                numGlobalEquations, numVnodeEquations,
                indexBase, comm, isDistributedGlobal);
        if(verbose) {
            if(ierr == 0) System.out.println("Checked OK");
            else System.out.println("Error code: "+ierr);
            System.out.println();
        }
        if(ierr != 0) System.exit(0);
    
        System.exit(0);
    }

    public static int checkmap(Jpetra.BlockMap map, int numGlobalElements, int numVnodeElements,
                               int [] globalElements, int elementSize, int [] elementSizeList,
                               int numGlobalEquations, int numVnodeEquations,
                               int indexBase, Jpetra.Comm comm, boolean isDistributedGlobal) {
    
        int i;
        
        if(elementSizeList == null) {
            if(!map.hasConstantElementSize()) return -1;
        }
        else {
            if(map.hasConstantElementSize()) return -1;
        }
        
        if(isDistributedGlobal != map.isDistributedGlobal()) return -3;
    
        int [] myElementSizeList;
    
        if(elementSizeList == null) {
            if(map.getElementSize() != elementSize) return -4;
            
            myElementSizeList = new int [numVnodeElements];
    
            if(map.getElementSizeList(myElementSizeList) != 0) return -5;
            for(i=0; i<numVnodeElements; i++) {
            if(myElementSizeList[i] != elementSize) return -5;
            }
        }
        else {
            myElementSizeList = new int[numVnodeElements];
            if(map.getElementSizeList(myElementSizeList) != 0) return -5;
            
            for(i=0; i<numVnodeElements; i++) {
            if(myElementSizeList[i] != elementSizeList[i]) return -5;
            }
        }
    
        Jpetra.SerialComm comm1 = (Jpetra.SerialComm)map.getComm();
        
        if(comm1.getNumVnodes() != comm.getNumVnodes()) return -6;
        if(comm1.getVnodeID() != comm.getVnodeID()) return -7;
        if(map.getIndexBase() != indexBase) return -8;
        if(!map.isLinearMap() && globalElements == null) return -9;
        if(map.isLinearMap() && globalElements != null) return -9;
        if(map.getMaxAllGID() != numGlobalElements-1+indexBase) return -10;
        if(map.getMaxElementSize() != elementSize) return -11;
        if(map.getMaxLID() != numVnodeElements-1+indexBase) return -12;
    
        int maxGid = (comm.getVnodeID()+1)*numVnodeElements-1+indexBase;
        if(comm.getVnodeID() > 2) maxGid += 3;
        if(!isDistributedGlobal) maxGid = numVnodeElements-1+indexBase;
        if(map.getMaxMyGID() != maxGid) return -13;
    
        if(map.getMinAllGID() != indexBase) return -14;
    
        if(elementSizeList == null) {
            if(map.getMinElementSize() != elementSize) return -15;
        }
        else if(map.getMinElementSize() != 2) return -15;
    
        if(map.getMinLID() != indexBase) return -16;
    
        int minGid = comm.getVnodeID()*numVnodeElements+indexBase;
        if(comm.getVnodeID() > 2) minGid += 3;
        if(!isDistributedGlobal) minGid = 0;
        if(map.getMinMyGID() != minGid) return -17;
    
        int [] globalElements1 = new int[numVnodeElements];
        globalElements1 = map.getMyGlobalElements();
        if(globalElements == null) {
            for(i=0; i<numVnodeElements; i++) {
            if(globalElements1[i] != minGid+i) return -18;
            }
        }
        else {
            for(i=0; i<numVnodeElements; i++) {
            if(globalElements1[i] != globalElements[i]) return -18;
            }
        }
    
        if(map.getNumGlobalElements() != numGlobalElements) return -19;
        if(map.getNumGlobalEquations() != numGlobalEquations) return -20;
        if(map.getNumMyElements() != numVnodeElements) return -21;
        if(map.getNumMyEquations() != numVnodeEquations) return -22;
    
        if(map.isLinearMap()) {
            int [] gidList = new int [3];
            int [] pidList = new int [3];
            int [] lidList = new int [3];
            int pid = map.getComm().getVnodeID();
    
            int numIDs = 0;
            //gidList[numIDs++] = map.maxAllGID()+1; // Should return -1 for both pid and lid
            if(map.getMinMyGID()-1 >= map.getMinAllGID()) gidList[numIDs++] = map.getMinMyGID()-1;
            if(map.getMaxMyGID()+1 <= map.getMaxAllGID()) gidList[numIDs++] = map.getMaxMyGID()+1;
    
            map.getRemoteIDList(numIDs, gidList, pidList, lidList);
    
            numIDs = 0;
    
            if(map.getMinMyGID()-1 >= map.getMinAllGID()) {
            if(pidList[numIDs++] != pid-1) {
                return -23;
            }
            }
            if(map.getMaxMyGID()+1 <= map.getMaxAllGID()) {
            if(pidList[numIDs] != pid-1) {
                return -24;
            }
            }
            if(map.getMaxMyGID()+1 <= map.getMaxAllGID()) {
            if(lidList[numIDs++] != 0) {
                return -25;
            }
            }
        }
    
        return 0;
    }

    
}
