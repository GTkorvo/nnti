/*
 * Directory.java
 *
 * Created on June 5, 2001, 2:03 PM
 */

package Jpetra;

/**
 *
 * @author  mwboldt
 * @version 
 */
public class Directory extends JpetraObject {

    private BlockMap map = null;
    private Map directoryMap = null;
    private int [] procList = null;
    private int [] localIndexList = null;
    private int [] allMinGIDs = null;
    
    /** Creates new Directory */
    public Directory(BlockMap map) 
    //throws JpetraException 
    {
        this.map = map;
        
        // Test for simple cases
        
        // Uniprocess and local map cases (nothing to set up)
        if(!map.isDistributedGlobal()) return;
        
        // Linear map case
        else if(map.isLinearMap()) {
            // Build a list of the Minimum global ids for all processes on each process.
            // Since the map is linear, we know that all GIDs are contiguous on each process
            // and can be found using the MinGIDs.
            
            int numProc = map.getComm().getNumVnodes();
            
            allMinGIDs = new int [numProc+1];
            
             
            int [] minProcessGID = new int [1];
            minProcessGID[0] = map.getMinVnodeGID();
           
            
            int[] tmpallMinGIDs = map.getComm().gatherAll(minProcessGID);
            System.arraycopy(tmpallMinGIDs, 0, allMinGIDs, 0, tmpallMinGIDs.length);
            allMinGIDs[numProc] = 1 + map.getMaxAllGID(); // Set max cap
        }
            
            // General case.  Need to build a directory via calls to communication functions
            else {
                int flag = generate();
                //if(flag != 0) throw  new JpetraException("error in Directory(BlockMap map)");
                if(flag != 0) {
                    System.out.println("error in Directory(BlockMap)");
                    System.exit(1);
                }
            }
    }
    
    public Directory(Directory directory) 
    //throws JpetraException 
    {
        map = directory.map;
        
        int i;
	if(directory.directoryMap != null) directoryMap = new Map(directory.directoryMap);

        
        int dirNumProcessElements = directoryMap.getNumVnodeElements();
        
        if(directory.procList != null) {
            procList = new int [dirNumProcessElements];
            for(i=0; i<dirNumProcessElements; i++) 
                procList[i] = directory.procList[i];
        }
        if(directory.localIndexList != null) {
            localIndexList = new int [dirNumProcessElements];
            for(i=0; i<dirNumProcessElements; i++) 
                localIndexList[i] = directory.localIndexList[i];
        }
    }
    
    private int generate() {
        // Do other stuff for MPI
        return 0;
    }
    
    public int getDirectoryEntries(int numEntries, int [] globalEntries,
    int [] procs, int []localEntries, int [] entrySizes) {
        int ierr = 0;
        int i;
        int myPID = map.getComm().getVnodeID();
        int numProc = map.getComm().getNumVnodes();
        int nOverP = map.getNumGlobalElements() / numProc;
        int remainder = map.getNumGlobalElements() % numProc;
        
        // Test for simple cases
        
        // Uniprocesor and local map cases
        if(!map.isDistributedGlobal()) {
            for(i=0; i<numEntries; i++) {
                if(localEntries != null)
                    localEntries[i] = map.getLID(globalEntries[i]);
                // If GID is not valied, return -1 for both the proc and local entry info
                if(localEntries != null) if(localEntries[i] == -1)
                    procs[i] = -1;
                else procs[i] = myPID;
            }
            if(entrySizes != null) {
                if(map.hasConstantElementSize()) {
                    int elementSize = map.getMaxElementSize();
                    for(i=0; i<numEntries; i++) entrySizes[i] = elementSize;
                }
                else {
                    int [] elementSizeList = map.getElementSizeList();
                    for (i=0; i<numEntries; i++)
                        if(localEntries[i] > -1)
                            entrySizes[i] = elementSizeList[localEntries[i]];
                }
            }
            return(ierr);
        }


        /* MPI: int j; */
        /*int j;
        // MPI:
          int minAllGID = map.getMinAllGID(); // Get min of all GID
          int maxAllGID = map.getMaxAllGID(); // Get max of all GID
          for(i=0; i<numEntries; i++) {
              int LID = -1; // Assume not found
              int proc = -1;
              int GID = globalEntries[i];
              if(GID < minAllGID) ierr = -1;
              else if (GID > maxAllGID) ierr = -1;
              else {
                  // Guess uniform distribution and start a little above it
                  int proc1 = Math.min(GID/Math.max(nOverP,1)+2, numProc-1);
                  boolean found = false;
                  while(proc1 >= 0 && proc1 < numProc) {
                      if(allMinGIDs[proc1] <= GID) {
                          if(GID < allMinGIDs[proc1+1]) {
                              found = true;
                              break;
                          }
                          else proc1++;
                      }
                      else proc1--;
                  }
                  if(found) {
                      proc = proc1;
                      LID = GID - allMinGIDs[proc];
                  }
              }
              procs[i] = proc;
              if(localEntries != null) localEntries[i] = LID;
          }
          
          if(entrySizes != null) {
              if(map.hasConstantElementSize()) {
                  int elementSize = map.getMaxElementSize();
                  for(i=0; i<numEntries; i++) entrySizes[i] = elementSize;
              }
              else {
                  int [] elementSizeList = map.getElementSizeList(); // We know this exists
         
                  int sizeMsgTag = 22761;
                  int sizeMsgTag2 = 22760;
         
                  /*
                   * Not sure what this is... my best guess is something that was used before Distributor came along...
                   * GSComm_Plan * GSPlan = new GSComm_Plan();
                   * GSComm_Comm * GSComm = new GSComm_Comm();
                   
         
                  int sizeNumSends;
                  int [] sizeSendGIDs = null;
                  int [] sizeSendProcs = null;
         
                  /*boolean commFlag = Size_GSPlan.CreateFormRecvs(numEntries, globalEntries,
                      procs, map.getComm(), sizeNumSends, sizeSendGIDs, sizeSendProcs);
                   
                  
                  Distributor myDist = map.getComm().createDistributor();
                  
                  DistObject myDistObj = myDist.createFromReceives(numEntries, globalEntries, dirProcs, false);
          
                  sizeNumSends = myDistObj.getNumExportIDs();
                  sizeSendGIDs = myDistObj.getExportGIDs();
                  sizeSendProcs = myDistObj.getExportPIDs();
                  
                  int [] sizeExports = null;
                  int [] sizeImports = null;
                  if(sizeNumSends > 0) {
                      sizeExports = new int [2*sizeNumSends];
                      for(i=0; i<sizeNumSends; i++) {
                          int sizeCurrGID = sizeSendGIDs[i];
                          int sizeCurrLID = map.getLID(sizeCurrGID);
                          if(sizeCurrLID == -1) throw JpetraException("error in function generate");
                          sizeExports[2*i] = sizeCurrGID;
                          int sizeCurrSize = elementSizeList[sizeCurrLID];
                          sizeExports[2*i+1] = sizeCurrSize;
                      }
                  }
         
                  if(numEntries > 0) sizeImports = new int [2*numEntries];
                  /*Size_GSComm.Do(Size_GSPlan, sizeMsgTag2,
                      reinterpret_cast<char []> (sizeExports),
                      2 * sizeof(int), reinterpret_cast<char []> (sizeImports));
                   
		  
		  sizeImports = myDist.Do();
		  for(i=0; i<numEntries; i++) {
                      // Need to change !!!!
                      // boolean found = false
                      int sizeCurrLID = sizeImports[2*i];
                      for(j=0; j<numEntries; j++)
                          if(sizeCurrLID == globalEntries[j]) {
                              entrySizes[j] = sizeImports[2*i+1];
                              // found = true;
                              break;
                          }
                  }
              }
              return ierr;
          }
          
          // General case (need to set up an actual directory structure)
          int [] elemnetSizeList = null;
          int packetSize = 2   int PacketSize = 2; // We will send at least the GID and PID.  Might also send LID and Size info
          boolean doSizes = false;
          if(entrySizes != null) { 
              if(map.hasConstantElementSize()) {
                  int elementSize = map.getMaxElementSize();
                  for(i=0; i<numEntries; i++) entrySizes[i] = elementSize;
              }
              else {
                  elementSizeList = map.getElementSizeList(); // We know this exists
                  doSizes = true;
                  packetSize++; // Sending size info
              }
          }
         
          boolean doLIDs = (localEntries != null); // Do LIDs?
          if(doLIDs) packetSize++; // Sending LIDs also
         
          int msgTag = 22763;
          int msgTag2 = 22762;
         
         /*
          * Not sure what this is... my best guess is something that was used before Distributor came along...
          * GSComm_Plan * GSPlan = new GSComm_Plan();
          * GSComm_Comm * GSComm = new GSComm_Comm();
          
          
          Distributor myDist = map.getComm().createDistributor();
          
          int [] dirProcs = null;
          if(numEntries > 0) dirProcs = new int [numEntries];
         
          // Get directory locations for the requested list of entries
          directoryMap.getRemoteIDList(numEntries, globalEntries, dirProcs, 0);
          int numSends;
          int [] sendGIDs = null;
          int [] sendProcs = null;
         
          /*boolean commFlag = GSPlan.CreateFromRecvs(numEntries, globalEntries, dirProcs,
              directoryMap.getComm().Comm(), msgTag, true, numSends, sendGIDs, sendProcs);
           
           
          DistObject myDistObj = myDist.createFromReceives(numEntries, globalEntries, dirProcs, false);
          
          numSends = myDistObj.getNumExportIDs();
          sendGIDs = myDistObj.getExportGIDs();
          sendProcs = myDistObj.getExportPIDs();
          
          int currLID;
          int [] exports = null;
          int [] imports = null;
          if(numSends > 0) {
              exports = new int [packetSize*numSends];
              for(i=0; i<numSends; ;) {
                  int currGID = sendGIDs[i];
                  exports[i++] = currGID;
                  currLID = directoryMap.getLID(currGID);
                  if(currLID == -1) throw new JpetraException("error in");
                  exports[i+1] = procList[currLID];
                  if(doLIDs) exports[2+i++] = localIndexList[currLID];
                  if(doSizes) exports[2+i++] = elementSizeList[currLID];
              }
          }
         
          if(numEntries > 0) imports = new int [packetSize*numEntries];
          
          GSComm.Do(GSPlan, msgTag2, reinterpret_cast<char *> (exports),
              pcketSize * sizeof(int), reinterpret_cast<char *> (imports));
          
          
          for(i=0; i<numEntries; i++) {
              // boolean found = false;
              currLID = imports[i++];
              for(j=0; j<numEntries; j++)
                  if(currLID == globalEntries[j]) {
                      procs[j] = imports[i++];
                      if(doLIDs) localEntries[j] = imports[j++];
                      if(doSizes) localEntries[j] = imports[j++];
                      // found = true;
                      break;
                  }
                  // if (!found) cout << "Internal error:  Petra_Directory::GetDirectoryEntries: Global Index " << curr_LID
                  //	     << " not on process " << MyPID << endl; abort();
         
                  return (commFlag?1:0);
         */
         return 0;
    }
}
