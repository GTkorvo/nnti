// @HEADER
// ***********************************************************************
// 
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

package Jpetra;

public class Directory extends JpetraObject {

    private BlockMap map = null;
    private Map DirectoryMap = null;
    private int [] ProcList = null;
    private int [] LocalIndexList = null;
    private int [] SizeList = null;
    private int [] AllMinGIDs = null;
    
    public Directory(BlockMap map) {
        this.map = map;
        this.map = map;
        
        // Test for simple cases
        
        // Uniprocess and local map cases (nothing to set up)
        if(!map.isDistributedGlobal()) {
            System.out.println("Map is local.");
            return;
        }
        // Linear map case
        else if(map.isLinearMap()) {
            // Build a list of the Minimum global ids for all processes on each process.
            // Since the map is linear, we know that all GIDs are contiguous on each process
            // and can be found using the MinGIDs.
            System.out.println("Map is linear.");
            int numVnodes = map.getComm().getNumVnodes();
            
            int minMyGID = map.getMinMyGID();
            int[] tmpAllMinGIDs = map.getComm().gatherAll(minMyGID); 
                
            this.AllMinGIDs = new int [numVnodes+1];
            System.arraycopy(tmpAllMinGIDs, 0, this.AllMinGIDs, 0, tmpAllMinGIDs.length);
            
            /*
            debug code
            for(int i=0; i<this.AllMinGIDs.length; i++) {
                System.out.println("AllMinGIDS: " + i + "=" + this.AllMinGIDs[i]);
            }
            */
            this.AllMinGIDs[numVnodes] = 1 + map.getMaxAllGID(); // Set max cap
            
            /*debug code
            for(int i=0; i<this.AllMinGIDs.length; i++) {
                System.out.println("AllMinGIDS: " + i + "=" + this.AllMinGIDs[i]);
            }
            */
        }
            
        // General case.  Need to build a directory via calls to communication functions
        else {
        System.out.println("Creating DirectoryMap: Map is a general case.");
            int flag = generate();
            //if(flag != 0) throw  new JpetraException("error in Directory(BlockMap map)");
            if(flag != 0) {
                System.out.println("error in Directory(BlockMap)");
                System.exit(1);
            }
        }
    }
    
    public int generate() {
        int i;
        boolean SizeIsConst = map.hasConstantElementSize();
        int minAllGID = map.getMinAllGID();
        int maxAllGID = map.getMaxAllGID();
        // DirectoryMap will have a range of elements from the minimum to the maximum
        // GID of the user map, and an IndexBase of minAllGID from the user map
        int dir_NumGlobalElements = maxAllGID - minAllGID + 1;
        
        // Create a uniform linear map to contain the directory
        DirectoryMap = new Map( dir_NumGlobalElements, minAllGID, map.getComm() );
        
        int dir_NumMyElements = DirectoryMap.getNumMyElements(); // Get NumMyElements
        
        
        
        // Allocate Processor list and Local Index List.  Initialize to -1s.
        
        if (dir_NumMyElements>0) {
            ProcList = new int[ dir_NumMyElements ];
            LocalIndexList = new int[ dir_NumMyElements ];
            if (!SizeIsConst) SizeList = new int[ dir_NumMyElements ];
                // Initialize values to -1 in case the user global element list does
                // fill all IDs from minAllGID to maxAllGID (e.g., allows global indices to be 
                // all even integers.
                for (i=0; i<dir_NumMyElements; i++) {
                    ProcList[i] = -1;
                    LocalIndexList[i] = -1;
                    if (!SizeIsConst) SizeList[i] = -1;
                }
        }
        
        
        // Get list of processors owning the directory entries for the Map GIDs
        
        int MyPID = map.getComm().getVnodeID();
        
        int Map_NumMyElements = map.getNumMyElements();
        int[] send_procs = null;
        if (Map_NumMyElements>0) send_procs = new int[Map_NumMyElements];
        int[] Map_MyGlobalElements = map.getMyGlobalElements();
        
        /*
        debug code
        for(int z=0; z<Map_MyGlobalElements.length; z++) {
	        System.out.println("Map_MyGlobalElements: " + "i=" + Map_MyGlobalElements[z]);
	    }
	    */
	    
        /*assert(Directorymap.RemoteIDList(Map_NumMyElements, Map_MyGlobalElements, 
                         send_procs, 0)==0); */
        
        DirectoryMap.getRemoteIDList(Map_NumMyElements, Map_MyGlobalElements, 
                         send_procs, null);
        
        boolean det_flag = true;
        
        int num_recvs=0;
        
        Distributor distor = map.getComm().createDistributor();
        
        num_recvs = distor.createFromSends( Map_NumMyElements, send_procs, det_flag );

        /*if (Map_NumMyElements>0) delete [] send_procs;*/
        
        int[] export_elements = null;
        int[] import_elements;
        int[] ElementSizeList = null;
        int packetSize = 3; // Assume we will send GIDs, PIDs and LIDs (will increase to 4 if also sending sizes)
        if (!SizeIsConst) packetSize++; // Must send element size info also
        
        if (Map_NumMyElements>0) {
            if (!SizeIsConst) {
                ElementSizeList = map.getElementSizeList();
            }
            
            export_elements = new int[ packetSize * Map_NumMyElements ];
            /*int * ptr = export_elements;*/
            int index = 0;
            for( i = 0; i < Map_NumMyElements; i++ ) {
                export_elements[index++] = Map_MyGlobalElements[i];
                export_elements[index++] = MyPID;
                export_elements[index++] = i;
                if (!SizeIsConst) {
                    export_elements[index++] = ElementSizeList[i];
                }
            }
        }
        
        //if (num_recvs>0) import_elements = new int[ packetSize * num_recvs ];
        //for (i=0; i< packetSize*num_recvs; i++) import_elements[i] = 0;
        
        //distor->do(reinterpret_cast<char *> (export_elements), 
        //      packetSize * sizeof( int ),
        //      reinterpret_cast<char *> (import_elements) );
        
        
        // old in = packetSize * num_recvs
        import_elements = (int[]) distor.distDo(packetSize, num_recvs, export_elements);
        
        //bool MYPID = (map.Comm().MyPID()==0);
        int curr_LID;
        //if (MYPID) cout << "Processor " << map.Comm().MyPID()<< "  num_recvs = "<< num_recvs << endl << flush;
        /*int * ptr = import_elements;*/
        int index = 0;
        for( i = 0; i < num_recvs; i++ ) {
            /*
            debug code
            System.out.println("incoming GID: " + import_elements[index]);
            */
            
            curr_LID = DirectoryMap.getLID(import_elements[index++]); // Convert incoming GID to Directory LID
            //if (MYPID) cout << " Receive ID = " << i << "  GID = " << import_elements[3*i] << "  LID = " << curr_LID << endl << flush;
            
            /*assert(curr_LID !=-1); // Internal error*/
            if (curr_LID == -1) {
                System.err.println("curr_LID is -1, there was an error!");
                System.exit(1);
            }
            
            /*
            debug code
            System.out.println("ProcList: " + ProcList.length);
            System.out.println("curr_LID: " + curr_LID);
            */
            
            ProcList[ curr_LID ] = import_elements[index++];
            LocalIndexList[ curr_LID ] = import_elements[index++];
            if (!SizeIsConst) SizeList[ curr_LID ] = import_elements[index++];        
        }
        
        /*if (import_elements!=0) delete [] import_elements;
        if (export_elements!=0) delete [] export_elements;
        
        delete Distor;*/
        
        return(0);
    }

    public int getDirectoryEntries(int NumEntries, int [] GlobalEntries,
    int [] Procs, int []LocalEntries, int [] EntrySizes) {
       //original code, can probably delete
       /* int ierr = 0;
        /* MPI: int j; 
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
   return(0);*/



    int ierr = 0;
    int j;
    int i;
    int MyPID = map.getComm().getVnodeID();
    int NumProc = map.getComm().getNumVnodes();
    int n_over_p = map.getNumGlobalElements() / NumProc;
    
    System.out.println("Doing getDirectoryEntries");
    // Test for simple cases
    
    // Uniprocessor and local map cases
    
    if (!map.isDistributedGlobal()) {
        System.out.println("GDE: Map is not global.");
        int ElementSize = 0;
        int[] ElementSizeList = null;
        boolean ConstantElementSize = map.hasConstantElementSize();
        if (ConstantElementSize)
          ElementSize = map.getMaxElementSize();
        else
          ElementSizeList = map.getElementSizeList();
        for (i=0; i<NumEntries; i++) {
          int LID = map.getLID(GlobalEntries[i]); // Get LID
          // Procs[i] will be MyPID, or -1 if the GID is not owned by this map
          if (LID==-1) {
        Procs[i] = -1; 
        ierr = 1; // Send warning error back that one of the GIDs is not part of this map
          }
          else Procs[i] = MyPID;
        
          // Put LID in return array if needed
          if (LocalEntries!=null) LocalEntries[i] = LID;
          
          // Fill EntrySizes if needed
          if (EntrySizes!=null) {
        if (ConstantElementSize)
          EntrySizes[i] = ElementSize;
        else if (LID>-1) 
          EntrySizes[i] = ElementSizeList[LID];
        else
          EntrySizes[i] = 0;
          }
        }
        // EPETRA_CHK_ERR(ierr);
        if (ierr != 0 ) {
            System.err.println("Error " + ierr + "in Directory.getEntries");
            System.exit(1);            
        }
        
        return(0);
    }
    
    // Linear Map case
    // Procs are the processors sending to
    if (map.isLinearMap()) {
        System.out.println("GDE: Map is linear");
        int MinAllGID = map.getMinAllGID(); // Get Min of all GID
        int MaxAllGID = map.getMaxAllGID(); // Get Max of all GID
        for (i=0; i<NumEntries; i++) {
              int LID = -1; // Assume not found
              int Proc = -1;
              int GID = GlobalEntries[i];
              if (GID<MinAllGID) ierr = 1;
              else if (GID>MaxAllGID) ierr = 1;
              else {
                // Guess uniform distribution and start a little above it
                int Proc1 = Math.min(GID/Math.max(n_over_p,1) + 2, NumProc-1);
                boolean found = false;
                while (Proc1 >= 0 && Proc1< NumProc) {
                  if (AllMinGIDs[Proc1]<=GID) {
                    if (GID <AllMinGIDs[Proc1+1]) {
                    found = true;
                    break;
                    }
                    else Proc1++;
                  }
                  else Proc1--;
                }
                if (found) {
                  Proc = Proc1;
                  LID = GID - AllMinGIDs[Proc];
                }
              }
              Procs[i] = Proc;
              if (LocalEntries!=null) LocalEntries[i] = LID;
        }
        if (EntrySizes!=null) {
              if (map.hasConstantElementSize()) {
            int ElementSize = map.getMaxElementSize();
            for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
              }
              else {
            int[] ElementSizeList = map.getElementSizeList(); // We know this exists
            
            
            Distributor Size_Distor = map.getComm().createDistributor();
            
            int Size_num_sends = -1;
            int[] Size_send_gids = null;
            int[] Size_send_procs;
            
            
            DistObject out = Size_Distor.createFromReceives( NumEntries, GlobalEntries, Procs, true);
            // out => Size_num_sends, Size_send_gids, Size_send_procs
            Size_num_sends = out.getNumExportIDs();
            Size_send_gids = out.getExportGIDs();
            Size_send_procs = out.getExportPIDs();
            
            int[] Size_exports = null;
            int[] Size_imports = null;
            if (Size_num_sends>0) {
              Size_exports = new int[ 2 * Size_num_sends ];
              for( i = 0; i < Size_num_sends; i++ )
                {
                  int Size_curr_GID = Size_send_gids[i];
                  int Size_curr_LID = map.getLID(Size_curr_GID);
                  // assert(Size_curr_LID!=-1); // Internal error
                  if (Size_curr_LID == -1) {
                        System.err.println("Size_curr_LID == -1 in getDirectoryEnteries Linear Map case.");
                        System.exit(1);
                  }
                  Size_exports[2*i] = Size_curr_GID;
                  int Size_curr_size = ElementSizeList[Size_curr_LID];
                  Size_exports[2*i+1] = Size_curr_size;
                }
            }
            
            /* if (NumEntries>0) Size_imports = new int[ 2 * NumEntries ]; */
            
            Size_imports = Size_Distor.distDo(2, NumEntries, Size_exports);
            
            for( i = 0; i < NumEntries; i++ )
              {
            
                // Need to change !!!!
                //bool found = false;
                int Size_curr_LID = Size_imports[2*i];
                for( j = 0; j < NumEntries; j++ )
                  if( Size_curr_LID == GlobalEntries[j] )
                {
                  EntrySizes[j] = Size_imports[2*i+1];
                  // found = true;
                  break;
                }
                //	if (!found) cout << "Internal error:  Epetra_BasicDirectory::GetDirectoryEntries: Global Index " << curr_LID
                //	     << " not on processor " << MyPID << endl; abort();
              }
            
            /*if( Size_send_gids != 0 ) delete [] Size_send_gids;
            if( Size_send_procs != 0 ) delete [] Size_send_procs;
            
            if( Size_imports != 0 ) delete [] Size_imports;
            if( Size_exports != 0 ) delete [] Size_exports;
            
            delete Size_Distor;*/
              }
        }
        if (ierr != 0 ) {
            System.err.println("Error " + ierr + "in Directory.getEntries");
            System.exit(1);
        }
        return(0);
    }
    
    // General case (need to set up an actual directory structure)
    System.out.println("GDE: General Case, need to setup directory structure.");
    int[] ElementSizeList;
    int PacketSize = 2; // We will send at least the GID and PID.  Might also send LID and Size info
    boolean DoSizes = false;
    if (EntrySizes!=null) {
        if (map.hasConstantElementSize()) {
          int ElementSize = map.getMaxElementSize();
        for (i=0; i<NumEntries; i++) EntrySizes[i] = ElementSize;
        }
        else {
          ElementSizeList = map.getElementSizeList(); // We know this exists
          DoSizes = true;
          PacketSize++; // Sending Size info
        }
    }
    
    boolean DoLIDs = (LocalEntries!=null); // Do LIDs?
    if (DoLIDs) PacketSize++; // Sending LIDs also
    
    
    Distributor Distor = DirectoryMap.getComm().createDistributor();
    
    
    int[] dir_procs = null;
    if (NumEntries>0) dir_procs = new int[ NumEntries ];
    
    // Get directory locations for the requested list of entries
    DirectoryMap.getRemoteIDList(NumEntries, GlobalEntries, dir_procs, null);
    
    //Check for unfound GlobalEntries and set cooresponding Procs to -1
    int NumMissing = 0;
    {for( i = 0; i < NumEntries; ++i )
    if( dir_procs[i] == -1 )
    {
      Procs[i] = -1;
      if (DoLIDs) LocalEntries[i] = -1;
      ++NumMissing;
    }}
    
    int num_sends = -1;
    int[] send_gids = null;
    int[] send_procs;
    
    DistObject out = Distor.createFromReceives( NumEntries, GlobalEntries, dir_procs, true);
    // out => num_sends, send_gids, send_procs
    num_sends = out.getNumExportIDs();
    send_gids = out.getExportGIDs();
    send_procs = out.getExportPIDs();
    
    /* if (NumEntries>0) delete [] dir_procs; */
    
    
    int curr_LID;
    int[] exports = null;
    int[] imports;
    if (num_sends>0) {
        exports = new int[ PacketSize * num_sends ];
        int offSet = 0;
        for( i = 0; i < num_sends; i++ )
          {
        int curr_GID = send_gids[i];
        exports[offSet++] = curr_GID;
        curr_LID = DirectoryMap.getLID(curr_GID);
        // assert(curr_LID!=-1); // Internal error 
        if (curr_LID == -1) {
            System.err.println("curr_LID is -1 while setting up a general case map.");
            System.exit(1);
        }
        exports[offSet++] = ProcList[ curr_LID ];
        if (DoLIDs) exports[offSet++] = LocalIndexList[curr_LID];
        if (DoSizes) exports[offSet++] = SizeList[curr_LID];
          }
    }
    
    int NumRecv = NumEntries - NumMissing;
    
    /* if (NumRecv>0) imports = new int[PacketSize*NumRecv]; */
    
    imports = Distor.distDo(PacketSize, NumRecv, exports);
    
    // int * ptr = imports;
    int offSet = 0;
    for( i = 0; i < NumRecv; i++ ) {
    curr_LID = imports[offSet++];
    for( j = 0; j < NumEntries; j++ )
      if( curr_LID == GlobalEntries[j] ) {
        Procs[j] = imports[offSet++];
        if (DoLIDs) LocalEntries[j] = imports[offSet++];
        if (DoSizes) EntrySizes[j] = imports[offSet++];
        break;
      }
    }
    
    /*if( send_gids ) delete [] send_gids;
    if( send_procs ) delete [] send_procs;
    
    if( imports ) delete [] imports;
    if( exports ) delete [] exports;
    
    delete Distor;*/
    return(0);
    }
}
