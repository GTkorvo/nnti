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

public class Importer extends JpetraObject {

    BlockMap targetMap;
    BlockMap sourceMap;
    
    int  numSameIDs;
    int  numPermuteIDs;
    int[] permuteToLIDs;
    int[] permuteFromLIDs;
    int  numRemoteIDs;
    int[] remoteLIDs;
    
    int  numExportIDs;
    int[] exportLIDs;
    int[] exportPIDs;
    
    int numSend;
    int numRecv;
    
    Distributor Distor;
    
    public Importer(BlockMap targetMap, BlockMap sourceMap) {
        this.targetMap = targetMap;
        this.sourceMap = sourceMap;
        this.numSameIDs = 0;
        this.numPermuteIDs = 0;
        this.permuteToLIDs = null;
        this.permuteFromLIDs = null;
        this.numRemoteIDs = 0;
        this.remoteLIDs = null;
        this.numExportIDs = 0;
        this.exportLIDs = null;
        this.exportPIDs = null;
        this.numSend = 0;
        this.numRecv = 0;
        this.Distor = null;
        
      int i;
      
      // Build three ID lists:
      // numSameIDs - Number of IDs in targetMap and sourceMap that are identical, up to the first
      //              nonidentical ID.
      // numPermuteIDs - Number of IDs in sourceMap that must be indirectly loaded but are on this processor.
      // numRemoteIDs - Number of IDs that are in sourceMap but not in targetMap, and thus must be imported.
      
      int NumSourceIDs = sourceMap.getNumMyElements();
      int NumTargetIDs = targetMap.getNumMyElements();
      
      int[] TargetGIDs = null;
      if (NumTargetIDs>0) {
        TargetGIDs = new int[NumTargetIDs];
        targetMap.getMyGlobalElements(TargetGIDs);
      }
      
      int[] SourceGIDs = null;
      if (NumSourceIDs>0) {
        SourceGIDs = new int[NumSourceIDs];
        sourceMap.getMyGlobalElements(SourceGIDs);
      }
      
      int MinIDs = Math.min(NumSourceIDs, NumTargetIDs);
      
      
      numSameIDs = 0;
      for (i=0; i< MinIDs; i++) if (TargetGIDs[i]==SourceGIDs[i]) numSameIDs++; else break;
      
      
      // Find count of Target IDs that are truly remote and those that are local but permuted
    
      numPermuteIDs = 0;
      numRemoteIDs = 0;
      for (i=numSameIDs; i< NumTargetIDs; i++) 
        if (sourceMap.isMyGID(TargetGIDs[i])) numPermuteIDs++; // Check if Target GID is a local Source GID
        else numRemoteIDs++; // If not, then it is remote
      
      
      
      // Define remote and permutation lists
      
      int[] RemoteGIDs=null;
      remoteLIDs = null;
      if (numRemoteIDs>0) {
        remoteLIDs = new int[numRemoteIDs];
        RemoteGIDs = new int[numRemoteIDs];
      }
      if (numPermuteIDs>0)  {
        permuteToLIDs = new int[numPermuteIDs];
        permuteFromLIDs = new int[numPermuteIDs];
      }
      
      numPermuteIDs = 0;
      numRemoteIDs = 0;
      for (i=numSameIDs; i< NumTargetIDs; i++) {
        if (sourceMap.isMyGID(TargetGIDs[i])) {
          permuteToLIDs[numPermuteIDs] = i;
          permuteFromLIDs[numPermuteIDs++] = sourceMap.getLID(TargetGIDs[i]);
        }
        else {
          //numRecv +=targetMap.ElementSize(i); // Count total number of entries to receive
          numRecv +=targetMap.getMaxElementSize(); // Count total number of entries to receive (currently need max)
          RemoteGIDs[numRemoteIDs] = TargetGIDs[i];
          remoteLIDs[numRemoteIDs++] = i;
        }
      }
    
      if( numRemoteIDs>0 && !sourceMap.isDistributedGlobal() )
        System.err.println("Warning in Epetra_Import: Serial Import has remote IDs. (Importing to Subset of Target Map)");
      
      // Test for distributed cases
      
      int[] RemotePIDs = null;
    
      if (sourceMap.isDistributedGlobal()) {
        
        if (numRemoteIDs>0)  RemotePIDs = new int[numRemoteIDs];
        int ierr = sourceMap.getRemoteIDList(numRemoteIDs, RemoteGIDs, RemotePIDs, null); // Get remote PIDs
        // if (ierr) throw ReportError("Error in sourceMap.RemoteIDList call", ierr);
        if (ierr != 0) {
            System.err.println("Error in sourceMap.RemoteIDList call: " + ierr);
        }
        
        //Get rid of IDs that don't exist in sourceMap
        if(numRemoteIDs>0) {
          int cnt = 0;
          for( i = 0; i < numRemoteIDs; ++i )
            if( RemotePIDs[i] == -1 ) ++cnt;
          if( cnt != 0 ) {
            int[] NewRemoteGIDs = new int[numRemoteIDs-cnt];
            int[] NewRemotePIDs = new int[numRemoteIDs-cnt];
            cnt = 0;
            for( i = 0; i < numRemoteIDs; ++i )
              if( RemotePIDs[i] != -1 ) {
                NewRemoteGIDs[cnt] = RemoteGIDs[i];
                NewRemotePIDs[cnt] = RemotePIDs[i];
                ++cnt;
              }
            numRemoteIDs = cnt;
            /*delete [] RemoteGIDs;
            delete [] RemotePIDs;*/
            RemoteGIDs = NewRemoteGIDs;
            RemotePIDs = NewRemotePIDs;
            System.err.println("Warning in Epetra_Import: Serial Import has remote IDs. (Importing to Subset of Target Map)");
          }
        }
    
        Distor = sourceMap.getComm().createDistributor();
        
        // Construct list of exports that calling processor needs to send as a result
        // of everyone asking for what it needs to receive.
        
        boolean Deterministic = true;
        DistObject out = Distor.createFromReceives( numRemoteIDs, RemoteGIDs, RemotePIDs, Deterministic);
        numExportIDs = out.getNumExportIDs();
        exportLIDs = out.getExportGIDs();
        exportPIDs = out.getExportPIDs();
        
        
        /* ierr = Distor->CreateFromRecvs( numRemoteIDs, RemoteGIDs,
                         RemotePIDs, Deterministic,
                         numExportIDs, exportLIDs, exportPIDs );
        
        if (ierr!=0) throw ReportError("Error in Epetra_Distributor.CreateFromRecvs()", ierr);*/
        
        // Use comm plan with Export GIDs (stored in exportLIDs) to
        // get proper ordering of GIDs for remote entries 
        // (that we will convert to LIDs when done).
        
        RemoteGIDs = Distor.distDo(1, RemoteGIDs.length, exportLIDs);
        
        /*ierr = Distor->Do( reinterpret_cast<char *> (exportLIDs), 
            sizeof( int ),
            reinterpret_cast<char *> (RemoteGIDs));
        
        if (ierr!=0) throw ReportError("Error in Epetra_Distributor.Do()", ierr);*/
    
            
        // Export IDs come in as GIDs, convert to LIDs
        for (i=0; i< numExportIDs; i++) {
          // if (exportPIDs[i] < 0) throw ReportError("targetMap requested a GID that is not in the sourceMap.", -1);
          if (exportPIDs[i] < 0) System.err.println("targetMap requested a GID that is not in the sourceMap.");
          
          exportLIDs[i] = sourceMap.getLID(exportLIDs[i]);
          //numSend += sourceMap.ElementSize(exportLIDs[i]); // Count total number of entries to send
          numSend += sourceMap.getMaxElementSize(); // Count total number of entries to send (currently need max)
        }
        
        // Remote IDs come in as GIDs, convert to LIDs in proper order
    
        // for (i=0; i< numRemoteIDs; i++) remoteLIDs[i] = targetMap.getLID(RemoteGIDs[i]); // Only works when target map has no repeated GIDs
        
        if (numRemoteIDs>0) {
          int[] ReorderedremoteLIDs = RemotePIDs; // Reuse some temp space
          for (i=0; i< numRemoteIDs; i++) {
        int CurrentGID = RemoteGIDs[i];
        boolean Found = false;
        for (int j=0; j < numRemoteIDs; j++) {
          if (remoteLIDs[j]!= -1) {
            if (CurrentGID==TargetGIDs[remoteLIDs[j]]) {
              ReorderedremoteLIDs[i] = remoteLIDs[j];
              remoteLIDs[j] = -1;
              Found = true;
              break;
            }
          }
        }
        // if (!Found) throw ReportError("Internal error.  Cannot map incoming GID to Target Map", -2);
        if (!Found) System.err.println("Internal error.  Cannot map incoming GID to Target Map");
          }
          
          // Clean up and leave....
          /*delete [] remoteLIDs;*/
          remoteLIDs = ReorderedremoteLIDs;
        }
      }
    
      /*if( RemoteGIDs ) delete [] RemoteGIDs;
    
      if (NumTargetIDs>0) delete [] TargetGIDs;
      if (NumSourceIDs>0) delete [] SourceGIDs;*/
      
      return;
    }
    
    public Importer(Importer importer) {
        this.targetMap = importer.targetMap;
        this.sourceMap = importer.sourceMap;
        this.numSameIDs = importer.numSameIDs;
        this.numPermuteIDs = importer.numPermuteIDs;
        this.numRemoteIDs = importer.numRemoteIDs;
        this.numExportIDs = importer.numExportIDs;
        this.numSend = importer.numSend;
        this.numRecv = importer.numRecv;
        
        int i;
        
        if(numPermuteIDs > 0) {
            permuteToLIDs = new int [numPermuteIDs];
            permuteFromLIDs = new int [numPermuteIDs];
            System.arraycopy(importer.permuteToLIDs, 0, permuteToLIDs, 0, numPermuteIDs);
            System.arraycopy(importer.permuteFromLIDs, 0, permuteFromLIDs, 0, numPermuteIDs);
        }
        
        if(numRemoteIDs > 0) {
            remoteLIDs = new int [numRemoteIDs];
            System.arraycopy(importer.remoteLIDs, 0, remoteLIDs, 0, numRemoteIDs);
        }
        
        if(numExportIDs > 0) {
            exportLIDs = new int [numExportIDs];
            exportPIDs = new int [numExportIDs];
            System.arraycopy(importer.exportLIDs, 0, exportLIDs, 0, numExportIDs);
            System.arraycopy(importer.exportPIDs, 0, exportPIDs, 0, numExportIDs);
        }
        
        
        // what to do about distributor?
        // if (Importer.Distor_!=0) Distor_ = Importer.Distor_->Clone();
        this.Distor = importer.Distor;
        
        /*
        #ifdef PETRA_MPI ...
        */
    }
    
    public int getNumSameIDs() { return numSameIDs; }
    
    public int getNumPermuteIDs() { return numPermuteIDs; }
    
    public int [] getPermuteFromLIDs() { return permuteFromLIDs; }
    
    public int [] getPermuteToLIDs() { return permuteToLIDs; }
    
    public int getNumRemoteIDs() { return numRemoteIDs; }
    
    public int [] getRemoteLIDs() { return remoteLIDs; }
    
    public int getNumExporterIDs() { return numExportIDs; }
    
    public int [] getExporterLIDs() { return exportLIDs; }
    
    public int [] getExporterPIDs() { return exportPIDs; }
    
    public int getNumSend() { return numSend; }
    
    public int getNumRecv() { return numRecv; }
    
    public BlockMap getSourceMap() { return sourceMap; }
    
    public BlockMap getTargetMap() { return targetMap; }
}
