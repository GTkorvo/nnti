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

/*
 * Exporter.java
 *
 * Created on June 19, 2001, 9:55 AM
 */

package Jpetra;

/**
 *
 * @author  mwboldt
 * @version 
 */
public class Exporter extends JpetraObject {

    private BlockMap targetMap = null;
    private BlockMap sourceMap = null;
    
    private int numSameIDs = 0;
    private int numPermuteIDs = 0;
    private int [] permuteToLIDs = null;
    private int [] permuteFromLIDs = null;
    private int numRemoteIDs = 0;
    private int [] remoteLIDs = null;
    private int numExporterIDs = 0;
    private int [] exportLIDs = null;
    private int [] exportPIDs = null;
    private int numSend = 0;
    private int numRecv = 0;
    
    
    private Distributor distor;
    // #ifdef PETRA_MPI ...
    
    /** Creates new Exporter */
    public Exporter(BlockMap sourceMap, BlockMap targetMap)
    //throws JpetraException
    {
        this.targetMap = targetMap;
        this.sourceMap = sourceMap;
        
        int i;
        
        int numSourceIDs = sourceMap.getNumMyElements();
        int numTargetIDs = targetMap.getNumMyElements();
        
        int [] targetGIDs = null;
        if(numTargetIDs > 0) {
            targetGIDs = new int [numTargetIDs];
            targetMap.getMyGlobalElements(targetGIDs);
        }
        
        int [] sourceGIDs = null;
        if(numSourceIDs > 0) {
            sourceGIDs = new int [numSourceIDs];
            sourceMap.getMyGlobalElements(sourceGIDs);
        }
        
        int minIDs = Math.min(numSourceIDs, numTargetIDs);
        
        numSameIDs = 0;
        for(i=0; i<minIDs; i++)
            if(targetGIDs[i] == sourceGIDs[i]) numSameIDs++;
            else break;
        
        for(i=numSameIDs; i<numSourceIDs; i++)
            if(targetMap.isMyGID(sourceGIDs[i])) numPermuteIDs++;
            else numExporterIDs++;
        
        int [] exportGIDs = null;
        if(numExporterIDs > 0) {
            exportLIDs = new int [numExporterIDs];
            exportGIDs = new int [numExporterIDs];
        }
        
        if(numPermuteIDs > 0) {
            permuteToLIDs = new int[numPermuteIDs];
            permuteFromLIDs = new int [numPermuteIDs];
        }
        
        numPermuteIDs = 0;
        numExporterIDs = 0;
        for(i=numSameIDs; i<numSourceIDs; i++) {
            if(targetMap.isMyGID(sourceGIDs[i])) {
                permuteFromLIDs[numPermuteIDs] = i;
                permuteToLIDs[numPermuteIDs++] = targetMap.getLID(sourceGIDs[i]);
            }
            else {
                numSend += sourceMap.getElementSize(i);
                exportGIDs[numExporterIDs] = sourceGIDs[i];
                exportLIDs[numExporterIDs++] = i;
            }
        }
        
        //if(numExporterIDs > 0 && !sourceMap.isDistributedGlobal())
        //    throw new JpetraException("serial Exporter has remote IDs");
        if(numExporterIDs > 0 && !sourceMap.isDistributedGlobal()) {
            System.out.println("serial Exporter has remote IDs");
            System.exit(1);
        }
        
        // Test for distributed cases
        int ierr = 0;
      
        if (sourceMap.isDistributedGlobal()) {
      
          if (numExporterIDs>0) exportPIDs = new int[numExporterIDs];
          ierr = targetMap.getRemoteIDList(numExporterIDs, exportGIDs, exportPIDs, null); // Get remote PIDs
          // if( ierr ) throw ReportError("Error in EpetraBlockMap::RemoteIDList", ierr);
          if( ierr != 0) System.err.println("Error in BlockMap.getRemoteIDList " + ierr);
          
          //Get rid of IDs not in Target Map
          if(numExporterIDs>0) {
            int cnt = 0;
            for( i = 0; i < numExporterIDs; ++i )
          if( exportPIDs[i] == -1 ) ++cnt;
            if( cnt != 0 ) {
          int[] newexportGIDs = null;
          int[] newexportPIDs = null;
          int cnt1 = numExporterIDs-cnt;
          if (cnt1 != 0) {
            newexportGIDs = new int[cnt1];
            newexportPIDs = new int[cnt1];
          }
          cnt = 0;
          for( i = 0; i < numExporterIDs; ++i )
            if( exportPIDs[i] != -1 ) {
              newexportGIDs[cnt] = exportGIDs[i];
              newexportPIDs[cnt] = exportPIDs[i];
              ++cnt;
                }
          // assert(cnt==cnt1); // Sanity test
          if (cnt!=cnt1) System.err.println("Failed sanity test cnt==cnt1 in Exporter");
          
          numExporterIDs = cnt;
          /*delete [] exportGIDs;
          delete [] exportPIDs;*/
          exportGIDs = newexportGIDs;
          exportPIDs = newexportPIDs;
          System.err.println("Warning in Exporter: Source IDs not found in Target Map (Do you want to export from subset of Source Map?)");
            }
          }
          
          distor = sourceMap.getComm().createDistributor();
          
          // Construct list of exports that calling processor needs to send as a result
          // of everyone asking for what it needs to receive.
          
          boolean Deterministic = true;
          // ierr = distor->CreateFromSends( numExporterIDs, exportPIDs,Deterministic, NumRemoteIDs);
          numRemoteIDs = distor.createFromSends( numExporterIDs, exportPIDs, Deterministic );
          
          
          // if (ierr!=0) throw ReportError("Error in EpetraDistributor.CreateFromSends()", ierr);
          // if (ierr!=0) System.err.println("Error in EpetraDistributor.CreateFromSends() " + ierr);
          
          // Use comm plan with exportGIDs to find out who is sending to us and
          // get proper ordering of GIDs for remote entries 
          // (that we will convert to LIDs when done).
          
          if (numRemoteIDs>0) remoteLIDs = new int[numRemoteIDs]; // Allocate space for LIDs in target that are
          // going to get something from off-processor.
          // ierr = distor->Do(reinterpretcast<char *> (exportGIDs), 
          //    sizeof( int ),
          //    reinterpretcast<char *> (RemoteLIDs));
          // if (ierr) throw ReportError("Error in EpetraDistributor.Do()", ierr);
          remoteLIDs = distor.distDo(1, numRemoteIDs, exportGIDs);
          
          // Remote IDs come in as GIDs, convert to LIDs
          for (i=0; i< numRemoteIDs; i++) {
            remoteLIDs[i] = targetMap.getLID(remoteLIDs[i]);
            //NumRecv += TargetMap.ElementSize(RemoteLIDs[i]); // Count total number of entries to receive
            numRecv += targetMap.getMaxElementSize(); // Count total number of entries to receive (currently need max)
          }
      
          /*if (numExporterIDs>0) delete [] exportGIDs;*/
        }
        /*if (NumTargetIDs>0) delete [] TargetGIDs;
        if (NumSourceIDs>0) delete [] SourceGIDs;*/
        
        return;
    }
    
    public Exporter(Exporter exporter) {
        targetMap = exporter.targetMap;
        sourceMap = exporter.sourceMap;
        numSameIDs = exporter.numSameIDs;
        numPermuteIDs = exporter.numPermuteIDs;
        numRemoteIDs = exporter.numRemoteIDs;
        numExporterIDs = exporter.numExporterIDs;
        numSend = exporter.numSend;
        numRecv = exporter.numRecv;
        
        int i;
        
        if(numPermuteIDs > 0) {
            permuteToLIDs = new int [numPermuteIDs];
            permuteFromLIDs = new int [numPermuteIDs];
            System.arraycopy(exporter.permuteToLIDs, 0, permuteToLIDs, 0, numPermuteIDs);
            System.arraycopy(exporter.permuteFromLIDs, 0, permuteFromLIDs, 0, numPermuteIDs);
        }
        
        if(numRemoteIDs > 0) {
            remoteLIDs = new int [numRemoteIDs];
            System.arraycopy(exporter.remoteLIDs, 0, remoteLIDs, 0, numRemoteIDs);
        }
        
        targetMap.getComm().barrier();
        if(numExporterIDs > 0) {
            exportLIDs = new int [numExporterIDs];
            exportPIDs = new int [numExporterIDs];
            System.arraycopy(exporter.exportLIDs, 0, exportLIDs, 0, numExporterIDs);
            System.arraycopy(exporter.exportPIDs, 0, exportPIDs, 0, numExporterIDs);
        }
        
        // #ifdef PETRA_MPI ...
    }
    
    public int getNumSameIDs() { return numSameIDs; }
    
    public int getNumPermuteIDs() { return numPermuteIDs; }
    
    public int [] getPermuteFromLIDs() { return permuteFromLIDs; }
    
    public int [] getPermuteToLIDs() { return permuteToLIDs; }
    
    public int getNumRemoteIDs() { return numRemoteIDs; }
    
    public int [] getRemoteLIDs() { return remoteLIDs; }
    
    public int getNumExporterIDs() { return numExporterIDs; }
    
    public int [] getExporterLIDs() { return exportLIDs; }
    
    public int [] getExporterPIDs() { return exportPIDs; }
    
    public int getNumSend() { return numSend; }
    
    public int getNumRecv() { return numRecv; }
    
    public BlockMap getSourceMap() { return sourceMap; }
    
    public BlockMap getTargetMap() { return targetMap; }
    
    // #ifdef PETRA_MPI ...
}
