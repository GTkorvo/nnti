/*
 * Importer.java
 *
 * Created on June 19, 2001, 9:19 AM
 */

package Jpetra;

/**
 *
 * @author  mwboldt
 * @version 
 */
public class Importer extends JpetraObject {

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
    
    /*
    #ifdef PETRA_MPI
        GSComm_Plan * GSPlan;
    #endif
    */
    
    /** Creates new Importer */
    public Importer(BlockMap targetMap, BlockMap sourceMap) 
	// throws JpetraException
    {
        this.targetMap = targetMap;
        this.sourceMap = sourceMap;
        
        int i;
        int numSourceIDs = sourceMap.getNumProcessElements();
        int numTargetIDs = targetMap.getNumProcessElements();
        
        int [] targetGIDs = null;
        if(numTargetIDs > 0) {
            targetGIDs = new int [numTargetIDs];
            targetMap.getGlobalElements(targetGIDs);
        }
        
        int [] sourceGIDs = null;
        if(numSourceIDs > 0) {
            sourceGIDs = new int [numSourceIDs];
            sourceMap.getGlobalElements(sourceGIDs);
        }
        
        int minIDs = Math.min(numSourceIDs, numTargetIDs);
        
        for(i=0; i<minIDs; i++) 
            if (targetGIDs[i] == sourceGIDs[i]) numSameIDs++;
            else break;
        
        // Find count of target IDs that are truely remote and thos that are local but permuted
        for(i=numSameIDs; i<numTargetIDs; i++)
            if(sourceMap.isMyGID(targetGIDs[i])) numPermuteIDs++;
            else numRemoteIDs++;
        
        // Define remote and permutation lists
        int [] remoteGIDs = null;
        if(numRemoteIDs > 0) {
            remoteLIDs = new int[numRemoteIDs];
            remoteGIDs = new int [numRemoteIDs];
        }
        
        if(numPermuteIDs > 0) {
            permuteToLIDs = new int [numPermuteIDs];
            permuteFromLIDs = new int [numPermuteIDs];
        }
        
        numPermuteIDs = 0;
        numRemoteIDs = 0;
        for(i=numSameIDs; i<numTargetIDs; i++) {
            if(sourceMap.isMyGID(targetGIDs[i])) {
                permuteToLIDs[numPermuteIDs] = i;
                permuteFromLIDs[numPermuteIDs++] = sourceMap.getLID(targetGIDs[i]);
            }
            else {
                numRecv += targetMap.getElementSize(i);
                remoteGIDs[numRemoteIDs] = targetGIDs[i];
                remoteLIDs[numRemoteIDs++] = i;
            }
        }
        
        //if(numRemoteIDs > 0 && !sourceMap.isDistributedGlobal())
        //    throw new JpetraException("serial Importer has remote IDs");

	if(numRemoteIDs > 0 && !sourceMap.isDistributedGlobal()) {
            System.out.println("serial Importer has remote IDs");
	    System.exit(1);
	}
        
        /*
        #ifdef PETRA_MPI ...
        */
    }
    
    public Importer(Importer importer) {
        targetMap = importer.targetMap;
        sourceMap = importer.sourceMap;
        numSameIDs = importer.numSameIDs;
        numPermuteIDs = importer.numPermuteIDs;
        numRemoteIDs = importer.numRemoteIDs;
        numExporterIDs = importer.numExporterIDs;
        numSend = importer.numSend;
        numRecv = importer.numRecv;
        
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
        
        if(numExporterIDs > 0) {
            exportLIDs = new int [numExporterIDs];
            exportPIDs = new int [numExporterIDs];
            System.arraycopy(importer.exportLIDs, 0, exportLIDs, 0, numExporterIDs);
            System.arraycopy(importer.exportPIDs, 0, exportPIDs, 0, numExporterIDs);
        }
        
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
    
    public int getNumExporterIDs() { return numExporterIDs; }
    
    public int [] getExporterLIDs() { return exportLIDs; }
    
    public int [] getExporterPIDs() { return exportPIDs; }
    
    public int getNumSend() { return numSend; }
    
    public int getNumRecv() { return numRecv; }
    
    public BlockMap getSourceMap() { return sourceMap; }
    
    public BlockMap getTargetMap() { return targetMap; }
    
    // #ifdef PETRA_MPI ...
}
