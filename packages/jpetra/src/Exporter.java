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
    
    // #ifdef PETRA_MPI ...
    
    /** Creates new Exporter */
    public Exporter(BlockMap sourceMap, BlockMap targetMap)
    //throws JpetraException
    {
        this.targetMap = targetMap;
        this.sourceMap = sourceMap;
        
        int i;
        
        int numSourceIDs = sourceMap.getNumVnodeElements();
        int numTargetIDs = targetMap.getNumVnodeElements();
        
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
        
        // #ifdef PETRA_MPI ...
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
