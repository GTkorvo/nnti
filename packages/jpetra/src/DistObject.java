package Jpetra;

public class DistObject extends JpetraObject {
   // int err;
    
    int numRemoteIDs;
    
    int numExportIDs;
    int[] exportGIDs;
    int[] exportPIDs;
    
    public DistObject (int numRemoteIDs) {
        //this.err  = err;
        this.numRemoteIDs = numRemoteIDs;   
    }

    public DistObject (int numExportIDs, int[] exportGIDs, int[] exportPIDs) {
        //this.err  = err;
        this.numExportIDs = numExportIDs;
        this.exportGIDs = exportGIDs;
        this.exportPIDs = exportPIDs;
    }
    
   // public int getErr () {
   //     return this.err;
   // }
    
    public int getNumExportIDs() {
        return numExportIDs;
    }
    
        public int[] getExportGIDs() {
        return exportGIDs;
    }
    
        public int[] getExportPIDs() {
        return exportPIDs;
    }
    
        public int getNumRemoteIDs() {
        return numRemoteIDs;
    }
    
}