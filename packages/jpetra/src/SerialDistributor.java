package Jpetra;

import java.io.Serializable;

public class SerialDistributor implements Distributor {
   //DistObject output;
    Comm myComm;
    
    public SerialDistributor (Comm myComm) {
    	this.myComm = myComm;
    }
    
    public DistObject CreateFromSends(int NumExportIDs, int[] ExportPIDs, boolean Deterministic) {
        int NumRemoteIDs = 0;
        return new DistObject(NumRemoteIDs);
    }

    public DistObject CreateFromReceives(int NumRemoteIDs, int RemoteGIDs, int[] RemotePIDs, boolean Deterministic) { //int NumExportIDs, int[] ExportGIDs, int[] ExportPIDs
        int NumExportIDs = 0;
        int[] ExportGIDs = null;
        int[] ExportPIDs = null;
        
        return new DistObject(NumExportIDs, ExportGIDs, ExportPIDs);

    }
    
    public Serializable[] Do (Serializable[] exportObjects) {
    	return exportObjects;
    }
}
