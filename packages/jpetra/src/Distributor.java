package Jpetra;

import java.io.Serializable;

public interface Distributor {
    
    public int createFromSends(int numExportIDs, int[] exportPIDs, boolean deterministic); //, int NumRemoteIDs

    public DistObject createFromReceives(int numRemoteIDs, int[] remoteGIDs, int[] remotePIDs, boolean deterministic); //int NumExportIDs, int[] ExportGIDs, int[] ExportPIDs
    
    // note that do is a reserved word, so must use something else, hence distDo
    public int[] distDo (int importSize, int numRecvs, int[] exportObjects);
}
