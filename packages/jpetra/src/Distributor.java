package Jpetra;

public interface Distributor {
    public DistObject CreateFromSends(int NumExportIDs, int[] ExportPIDs, boolean Deterministic); //, int NumRemoteIDs

    public DistObject CreateFromReceives(int NumRemoteIDs, int RemoteGIDs, int[] RemotePIDs, boolean Deterministic); //int NumExportIDs, int[] ExportGIDs, int[] ExportPIDs

    public Serilizable Do (Serilizable[] exportObjects);
}
