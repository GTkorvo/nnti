import Jpetra.*;
import java.util.Arrays;

public class CcjCommBroadcastTest {
    
    public static void main (String[] args) {
        new CcjCommBroadcastTest(args);
     }
     
     public CcjCommBroadcastTest(String[] args) {
        Jpetra.CcjComm comm = new Jpetra.CcjComm("ccjhosts.txt");
        Jpetra.JpetraObject output = new JpetraObject();
        
        if(args.length > 0 && args[0].equals("-v")) output.setRootPrint("VERBOSE", true);
     
        output.println("VERBOSE", "Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
        
        double[] toSend = new double[3];
        if (comm.getVnodeID() == 0) {
            toSend = new double[]{1.5,2.6,3.7};
        }       
        
        
        double[] toRecieve = (double[]) comm.broadcast(toSend,0);
        
        if (!java.util.Arrays.equals(toRecieve, new double[]{1.5,2.6,3.7})) {output.println("ERR", "Test Failed.");}
        
        output.println("VERBOSE", toRecieve[0] + " " + toRecieve[1] + " " + toRecieve[2]);
        
        System.exit(1);
    }
}