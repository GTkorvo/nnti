import CCJ.*;

import java.io.Serializable;

public class testR extends ColMember {
    int gSize;
    ColGroup group;
    
    public static void main(String[] args) {
    
        try {
            ColGroupMaster groupMaster = new ColGroupMaster("ccjhosts.txt");
            new testR(groupMaster);
        } catch (CCJException e) {
            System.err.println("Error in PingPong constructor: " + e);
            e.printStackTrace();
            System.exit(1);
        } 
    }
    
    public testR(ColGroupMaster groupMaster)  throws CCJException {
        super();
        groupMaster.addMember("myGroup", this);

	    gSize = groupMaster.getNumberOfCpus();

	    group = groupMaster.getGroup("myGroup",gSize);
	    
	    begin();
        
    }
    
    public void run() {
        String data;
        String receive = new String("blank: ");
        int rank = -1;
        
        try {
        
        rank = group.getRank(this);
        
        int[] records = new int[]{-1};
        boolean record=false;
        
        if (rank == 0) {
            data = new String("hello number one");
            record=true;
            setupRecords(2);
            barrier(group);
            send_async(group, data, 1);
            //receive = (String) receive(group, 1);
            send_async(group, data, 2);
            //receive = new String(receive + (String) receive(group, 2));
        }
        else if (rank == 1) {
            setupRecords(1);
            record=true;
            barrier(group);
            data = new String("hello root! (from " + rank + ")");
            send_async(group, data, 0);
            //receive = (String) receive(group, 0);
        }
        else {
            setupRecords(1);
            barrier(group);
            record=true;
            data = new String("hello root! (from " + rank + ")");
            send_async(group, data, 0);
            //receive = (String) receive(group, 0);
        }
        
        if (record) {
            records = getRecords();
        }
        
        for(int i = 0; i < records.length; i++) {
            System.out.println("-->#" + records[i]);
        }
        for(int i = 0; i < records.length; i++) {
            if (records[i] != -1) {
                receive = new String(receive + (String) receive(group, records[i]));
            }
        }
        
        endRecords();
        
        } catch (CCJException e) {
	    System.out.println(rank + ": Test : oops, something is wrong !");
	    }
        System.out.println("Recieved: "+ receive);
        System.exit(1);	    
    }

}