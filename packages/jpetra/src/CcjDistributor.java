package Jpetra;

import java.io.Serializable;
import CCJ.CCJException;

public class CcjDistributor implements Distributor {
    CcjComm myComm;
    CcjLink myCcjLink;
    CCJ.ColGroup group;
    int self_msg;
    int nsends;
    int nrecvs;
    int[] indices_to;
    int[] lengths_to;
    int[] indices_from;
    int[] lengths_from;
    int[] procs_to;
    int[] procs_from;
    int[] recv_array;
    int max_send_length;
    int total_recv_length;
    int size_indices_to;
    int size_indices_from;
    
    int[] tmpnrecvs;
    
    public CcjDistributor (CcjComm myComm) {
    	this.myComm = myComm;
    	this.myCcjLink = myComm.getCcjLink();
    	this.group = myCcjLink.getGroup();
    	System.out.println("Distributor created.");
    }
    
    public int createFromSends(int NumExportIDs, int[] ExportPIDs, boolean Deterministic) {
        //int NumRemoteIDs = 0;
        //return new DistObject(numRemoteIDs);
        System.out.println("Starting createFromSends...");
        int i;
        
        int my_proc = myComm.getVnodeID();
        int nprocs = myComm.getNumVnodes();
                
        // Check to see if items are grouped by processor w/o gaps
        // If so, indices_to -> 0
        
        // Setup data structures for quick traversal of arrays
        int[] starts = new int[ nprocs + 1 ];
        
        /*for( i = 0; i < nprocs; i++ )
        starts[i] = 0;*/
        
        int nactive = 0;
        
        for( i = 0; i < NumExportIDs; i++ ) {
            if( ExportPIDs[i] >= 0 ) {
              ++starts[ ExportPIDs[i] ];
              ++nactive;
            }
        }
        
        self_msg=0;
        if ( starts[my_proc] != 0 ) {
            self_msg = 1;
        }
        
        nsends = 0;
        if( starts[0] != 0 ) nsends = 1;
        
        for( i = 1; i < nprocs; i++ )
        {
        if( starts[i] != 0 ) ++nsends;
        starts[i] += starts[i-1];
        }
        
        for( i = nprocs-1; i != 0; i-- )
        starts[i] = starts[i-1];
        
        starts[0] = 0;
        
        if (nactive>0) {
        indices_to = new int[ nactive ];
        size_indices_to = nactive;
        }
        
        for( i = 0; i < NumExportIDs; i++ )
        if( ExportPIDs[i] >= 0 )
        {
          indices_to[ starts[ ExportPIDs[i] ] ] = i;
          ++starts[ ExportPIDs[i] ];
        }
        
        //Reconstuct starts array to index into indices_to.
        
        for( i = nprocs-1; i != 0; i-- )
        starts[i] = starts[i-1];
        starts[0] = 0;
        starts[nprocs] = nactive;
        
        if (nsends>0) {
        lengths_to = new int[ nsends ];
        procs_to = new int[ nsends ];
        }
        
        int j = 0;
        max_send_length = 0;
        
        for( i = 0; i < nprocs; i++ )
        if( starts[i+1] != starts[i] )
        {
          lengths_to[j] = starts[i+1] - starts[i];
          if( ( i != my_proc ) && ( lengths_to[j] > max_send_length ) )
            max_send_length = lengths_to[j];
          procs_to[j] = i;
          j++;
        }
        
        /*delete [] starts;*/
        
        //Invert map to see what msgs are received and what length
        computeRecvs( my_proc, nprocs, Deterministic );
        
        total_recv_length = 0;
        for( i = 0; i < nrecvs; i++ )
        total_recv_length += lengths_from[i];
        
        /*
        if (nrecvs>0) {
        request_ = new MPI_Request[ nrecvs ];
        status_ = new MPI_Status[ nrecvs ];
        }*/
        
        // NumRemoteIDs = total_recv_length;
        
        System.out.println("createFromSends done.");
        return total_recv_length;
    }

public boolean computeRecvs( int my_proc, int nprocs, boolean Deterministic ) {
        int[][] msg_count = new int[ nprocs ][];
        //int[] counts = new int[ nprocs ];
        
        int i;
        
        System.out.println("Starting ComputeRecvs...");
        
        /*MPI_Status status;*/
        
        /*for( i = 0; i < nprocs; i++ )
        {
        msg_count[i] = 0;
        counts[i] = 1;
        }*/
        
        for( i = 0; i < nsends; i++ ) {
            //msg_count[ procs_to[i] ] = 1;
            msg_count[ procs_to[i] ] = new int[]{my_proc};
            System.out.println("Going to send a message to " + procs_to[i]);
        }
        
        
        
        /*MPI_Reduce_scatter( msg_count, &nrecvs_, counts, MPI_INT, MPI_SUM, comm_ );*/
        
        CcjDistributorReduce idCombiner = new CcjDistributorReduce();
        
        try {
            msg_count = (int[][]) myCcjLink.reduce(group, msg_count, idCombiner, 0);
        } 
        catch (CCJException e) {
            System.err.println("Error in ComputeRecvs:CCJ reduce: " + e);
        }
        
        CcjScatter2DArray scatterData = new CcjScatter2DArray(msg_count);
        
        //tmpnrecvs = new int[1];
        try {
            tmpnrecvs = (int[]) myCcjLink.scatter(group, scatterData, 0);
        } catch (CCJException e) {
            System.err.println("Error in ComputeRecvs:CCJ scatter: " + e);
        }
        
        System.out.println("reduce and scatter done");
        nrecvs=0;
        if (tmpnrecvs != null ) {nrecvs = tmpnrecvs.length;}
        System.out.println("Total number of recvs: " + nrecvs);
        for(int z=0; z<nrecvs; z++) {
            System.out.println("going to receive from "+ tmpnrecvs[z]);
        }
        
        /*delete [] msg_count;
        delete [] counts;*/
        
        if (nrecvs>0) {
            lengths_from = new int[ nrecvs + self_msg ];
            procs_from = new int[ nrecvs + self_msg ];
            
            /*for(i=0; i<nrecvs+self_msg; ++i) {
              lengths_from[i] = 0;
              procs_from[i] = 0;
            }*/
            
        }

        System.out.println("Doing sends:");
        for( i = 0; i < nsends; i++ )
        if( procs_to[i] != my_proc ) {
          /*MPI_Send( &(lengths_to_[i]), 1, MPI_INT, procs_to_[i], tag_, comm_ );*/
            
            try {
                System.out.println("sending to " + procs_to[i] + " lengths_to[i]: " + lengths_to[i]);
                myCcjLink.send_async(group, new int[]{lengths_to[i]}, procs_to[i]);
            } catch (CCJException e) {
                System.err.println("Error in ComputeRecvs:CCJ send_async: " + e);
            }
        }
        else
        {
            System.out.println("Found myself, sending to myself...");
          /*assert(nrecvs_>0);*/
          lengths_from[nrecvs-1] = lengths_to[i];
          procs_from[nrecvs-1] = my_proc;
        }
        System.out.println("Sending done.  Starting to receive...");
        int[] tmprec;
        for( i = 0; i < nrecvs; i++ )
        {
            /*MPI_Recv( &(lengths_from_[i]), 1, MPI_INT, MPI_ANY_SOURCE, tag_, comm_, &status );
            procs_from_[i] = status.MPI_SOURCE;*/
            if (tmpnrecvs[i] != my_proc) {
                try {
                    System.out.println("receiving from: " + tmpnrecvs[i]);
                    tmprec = (int[]) myCcjLink.receive(group, tmpnrecvs[i]);
                    lengths_from[i] = tmprec[0];
                } catch (CCJException e) {
                    System.err.println("Error in ComputeRecvs:CCJ receive: " + e);
                }
            }
        }
        
        System.out.println("Receives done, starting barrier...");
        
        myComm.barrier();
        /*MPI_Barrier( comm_ );*/
        System.out.println("barrier done.");
        
        if( Deterministic )
        {
        int j;
        int temp;
        
        for( i = 1; i < ( nrecvs - self_msg ); i++ )
        {
          j = i;
        
          while( ( j > 0 ) && ( procs_from[j] < procs_from[j-1] ) )
          {
            temp = procs_from[j];
            procs_from[j] = procs_from[j-1];
            procs_from[j-1] = temp;
        
            temp = lengths_from[j];
            lengths_from[j] = lengths_from[j-1];
            lengths_from[j-1] = temp;
        
            j--;
          }
        }
        }
        
        // Compute indices_from_
        
        size_indices_from = 0;
        for( i = 0; i < nrecvs; i++ )  size_indices_from += lengths_from[i];
        indices_from = new int[ size_indices_from ];
        
        for (i=0; i<size_indices_from; i++) indices_from[i] = i;
        
        myComm.barrier();
        
        /*MPI_Barrier( comm_ );*/
        
        
        return false;
}

    public DistObject createFromReceives(int NumRemoteIDs, int[] RemoteGIDs, int[] RemotePIDs, boolean Deterministic) { //int numExportIDs, int[] exportGIDs, int[] exportPIDs
        /*int numExportIDs = 0;
        int[] exportGIDs = null;
        int[] exportPIDs = null;
        
        return new DistObject(numExportIDs, exportGIDs, exportPIDs);*/
        
        int i;

        int my_proc = myComm.getVnodeID();
        int nprocs = myComm.getNumVnodes();;
        
        //computeSends( NumRemoteIDs, RemoteGIDs, RemotePIDs, my_proc);
        // out => NumExportIDs, ExportGIDs, ExportPIDs
        // BEGIN computeSends
             /*
         out NumExportIDs == num_exports
         out ExportGIDs == export_ids
         out ExportPIDs == export_procs
         */
         
        CcjDistributor tmp_plan = new CcjDistributor(myComm);
        //setup input
        int num_imports = NumRemoteIDs;
        int[] import_ids = RemoteGIDs;
        int[] import_procs = RemotePIDs;
        //end
        
        // define new
        int num_exports = 0;
        int[] export_ids = null;
        int[] export_procs = null;
        // end
        
        int[] proc_list = null;
        int[] import_objs = null;
        int[] export_objs = null;
        
        if( num_imports > 0 ) {
            proc_list = new int[ num_imports ];
            import_objs = new int[ 2 * num_imports ];
            
            for( i = 0; i < num_imports; i++ ) {
              proc_list[i] = import_procs[i];
            
              import_objs[2*i] = import_ids[i];
              import_objs[2*i+1] = my_proc;
            }
        }
        
        num_exports = tmp_plan.createFromSends( num_imports, proc_list, true); 
        if( num_exports > 0 ) {
            export_objs = new int[ 2 * num_exports ];
            export_ids = new int[ num_exports ]; // Note: export_ids and export_procs must be deleted by the calling routine
            export_procs = new int[ num_exports ];
        }
        else {
            export_ids = null;
            export_procs = null;
        }
        
        export_objs = tmp_plan.distDo(2, NumRemoteIDs, import_objs);
        
        
        for( i = 0; i < num_exports; i++ ) {
            export_ids[i] = export_objs[2*i];
            export_procs[i] = export_objs[2*i+1];
        }
        
        
        // END computeSends
        // setup output
        int NumExportIDs = num_exports;
        int[] ExportGIDs = export_ids;
        int[] ExportPIDs = export_procs;
        // end
        
         // Setup data structures for quick traversal of arrays
        int[] starts = new int[ nprocs + 1 ];
        
        /*for( i = 0; i < nprocs; i++ )
        starts[i] = 0;*/
        
        int nactive = 0;
        
        for( i = 0; i < NumExportIDs; i++ ) {
            if( ExportPIDs[i] >= 0 ) {
              ++starts[ ExportPIDs[i] ];
              ++nactive;
            }
        }
        
        self_msg=0;
        if ( starts[my_proc] != 0 ) {
            self_msg = 1;
        }
        
        nsends = 0;
        if( starts[0] != 0 ) nsends = 1;
        
        for( i = 1; i < nprocs; i++ )
        {
        if( starts[i] != 0 ) ++nsends;
        starts[i] += starts[i-1];
        }
        
        for( i = nprocs-1; i != 0; i-- )
        starts[i] = starts[i-1];
        
        starts[0] = 0;
        
        if (nactive>0) {
        indices_to = new int[ nactive ];
        size_indices_to = nactive;
        }
        
        for( i = 0; i < NumExportIDs; i++ )
        if( ExportPIDs[i] >= 0 )
        {
          indices_to[ starts[ ExportPIDs[i] ] ] = i;
          ++starts[ ExportPIDs[i] ];
        }
        
        //Reconstuct starts array to index into indices_to.
        
        for( i = nprocs-1; i != 0; i-- )
        starts[i] = starts[i-1];
        starts[0] = 0;
        starts[nprocs] = nactive;
        
        if (nsends>0) {
        lengths_to = new int[ nsends ];
        procs_to = new int[ nsends ];
        }
        
        int j = 0;
        max_send_length = 0;
        
        for( i = 0; i < nprocs; i++ )
        if( starts[i+1] != starts[i] )
        {
          lengths_to[j] = starts[i+1] - starts[i];
          if( ( i != my_proc ) && ( lengths_to[j] > max_send_length ) )
            max_send_length = lengths_to[j];
          procs_to[j] = i;
          j++;
        }
        
        /*delete [] starts;*/
        
        //Invert map to see what msgs are received and what length
        computeRecvs( my_proc, nprocs, Deterministic );
        
        total_recv_length = 0;
        for( i = 0; i < nrecvs; i++ )
        total_recv_length += lengths_from[i];
        
        /*
        if (nrecvs>0) {
        request_ = new MPI_Request[ nrecvs ];
        status_ = new MPI_Status[ nrecvs ];
        }*/
        
        NumRemoteIDs = total_recv_length;
        
        System.out.println("createFromReceives done.");
        return new DistObject(NumExportIDs, ExportGIDs, ExportPIDs);
    }
    
    public int[] distDo (int importSize, int numRecvs, int[] exportObjects) {
      int i, j, k;
    
    int obj_size = importSize;
    int[] export_objs = exportObjects;
    int self_recv_address = 0;
    
    int my_proc = myComm.getVnodeID();
        
    //recv_array_ = import_objs;
    int[] recv_array = new int[numRecvs * obj_size];
    int[] send_array;
    j = 0;
    k = 0;
    
    int self_num = 0;
    int self_index = 0;
    
    //if (max_send_length>0) 
    // send_array = new int[ max_send_length ];
    
    for( i = 0; i < nsends; i++ )
    {
    if( procs_to[i] != my_proc )
    {
      int offset = 0;
      send_array = new int[lengths_to[i]*obj_size];
      for( k = 0; k < lengths_to[i]; k++ )
      {
        //memcpy( &(send_array_[offset]), 
      //&(export_objs[indices_to_[j]*obj_size]), obj_size );
        
        System.arraycopy(export_objs, indices_to[j]*obj_size, send_array, offset, obj_size);
        send_array[k] = export_objs[indices_to[j]];
        
        j++;
        offset += obj_size;
      }
      //   cout << "my_proc = " << my_proc << " length = " << lengths_to_[i] * obj_size 
      //   << " send to = " << procs_to_[i] << " tag = " << tag << endl;
      //cout << "Processor " << my_proc << "lengths_to_["<<i<<"] = " << lengths_to_[i] << " obj_size = " <<obj_size<<endl;
    /*  MPI_Rsend( send_array_, lengths_to_[i] * obj_size,
    MPI_CHAR, procs_to_[i], tag_, comm_ );*/
    
        try {
            myCcjLink.send_async(group, send_array, procs_to[i]);
	    /*
	    debug code
	    System.out.println("Sending to " + procs_to[i]);
	    for(int z=0; z<send_array.length; z++) {
	        System.out.println("send_array[" + z + "]=" + send_array[z]);
	    }
	    */
        } catch (CCJException e) {
            System.err.println("Error in distDo:CCJ send_async: " + e);
        }
    }
    else
    {
      self_num = i;
      self_index = j;
      j += lengths_to[i];
    }
    }
    
    myComm.barrier();
    
    j=0;
    for( i = 0; i < nrecvs; i++ )
    {
        if( procs_from[i] != my_proc )
        {
          //cout << "Processor " << my_proc << "lengths_from_["<<i<<"] = " << lengths_from_[i] << " obj_size = " <<obj_size<<endl;
            /*MPI_Irecv( &(recv_array_[j]), lengths_from_[i] * obj_size,
            MPI_CHAR, procs_from_[i], tag_, comm_,
            &(request_[k]) );*/
            int[] tmp;
            try {
                tmp = (int[]) myCcjLink.receive(group, tmpnrecvs[i]);
                System.arraycopy(tmp, 0, recv_array, j, tmp.length);
                //j += tmp.length;
            } catch (CCJException e) {
                System.err.println("Error in distDo:CCJ receive: " + e);
            }
            // k++;
        }
        else
          self_recv_address = j;

    j += lengths_from[i] * obj_size;
    }
    
    
    

    
    if( self_msg == 1) {
        for( k = 0; k < lengths_to[self_num]; k++ ) {
            /* memcpy( &(recv_array_[self_recv_address]),
            &(export_objs[indices_to_[self_index]*obj_size]),
             obj_size );*/
            
            /*
            debug code
             System.out.println("------------------------------");
             System.out.println("self_recv_address: " + self_recv_address);
             System.out.println("recv_array.length: " + recv_array.length);
             System.out.println("self_index: " + self_index);
             System.out.println("indices_to.length: " + indices_to.length);
             System.out.println("export_objs.length: " + export_objs.length);
             System.out.println("indices_to[self_index]: " + indices_to[self_index]);
             
             System.out.println("lengths_to: " + lengths_to[self_num]);
             */
             
             System.arraycopy(export_objs, indices_to[self_index]*obj_size, recv_array, self_recv_address, obj_size);
             self_index++;
             self_recv_address += obj_size;
        }
    }
    
    /*
    debug code
    System.out.println ("recv_array dump\n==========================");
    for(int z=0; z<recv_array.length; z++) {
        System.out.println("recv_array[" + z + "]=" + recv_array[z]);
    }
    */
    
    System.out.println("Made it through distDo.");
    
    /*// remove, fix, change?
    if (send_array!=null) {
    //delete [] send_array_;
    send_array = null;
    }*/
    
    return recv_array;
    }
}
