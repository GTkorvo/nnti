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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

package test;

import Jpetra.*;
import java.util.Arrays;

public class CcjCommBroadcastTest {
    
    public static void main(String[] args) {
        new CcjCommBroadcastTest(args);
    }
    
    public CcjCommBroadcastTest(String[] args) {
        Jpetra.JpetraObject output = new JpetraObject();
        output.initializeOutput();
        
        Jpetra.CcjComm comm = new Jpetra.CcjComm("test/ccjhosts.txt");

        if(args.length > 0 && args[0].equals("-v")) output.setRootPrint("VERBOSE", true);
        
        output.println("VERBOSE", "Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
        
        double[] toSend = new double[3];
        if (comm.getVnodeID() == 0) {
            toSend = new double[]{1.5,2.6,3.7};
        }
        
        
        double[] toRecieve = (double[]) comm.broadcast(toSend,0);
        
        if (!java.util.Arrays.equals(toRecieve, new double[]{1.5,2.6,3.7})) {
            output.println("ERR", "Test Failed.");
            System.exit(1);
        }
        
        output.println("VERBOSE", toRecieve[0] + " " + toRecieve[1] + " " + toRecieve[2]);
        
        System.exit(0);
    }
}
