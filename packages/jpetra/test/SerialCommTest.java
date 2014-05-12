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

class SerialCommTest {

    public static void main (String [] args) {
	
	int size = 1;
	int rank = 0;
	boolean verbose = false;
	
	if(args.length > 0 && args[0].equals("-v")) verbose = true;

	Jpetra.SerialComm comm = new Jpetra.SerialComm();
	int myPid = comm.getVnodeId();
	int numProce = comm.getNumVnodes();
	if(verbose) System.out.println("Processor "+comm.getVnodeId()+" of "+comm.getNumVnodes());

	System.exit(0);
    }
}
