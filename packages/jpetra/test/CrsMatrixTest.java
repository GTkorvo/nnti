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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

import Jpetra.*;

public class CrsMatrixTest {

    public static void main(String [] args) {

	int ierr=0, i, j;
	boolean debug = false;

	// #ifdef MPI...
	// #else
	
	int size = 1;
	int rank = 0;

	// #endif

	boolean verbose = false;

	// Should we print results to standard out?
	if(args.length > 0 && args[0].equals("-v")) verbose = true;

	// #ifdef MPI...
	// #else
	
	Jpetra.SerialComm comm = new Jpetra.SerialComm();

	// #endif

	int pid = comm.getVnodeID();
	int numProc = comm.getNumVnodes();
	if(verbose) System.out.println("Processor "+pid+" of "+numProc+" is alive.");

	boolean verbose1 = verbose;
	
	// Redefine verbose to only printon PE 0
	if(verbose && rank != 0) verbose = false;

	int numVnodeEquations = 10000;
	int numGlobalEquations = numVnodeEquations * numProc + Math.min(numProc, 3);
	if(pid < 3) numVnodeEquations++;
	int indexBase = 0;
	int elementSize = 7;
	boolean isDistributedGlobal = numGlobalEquations > numVnodeEquations;

	// Construct a map that puts approximately the same number of equations on each processor
	Jpetra.Map map = new Jpetra.Map(numGlobalEquations, numVnodeEquations, 0, comm);

	// Get and update list and number of local equations from newly created map
	int [] nodeGlobalElements = new int[map.getNumMyElements()];
	map.getMyGlobalElements(nodeGlobalElements);

	// Create an integer vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

	int [] numNz = new int[numVnodeEquations];

	// We are building a tridiagonal matrix where each row has (-1 2 -1)
	// So we need 2 off-diagonal terms (except for the first and last equation)
	for(i=0; i<numVnodeEquations; i++)
	    if(nodeGlobalElements[i] == 0 || nodeGlobalElements[i] == numGlobalEquations - 1)
		numNz[i] = 1;
	    else
		numNz[i] = 2;

	// Create a CrsMatrix
	Jpetra.CrsMatrix a = new Jpetra.CrsMatrix("Copy", map, numNz);
	if(a.indicesAreGlobal() || a.indicesAreLocal()) {
	    System.out.println("Error in Jpetra.CrsMatrix(String, Jpetra.Map, int[])");
	    System.exit(1);
	}
	// Add  rows one-at-a-time
	// Need some vectors to help
	// Off diagonal Values will always be -1

	double [] values = new double[2];
	values[0] = -1.0; values[1] = -1.0;
	int [] indices = new int[2];
	double two = 2.0;
	double [] twoa = {2.0};
	int numEntries;

	for(i=0; i<numVnodeEquations; i++) {
	    if(nodeGlobalElements[i] == 0) {
		indices[0] = 1;
		numEntries = 1;
	    }
	    else if(nodeGlobalElements[i] == numGlobalEquations-1) {
		indices[0] = numGlobalEquations - 2;
		numEntries = 1;
	    }
	    else {
		indices[0] = nodeGlobalElements[i] - 1;
		indices[1] = nodeGlobalElements[i] + 1;
		numEntries = 2;
	    }
	    ierr = a.insertGlobalValues(nodeGlobalElements[i], numEntries, values, indices);
	    if(ierr != 0) {
		System.out.println("insertGlobalValues: "+ierr);
		System.exit(0);
	    }
	    int [] tmp = new int[nodeGlobalElements.length-i];
	    System.arraycopy(nodeGlobalElements, i, tmp, 0, nodeGlobalElements.length-i);
	    ierr = a.insertGlobalValues(nodeGlobalElements[i], 1, twoa, tmp);

	    if(ierr <= 0) {
		System.out.println("insertGlobalValues: "+ierr);
		System.exit(0);
	    }
	}

	a.print();

	// Finish up
	String e = "";
	if(!a.indicesAreGlobal()) e = "indicesAreGlobal";
	if(a.transformToLocal() != 0) e = "transformToLocal";
	if(!a.indicesAreLocal()) e = "indicesAreLocal";
	if(a.isUpperTriangular()) e = "isUpperTriangular";
	if(a.isLowerTriangular()) e = "isLowerTriangular";
	if(!e.equals("")) {
	    System.out.println("Error: "+e);
	    System.exit(0);
	}

	int numVnodeNonzeros = 3 * numVnodeEquations;
	if(a.LRID(0) >= 0) numVnodeNonzeros--;
	if(a.LRID(numGlobalEquations-1) >= 0) numVnodeNonzeros--;

	if(check(a, numVnodeEquations, numGlobalEquations, numVnodeNonzeros, 3*numGlobalEquations-2,
		 nodeGlobalElements, verbose) != 0) {
	    System.out.println("Error in check");
	    System.exit(0);
	}
    }

    public static int check(Jpetra.CrsMatrix a, int numVnodeRows1, int numGlobalRows1, int numVnodeNonzeros1,
			    int numGlobalNonzeros1, int [] nodeGlobalElements, boolean verbose) {
	return 0;
    }

    public static int powerMethod(boolean transA, Jpetra.CrsMatrix a, Jpetra.Vector q,
				  Jpetra.Vector z, Jpetra.Vector resid, double [] lambda,
				  int niters, double tolerance, boolean verbose) {
	return 0;
    }

}
