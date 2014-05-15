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

import Jpetra.*;
import matvis.*;
import javax.swing.*;
import java.awt.*;
import JpetraVis.*;


public class CGDriver_play implements Runnable {
    public CGDriver_play() {}

    public void run() {
    }

    public static void main(String [] args) {
      try {
	int ierr = 0;
	int i;
	int j;

	int size = 1;
	int rank = 0;
	boolean verbose = true;
	Jpetra.SerialComm comm = new Jpetra.SerialComm();
	int pid = comm.getVnodeID();
	int numProc = comm.getNumVnodes();

	if(verbose) System.out.println("Processor "+pid+" of "+numProc+" is alive.");

	boolean verbose1 = verbose;
	if(verbose && rank != 0) verbose = false;
	if (args.length == 0) {
	    System.out.println("Enter path and name of a Harwell-Boeing matrix file");
	    System.exit(1);
	}

	matvis.HBReader hb = getMatrixFile(args[0]);
	
        int numEquations = hb.numMatrixRows;
	Jpetra.Map map = new Jpetra.Map(numEquations, 0, comm);
	Jpetra.CrsMatrix A = createMatrix(hb, map);
	System.out.println("Matrix in Matrix Market Format");
	String matrixMarketString = A.toMatrixMarketString();
	System.out.println(matrixMarketString);

	MatView myMatrix = new MatView(matrixMarketString);
	myMatrix.display();

	Jpetra.Vector x = new Jpetra.Vector(map);
	Jpetra.Vector b = new Jpetra.Vector(map);
	Jpetra.Vector xexact = new Jpetra.Vector(map);
	Jpetra.Vector residual = new Jpetra.Vector(map);

	xexact.putScalar(1.0); // fill exact solution with random numbers
	A.multiply(xexact,b); // Explicitly compute b to match A*xexact

        double [] normb = new double[1];
        double [] normxexact = new double[1];
        double [] normresidual = new double[1];
        b.norm2(normb);
	xexact.norm2(normxexact);
        
	System.out.println("Norm of b        " + normb[0]);
	System.out.println("Norm of x(exact) " + normxexact[0]);
	A.multiply(x,residual);
        residual.update(1.0, b, -1.0);
	residual.norm2(normresidual);
	System.out.println("Norm of residual prior to solution " + normresidual[0]);
	CGSolver solver = new CGSolver(A, x, b);
	
	for (int inum=0; inum<2; inum++) {
	    x = solver.iterate(1,1.0e-10);
	    VecView p_vec = new VecView(solver.getp(), "p");
	    VecView r_vec = new VecView(solver.getr(), "r");
	    p_vec.setColor(Color.magenta);
	    p_vec.setSize(250, 250);
	    p_vec.xygraph();
	    r_vec.setColor(Color.red);
	    r_vec.setSize(275, 275);
	    r_vec.xygraph();
	}
	VecView x_vec = new VecView(x, "x");
	x_vec.xygraph();
	
	    
	A.multiply(x,residual);
        residual.update(1.0, b, -1.0);
	residual.norm2(normresidual);
	System.out.println("Norm of residual after solution    " + normresidual[0]);

      } catch (Jpetra.JpetraException e) {
	  System.out.println(e);
      }
      
	//System.exit(0);
    }
    
    public static matvis.HBReader getMatrixFile(String filename) {
	matvis.HBReader hb = new matvis.HBReader(filename);
	System.out.println("hb.key: " + hb.key);
	System.out.println("ptrcrd: "+hb.linesForPtrs+
			   ", indcrd: "+hb.linesForRowIndices+
			   ", valcrd: "+hb.linesForNumericalValues+
			   ", rhscrd: "+hb.linesForRHSs);
	System.out.println("mxtype: "+hb.matrixType+
			   ", rows: "+hb.numMatrixRows+
			   ", cols: "+hb.numMatrixCols+
			   ", nnz: "+hb.numMatrixNonZeros);
	System.out.println("colsPerPtrEntry: "+hb.colsPerPtrEntry+
			   ", colsPerIndEntry: "+hb.colsPerIndEntry+
			   ", colsPerValEntry: "+hb.colsPerValEntry);
	return(hb);
    }
    public static Jpetra.CrsMatrix createMatrix(matvis.HBReader hb, Jpetra.Map map) {
     int[] ptrs = new int[hb.numMatrixCols+1];
     int[] rowInds = new int[hb.numMatrixNonZeros];
     float[] vals = new float[hb.numMatrixNonZeros];

     boolean success = hb.readData(ptrs, rowInds, vals);
     if (success) System.out.println("ptrs[0]: "+ptrs[0]+", ptrs["+
                         hb.numMatrixCols+"]: "+
                         ptrs[hb.numMatrixCols]+"\nrowInds[0]: "
                         +rowInds[0]+", rowInds["+
                         (hb.numMatrixNonZeros-1)+"]: "+
                         rowInds[hb.numMatrixNonZeros-1]+"\n"+
                         "vals[0]: "+vals[0]+", vals["+
                         (hb.numMatrixNonZeros-1)+"]: "+
                         vals[hb.numMatrixNonZeros-1]);
     else {
         System.out.println("hb.readData not successful");
         System.exit(1);
     }

     int numEquations = hb.numMatrixRows;
     int maxNumEntries = 0;
     int [] numEntriesVec = new int[numEquations];
     for(int i=0; i<numEquations; i++) {
	 int numEntries = ptrs[i+1] - ptrs[i];
         numEntriesVec[i] = numEntries;
	 if (numEntries > maxNumEntries) maxNumEntries = numEntries;
     }
	Jpetra.CrsMatrix A = new Jpetra.CrsMatrix("Copy", map, numEntriesVec);

     double [] values = new double[maxNumEntries];
     int [] indices = new int[maxNumEntries];
     int offset = 0;
     for(int i=0; i<numEquations; i++) {
	 int numEntries = ptrs[i+1] - ptrs[i];
	 for (int j=0; j < numEntries; j++) {
	     values[j] = vals[j+offset];
	     indices[j] = rowInds[j+offset]-1;
	 }
	 A.insertGlobalValues(i, numEntries, values, indices);
	 offset += numEntries;
     } 
     A.transformToLocal();
     return(A);
    }
    public static java.net.InetAddress getLocalInetAddress() {
        java.net.InetAddress iadd = null;
        try {
            iadd = java.net.InetAddress.getLocalHost();
        }
        catch(java.net.UnknownHostException uhe) {
            System.out.println("!! unknown-host-exception !!");
         return(null);
        }

     return(iadd);
     }
}

