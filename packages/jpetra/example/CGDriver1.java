import Jpetra.*;
import matvis.*;
import JpetraVis.*;

import java.awt.*;

class CGDriver1 {

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
	String matrixString = A.toMatrixMarketString();
	System.out.println(matrixString);

	MatView viewA = new MatView(matrixString);
	viewA.display();
	
	int [] nonzeroCount = A.getNumEntriesPerRow();
	
	VecView viewNonzero = new VecView(nonzeroCount);
	viewNonzero.setColor(Color.green);
	viewNonzero.histogram();
	viewNonzero.xygraph();

	Jpetra.Vector x = new Jpetra.Vector(map);
	Jpetra.Vector b = new Jpetra.Vector(map);
	Jpetra.Vector xexact = new Jpetra.Vector(map);
	Jpetra.Vector residual = new Jpetra.Vector(map);

	xexact.random(); // fill exact solution with random numbers
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
	x = solver.iterate(200,1.0e-10);
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

}

