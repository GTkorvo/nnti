package test;

import Jpetra.*;

/**
 * Modified BasicPerTest from epetra/test.
 * @author  Jason Cross
 */
public class BasicPerfTestSerial extends JpetraObject {
    public static final boolean DO_TRANSPOSE = false;
    public BasicPerfTestSerial(boolean verbose, boolean summary, int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints) {
        int j;
        double elapsed_time;
        double total_flops;
        double MFLOPs;
        
        initializeOutput();
        this.outputStreams.put("SUMMARY", new Output("", summary, System.out, false, System.out));
        this.setRootPrint("VERBOSE", verbose);
        
        Comm comm = new SerialComm();
        
        if (verbose || (summary && comm.getNumVnodes()==1)) {
            System.out.println(" Number of local nodes in X direction  = " + numNodesX);
            System.out.println(" Number of local nodes in Y direction  = " + numNodesY);
            System.out.println(" Number of global nodes in X direction = " + numNodesX*numProcsX);
            System.out.println(" Number of global nodes in Y direction = " + numNodesY*numProcsY);
            System.out.println(" Number of local nonzero entries       = " + numNodesX*numNodesY*numPoints);
            System.out.println(" Number of global nonzero entries      = " + numNodesX*numNodesY*numPoints*numProcsX*numProcsY);
            System.out.println(" Number of Processors in X direction   = " + numProcsX);
            System.out.println(" Number of Processors in Y direction   = " + numProcsY);
            System.out.println(" Number of Points in stencil           = " + numPoints + "\n");
        }
        
        if (summary && comm.getNumVnodes()>1) {
            System.out.print("\n\n\n\n\n\n\n\n\n\n");
        }
        
        if (numProcsX * numProcsY != comm.getNumVnodes()) {
            System.err.println("Number of processors = " + comm.getNumVnodes());
            System.err.println(" is not the product of " + numProcsX + " and " + numProcsY + "\n");
            System.exit(-1);
        }
        
        if (numPoints!=5 && numPoints!=9) {
            System.err.println("Number of points specified = " + numPoints + "\n");
            System.err.println(" is not 5 or 9");
            System.exit(-1);
        }
        
        if (numNodesX * numNodesY <= 0) {
            System.err.println("Product of number of nodes is <= zero\n");
            System.exit(-1);
        }
        
        int[] Xoff;
        int[] Yoff;
        if (numPoints == 5) {
            // Generate a 5-point 2D Finite Difference matrix
            Xoff = new int[5];
            Yoff = new int[5];
            Xoff[0] = -1; Xoff[1] = 1; Xoff[2] = 0; Xoff[3] = 0;  Xoff[4] = 0;
            Yoff[0] = 0;  Yoff[1] = 0; Yoff[2] = 0; Yoff[3] = -1; Yoff[4] = 1;
        }
        else {
            // Generate a 9-point 2D Finite Difference matrix
            Xoff = new int[9];
            Yoff = new int[9];
            Xoff[0] = -1;  Xoff[1] =  0; Xoff[2] =  1;
            Yoff[0] = -1;  Yoff[1] = -1; Yoff[2] = -1;
            Xoff[3] = -1;  Xoff[4] =  0; Xoff[5] =  1;
            Yoff[3] =  0;  Yoff[4] =  0; Yoff[5] =  0;
            Xoff[6] = -1;  Xoff[7] =  0; Xoff[8] =  1;
            Yoff[6] =  1;  Yoff[7] =  1; Yoff[8] =  1;
        }
        
        
        
        Object[] objects = GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints, Xoff, Yoff, comm);
        
        VectorSpace map = (VectorSpace) objects[0];
        CisMatrix A = (CisMatrix) objects[1];
        MultiVector b = (MultiVector) objects[2];
        MultiVector bt = (MultiVector) objects[3];
        MultiVector xexact = (MultiVector) objects[4];
        
        int length = map.getNumMyGlobalEntries();
        MultiVector q = new MultiVector(map, new double[1][length]);
        MultiVector z = new MultiVector(map, new double[1][length]);
        MultiVector r = new MultiVector(map, new double[1][length]);
        
        //Timings
        FlopCounter flopcounter = new FlopCounter();
        A.setFlopCounter(flopcounter);
        Time timer = new Time(comm);
        
        for (j=0; j<2; j++) { // j = 0 is notrans, j = 1 is trans
            
            flopcounter.resetFlops();
            timer.resetStartTime();
            
            boolean TransA = (j==1);
            if (TransA) {
                TransA = CisMatrix.USE_TRANSPOSE_A;
            }
            else {
                TransA = CisMatrix.USE_A;
            }
            
            if ((DO_TRANSPOSE && (TransA == CisMatrix.USE_TRANSPOSE_A)) || (TransA == CisMatrix.USE_A)) {
                
                //10 matvecs
                for(int i=0; i < 10; i++) {
                    A.multiply(TransA, xexact, z); // Compute z = A*xexact or z = A'*xexact
                }
                
                elapsed_time = timer.getElapsedTime();
                total_flops = A.getGlobalMegaFlops(comm);
                MFLOPs = total_flops/elapsed_time;
                this.println("VERBOSE", "\n\nTotal MFLOPs for 10 MatVec's (Trans = " + TransA
                + ")     = " + MFLOPs);
                if (comm.getNumVnodes() == 1) {
                    
                    if (TransA == CisMatrix.USE_TRANSPOSE_A) {
                        this.print("SUMMARY", "TransMV\t");
                    }
                    else {
                        this.print("SUMMARY", "NoTransMV\t");
                    }
                }
                this.println("SUMMARY", MFLOPs + "");

                // Compute residual
                if (TransA == CisMatrix.USE_TRANSPOSE_A) {
                    r.update(-1.0, z, 1.0, bt, 0.0); // r = bt - z
                }
                else {
                    r.update(-1.0, z, 1.0, b, 0.0); // r = b - z
                }
                
                double[] rnorm = r.norm2();
                this.println("VERBOSE", "Norm of difference between computed and exact RHS = " + rnorm[0]);
            }
        }
        
        q.setFlopCounter(flopcounter);
        z.setFlopCounter(flopcounter);
        
        flopcounter.resetFlops();
        timer.resetStartTime();
        //10 norms
        double[] n_out;
        for(int i=0; i < 10; i++) {
            n_out = q.norm2();
        }
        
        elapsed_time = timer.getElapsedTime();
        total_flops = q.getGlobalMegaFlops(comm);
        MFLOPs = total_flops/elapsed_time;
        this.println("VERBOSE", "\n\nTotal MFLOPs for 10 Norm2's= " + MFLOPs);
        
        
        if (comm.getNumVnodes() == 1) {
            this.print("SUMMARY", "Norm2\t");
        }
        this.println("SUMMARY", MFLOPs + "");
        
        
        r.setFlopCounter(flopcounter);
        
        flopcounter.resetFlops();
        timer.resetStartTime();
        
        //10 dot's
        for(int i=0; i < 10; i++) {
            n_out = q.dot(z);
        }
        
        elapsed_time = timer.getElapsedTime();
        total_flops = q.getGlobalMegaFlops(comm);
        MFLOPs = total_flops/elapsed_time;
        this.println("VERBOSE", "\n\nTotal MFLOPs for 10 Dot's= " + MFLOPs);
        
        if (comm.getNumVnodes() == 1) {
            this.print("SUMMARY", "DotProd\t");
        }
        this.println("SUMMARY", MFLOPs + "");
        
        flopcounter.resetFlops();
        timer.resetStartTime();
        
        //10 dot's
        for(int i=0; i < 10; i++)
            q.update(1.0, z, 1.0, r, 0.0);
        
        elapsed_time = timer.getElapsedTime();
        total_flops = q.getGlobalMegaFlops(comm);
        MFLOPs = total_flops/elapsed_time;
        this.println("VERBOSE", "\n\nTotal MFLOPs for 10 Updates= " + MFLOPs);
        
        if (comm.getNumVnodes() == 1) {
            this.print("SUMMARY", "Update\t");
        }
        this.println("SUMMARY", MFLOPs + "");
        
        System.exit(0);
    }
    
    // Constructs a 2D PDE finite difference matrix using the list of x and y offsets.
    //
    // nx      (In) - number of grid points in x direction
    // ny      (In) - number of grid points in y direction
    //   The total number of equations will be nx*ny ordered such that the x direction changes
    //   most rapidly:
    //      First equation is at point (0,0)
    //      Second at                  (1,0)
    //       ...
    //      nx equation at             (nx-1,0)
    //      nx+1st equation at         (0,1)
    
    // numPoints (In) - number of points in finite difference stencil
    // xoff    (In) - stencil offsets in x direction (of length numPoints)
    // yoff    (In) - stencil offsets in y direction (of length numPoints)
    //   A standard 5-point finite difference stencil would be described as:
    //     numPoints = 5
    //     xoff = [-1, 1, 0,  0, 0]
    //     yoff = [ 0, 0, 0, -1, 1]
    
    // nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed
    
    // comm    (In) - an Epetra_Comm object describing the parallel machine (numProcs and my proc ID)
    // map    (Out) - Epetra_Map describing distribution of matrix and vectors/multivectors
    // A      (Out) - Epetra_CrsMatrix constructed for nx by ny grid using prescribed stencil
    //                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
    //                diagonal will be slightly diag dominant.
    // b      (Out) - Generated RHS.  Values satisfy b = A*xexact
    // bt     (Out) - Generated RHS.  Values satisfy b = A'*xexact
    // xexact (Out) - Generated exact solution to Ax = b and b' = A'xexact
    
    // Note: Caller of this function is responsible for deleting all output objects.
    
    Object[] GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints,
    int[] xoff, int[] yoff,
    Comm comm) {
        
        MultiVector b1 = null;
        MultiVector bt1 = null;
        MultiVector xexact1 = null;
        
        return GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints,
        xoff, yoff, 1, comm);
    }
    
    Object[] GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, int[] xoff, int[] yoff, int nrhs, Comm comm) {
        VectorSpace map;
        CisMatrix A;
        MultiVector b;
        MultiVector bt;
        MultiVector xexact;
        
        Time timer = new Time(comm);
        // Determine my global IDs
        int[] myGlobalElements = GenerateMyGlobalElements(numNodesX, numNodesY, numProcsX, numProcsY, comm.getVnodeId());
        
        int numMyEquations = numNodesX*numNodesY;
        
        ElementSpace tempES = new ElementSpace(myGlobalElements, comm);
        map = new VectorSpace(tempES); // Create map with 2D block partitioning.
        myGlobalElements = null;
        
        int numGlobalEquations = map.getNumGlobalEntries();
        
        A = new CisMatrix(map, CisMatrix.ROW_ORIENTED); // Construct matrix
        
        int[] indices = new int[numPoints];
        double[] values = new double[numPoints];
        
        double dnumPoints = numPoints;
        int nx = numNodesX*numProcsX;
        
        RandomNumberGenerator random = new RngPackRandom(RngPackRandom.MT, 0.0, 10.0);
        for (int i=0; i<numMyEquations; i++) {
            
            int rowID = map.getGlobalIndex(i);
            int numIndices = 0;
            
            
            for (int j=0; j<numPoints; j++) {
                int colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
                if (colID>-1 && colID<numGlobalEquations) {
                    indices[numIndices] = colID;
                    double value = - random.getDouble();
                    if (colID==rowID)
                        values[numIndices++] = dnumPoints - value; // Make diagonal dominant
                    else
                        values[numIndices++] = -value;
                }
            }
            //cout << "Building row " << rowID << endl;
            A.insertEntries(i, indices, values, DistObject.REPLACE);
        }
        
        indices = null;
        values = null;
        double insertTime = timer.getElapsedTime();
        timer.resetStartTime();
        A.fillComplete();
        double fillCompleteTime = timer.getElapsedTime();
        
        this.println("VERBOSE", "Time to insert matrix values = " + insertTime);
        this.println("VERBOSE", "Time to complete fill        = " + fillCompleteTime);
        if (comm.getNumVnodes() == 1) this.print("SUMMARY", "InsertTime\t");
        this.println("SUMMARY", insertTime + "");
        if (comm.getNumVnodes() == 1) this.print("SUMMARY", "FillCompleteTime\t");
        this.println("SUMMARY", fillCompleteTime + "");
        
        int length = map.getNumMyGlobalEntries();
        if (nrhs<=1) {
            b = new MultiVector(map, new double[1][length]);
            bt = new MultiVector(map, new double[1][length]);
            xexact = new MultiVector(map, new double[1][length]);
        }
        else {
            b = new MultiVector(map, new double[nrhs][length]);
            bt = new MultiVector(map, new double[nrhs][length]);
            xexact = new MultiVector(map, new double[nrhs][length]);
        }
        
        xexact.putRandom(random); // Fill xexact with random values
        
        A.multiply(CisMatrix.USE_A, xexact, b);
        if (DO_TRANSPOSE) A.multiply(CisMatrix.USE_TRANSPOSE_A, xexact, bt);
        
        return new Object[]{map, A, b, bt, xexact};
    }
    
    int[] GenerateMyGlobalElements(int numNodesX, int numNodesY, int numProcsX, int numProcs, int myPID) {
        
        int[] myGlobalElements = new int[numNodesX*numNodesY];
        int myProcX = myPID%numProcsX;
        int myProcY = myPID/numProcsX;
        int curGID = myProcY*(numProcsX*numNodesX)*numNodesY+myProcX*numNodesX;
        for (int j=0; j<numNodesY; j++) {
            for (int i=0; i<numNodesX; i++) {
                myGlobalElements[j*numNodesX+i] = curGID+i;
            }
            curGID+=numNodesX*numProcsX;
        }
        //for (int i=0; i<numNodesX*numNodesY; i++) cout << "MYPID " << myPID <<" GID "<< myGlobalElements[i] << endl;
        
        return myGlobalElements;
    }
    
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        boolean verbose = false;
        boolean summary = false;
        
        // Check if we should print verbose results to standard out
        if ((args.length > 5) && (args[6].charAt(0) == '-' ) && (args[6].charAt(1) == 'v')) {
            verbose = true;
        }
        
        // Check if we should print verbose results to standard out
        if ((args.length > 5) && (args[6].charAt(0) == '-') && (args[6].charAt(1) == 's')) {
            summary = true;
        }
        
        if(args.length < 5) {
            System.err.println("Usage: ");
            System.err.println(" NumNodesX NumNodesY NumProcX NumProcY NumPoints [-v|-s]");
            System.err.println("where:");
            System.err.println("NumNodesX         - Number of mesh nodes in X direction per processor");
            System.err.println("NumNodesY         - Number of mesh nodes in Y direction per processor");
            System.err.println("NumProcX          - Number of processors to use in X direction");
            System.err.println("NumProcY          - Number of processors to use in Y direction");
            System.err.println("NumPoints         - Number of points to use in stencil (5 or 9 only)");
            System.err.println("-v|-s             - (Optional) Run in verbose mode if -v present or summary mode if -s present");
            System.err.println(" NOTES: NumProcX*NumProcY must equal the number of processors used to run the problem. Example:");
            System.err.println("mpirun -np 32  << argv[0] <<  10 12 4 8 -v");
            System.err.println(" Run this program on 32 processors putting a 10 X 12 subgrid on each processor using 4 processors ");
            System.err.println(" in the X direction and 8 in the Y direction.  Total grid size is 40 points in X and 96 in Y.");
            System.exit(-1);
        }
        
        int numNodesX = Integer.parseInt(args[0]);
        int numNodesY = Integer.parseInt(args[1]);
        int numProcsX = Integer.parseInt(args[2]);
        int numProcsY = Integer.parseInt(args[3]);
        int numPoints = Integer.parseInt(args[4]);
        
        new BasicPerfTestSerial(verbose, summary, numNodesX, numNodesY, numProcsX, numProcsY, numPoints);
    }
}
