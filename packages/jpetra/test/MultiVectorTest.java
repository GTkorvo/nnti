import Jpetra.*;

class MultiVectorTest {

    public static void main(String [] args) {
      try {
	int ierr = 0;
	int i;
	int j;

	int size = 1;
	int rank = 0;
	boolean verbose = false;
	if(args.length > 0 && args[0].equals("-v")) verbose = true;

	Jpetra.SerialComm comm = new Jpetra.SerialComm();
	int pid = comm.getPID();
	int numProc = comm.getNumProc();

	if(verbose) System.out.println("Processor "+pid+" of "+numProc+" is alive.");

	boolean verbose1 = verbose;
	if(verbose && rank != 0) verbose = false;
	
	int numNodeElements = 10000;
	int numNodeElements1 = numNodeElements;
	int numGlobalElements = numNodeElements * numProc + Math.min(numProc, 3);
	if(pid < 3) numNodeElements++;
	int indexBase = 0;
	int elementSize = 7;
	boolean isDistributedGloal = (numGlobalElements > numNodeElements);
	int numVectors = 4;
	
	// Test LocalMap constructor
	if(verbose) System.out.println("Checking LocalMap(numNodeElements1, indexBase, comm)");
	Jpetra.LocalMap localMap = new Jpetra.LocalMap(numNodeElements1, indexBase, comm);

	// Test uniform linear linear distribution constructor
	if(verbose) System.out.println("Checking BlockMap(numGlobalElements, elementSize, indexBase, comm)");
	Jpetra.BlockMap blockMap = new Jpetra.BlockMap(numGlobalElements, elementSize, indexBase, comm);
	ierr = MultiVectorTests(blockMap, numVectors, verbose);
	if(verbose) {
	    if(ierr == 0) System.out.println("Checked OK");
	    else System.out.println("Error code: "+ierr);
	}
	if(ierr != 0) System.exit(0);

	ierr = MatrixTests(blockMap, localMap, numVectors, verbose);
	if(verbose) {
	    if(ierr == 0) System.out.println("Checked OK");
	    else System.out.println("Error code: "+ierr);
	}
	if(ierr != 0) System.exit(0);

	// Test user defined linear distribution constructor

      } catch (Jpetra.JpetraException e) {
	  System.out.println(e);
      }
	System.exit(0);
    }

    public static int MultiVectorTests(Jpetra.BlockMap map, int numVectors, boolean verbose)
    throws Jpetra.JpetraException {
	Jpetra.Comm comm = map.getComm();
	int ierr = 0;
	double [] residual = new double [numVectors];
	int pid = comm.getPID();
	Jpetra.MultiVector a = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector sqrtA = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector b = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector c = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector cAlphaA = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector cAlphaAPlusB = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector cPlusB = new Jpetra.MultiVector(map, numVectors);
	Jpetra.MultiVector weights = new Jpetra.MultiVector(map, numVectors);

	// Construct double vectors;
	double [] dotVecAB = new double [numVectors];
	double [] norm1A = new double [numVectors];
	double [] norm2SqrtA = new double [numVectors];
	double [] normInfA = new double [numVectors];
	double [] normWA = new double [numVectors];
	double [] minvalA = new double [numVectors];
	double [] maxvalA = new double [numVectors];
	double [] meanvalA = new double [numVectors];

	// Generate data
	c.random();
	double alpha = 2.0;

	BuildMultiVectorTests(c, alpha, a, sqrtA, b, cAlphaA,
			      cAlphaAPlusB, cPlusB, dotVecAB,
			      norm1A, norm2SqrtA, normInfA,
			      normWA, weights, minvalA,
			      maxvalA, meanvalA);

	// Test alpha * a
	if(verbose) System.out.println("XXXXX  Testing alph * a XXXXX");
	Jpetra.MultiVector alphaA = new Jpetra.MultiVector(a);
	ierr += alphaA.scale(alpha);
	ierr += alphaA.update(-1.0, cAlphaA, 1.0);
	ierr += alphaA.norm2(residual);
	if(ierr != 0 && verbose) {
	    System.out.println("Error in alpha * a");
	}
	if(ierr != 0) return -2;
	if(verbose) {
	    System.out.println("alpha = "+alpha);
	    for(int j=0; j<numVectors; j++) {
		System.out.println("Residual["+j+"] = "+residual[j]);
	    }
	}
	
	// Test alpha * a + b
	if(verbose) {
	    System.out.println("XXXXX  Testing c = alpha * a + b  XXXXX");
	}
	Jpetra.MultiVector alphaAPlusB = new Jpetra.MultiVector(a);
	ierr += alphaAPlusB.update(1.0, b, alpha, a, 0.0);
	ierr += alphaAPlusB.update(-1.0, cAlphaAPlusB, 1.0);
	ierr += alphaAPlusB.norm2(residual);
	if(ierr != 0 && verbose) {
	    System.out.println("Error in alpha * a + b");
	}
	if(ierr != 0) return -2;
	if(verbose) {
	    System.out.println("  alpha = "+alpha);
	    for(int j=0; j<numVectors; j++) {
		System.out.println("Residual["+j+"] = "+residual[j]);
	    }
	}
	if(badResidual(residual, numVectors)) return -1;


	return 0;
    }

    public static int BuildMultiVectorTests(Jpetra.MultiVector c, double alpha, Jpetra.MultiVector a,
				     Jpetra.MultiVector sqrtA, Jpetra.MultiVector b,
				     Jpetra.MultiVector cAlphaA, Jpetra.MultiVector cAlphaAPlusB,
				     Jpetra.MultiVector cPlusB, double [] dotvecAB,
				     double [] norm1A, double [] norm2SqrtA, double [] norminfA,
				     double [] normWA, Jpetra.MultiVector weights, double [] minvalA,
				     double [] maxvalA, double [] meanvalA) {
	int i;
	int j;
	double fi;
	double fj;
	int aRows = a.getNodeLength();
	int aCols = a.getNumVectors();
	int sqrtARows = sqrtA.getNodeLength();
	int sqrtACols = sqrtA.getNumVectors();
	int bRows = b.getNodeLength();
	int bCols = b.getNumVectors();

	double [][] ap;
	double [][] sqrtAp;
	double [][] bp;
	double [][] cp;
	double [][] cAlphaAp;
	double [][] cAlphaAPlusBp;
	double [][] cPlusBp;
	double [][] weightsp;

	ap = a.extractView();
	sqrtAp = sqrtA.extractView();
	bp = b.extractView();
	cp = c.extractView();
	cAlphaAp = cAlphaA.extractView();
	cAlphaAPlusBp = cAlphaAPlusB.extractView();
	cPlusBp = cPlusB.extractView();
	weightsp = weights.extractView();

	boolean aIsLocal = (a.getNodeLength() == a.getGlobalLength());
	boolean bIsLocal = (b.getNodeLength() == b.getGlobalLength());
	boolean cIsLocal = (c.getNodeLength() == c.getGlobalLength());

	int aIndexBase = a.getMap().getIndexBase();
	int bIndexBase = b.getMap().getIndexBase();

	Jpetra.Map aMap = new Jpetra.Map(-1, aRows, aIndexBase, a.getMap().getComm());
	Jpetra.Map bMap = new Jpetra.Map(-1, bRows, bIndexBase, b.getMap().getComm());
	int [] aGlobalElements = new int [aRows];
	aMap.getGlobalElements(aGlobalElements);
	int [] bGlobalElements = new int [bRows];
	aMap.getGlobalElements(bGlobalElements);

	if(c.getNodeLength() != aRows ||
	   aRows             != bRows ||
	   c.getNumVectors() != aCols ||
	   aCols             != bCols ||
	   sqrtARows         != aRows ||
	   sqrtACols         != aCols ||
	   c.getNodeLength() != cAlphaA.getNodeLength() ||
	   c.getNumVectors() != cAlphaA.getNumVectors() ||
	   c.getNodeLength() != cAlphaAPlusB.getNodeLength() ||
	   c.getNodeLength() != cPlusB.getNodeLength() ||
	   c.getNumVectors() != cPlusB.getNumVectors() ) {
	    return -2;
	}
	
	boolean case1 =  aIsLocal &&  bIsLocal &&  cIsLocal;
	boolean case2 = !aIsLocal && !bIsLocal && !cIsLocal;

	if(!(case1 || case2)) return -3;

	double sf = a.getGlobalLength();
	double sfinv = 1.0 / sf;

	// Define a
	if(aIsLocal) {
	    for(i=0; i<aCols; i++) {
		for(j=0; j<aRows; j++) {
		    fi = i+1;
		    fj = j+1;
		    
		    ap[i][j] = (fi*sfinv)*fj;
		    sqrtAp[i][j] = Math.sqrt(ap[i][j]);
		}
	    }
	}
	else {
	    for(j=0; j<aCols; j++) {
		for(i=0; i<aRows; i++) {
		    fi = aGlobalElements[i] + 1;
		    fj = j + 1;
		    ap[j][i] = (fi*sfinv)*fj;
		    sqrtAp[j][i] = Math.sqrt(ap[j][i]);
		}
	    }
	}

	// Define b
	if(bIsLocal) {
	    for(j=0; j<bRows; j++) {
		for(i=0; i<bCols; i++) {
		    fi = i + 1;
		    fj = j + 1;
		    bp[i][j] = 1.0 / ((fi*sfinv)*fj);
		}
	    }
	}
	else {
	    for(j=0; j<bRows; j++) {
		for(i=0; i<bCols; i++) {
		    fi = bGlobalElements[i] + 1;
		    fj = j + 1;
		    bp[i][j] = 1.0 / ((fi*sfinv)*fj);
		}
	    }
	}

	// Generate cAlphaA
	for(j=0; j<aRows; j++) {
	    for(i=0; i<aCols; i++) {
		cAlphaAp[i][j] = alpha * ap[i][j];
	    }
	}

	// Generate cAlphaAPlusB
	for(j=0; j<aRows; j++) {
	    for(i=0; i<aCols; i++) {
		cAlphaAPlusBp[i][j] = alpha * ap[i][j] + bp[i][j];
	    }
	}

	// Generate cPlusB
	for(j=0; j<aRows; j++) {
	    for(i=0; i<aCols; i++) {
		cPlusBp[i][j] = cp[i][j] + bp[i][j];
	    }
	}

	// Generate dotvecAB
	for(i=0; i<a.getNumVectors(); i++) {
	    dotvecAB[i] = c.getGlobalLength();
	}

	double result = c.getGlobalLength();
	result *= sfinv;
	result /= 2.0;
	result *= (double)(c.getGlobalLength()+1);
	
	// Generate norm1A
	for(i=0; i<a.getNumVectors(); i++) {
	    norm1A[i] = result * ((double)(i+1));
	}

	// Generate norm2SqrtA
	for(i=0; i<a.getNumVectors(); i++) {
	    norm2SqrtA[i] = Math.sqrt(result *((double)(i+1)));
	}

	// Generate norminfA, minvalA, maxvalA, meanvalA
	for(i=0; i<a.getNumVectors(); i++) {
	    norminfA[i] = (double)(i+1);
	    minvalA[i] = (double)(i+1) / (double)a.getGlobalLength();
	    maxvalA[i] = (double)(i+i);
	    meanvalA[i] = norm1A[i] / ((double)(a.getGlobalLength()));
	}

	// Define weights and expected weight norm
	for(i=0; i<a.getNumVectors(); i++) {
	    double ip1 = (double) i+1;
	    normWA[i] = Math.sqrt(ip1);
	    for(j=0; j<aRows; j++) weightsp[i][j] = ip1*ap[i][j];
	}

	return 0;
    }

    public static int MatrixTests(Jpetra.BlockMap map, Jpetra.LocalMap localMap, int numVectors, boolean verbose) throws Jpetra.JpetraException {
	Jpetra.SerialComm comm = (Jpetra.SerialComm)map.getComm();
	int ierr = 0;
	int indexBase = 0;
	int i;
	double [] residual = new double [numVectors];
	int pid = comm.getPID();
	
	Jpetra.MultiVector a = new Jpetra.MultiVector(localMap, numVectors);
	Jpetra.MultiVector b = new Jpetra.MultiVector(localMap, numVectors);
	Jpetra.LocalMap map2d = new Jpetra.LocalMap(numVectors, indexBase, comm);
	Jpetra.MultiVector c = new Jpetra.MultiVector(map2d, numVectors);
	Jpetra.MultiVector cGemm = new Jpetra.MultiVector(map2d, numVectors);

	double [][] app;
	double [][] bpp;
	double [][] cpp;

	Jpetra.MultiVector ap;
	Jpetra.MultiVector bp;
	Jpetra.MultiVector cp;

	app = new double [numVectors][];
	bpp = new double [numVectors][];
	cpp = new double [numVectors][];

	for(i=0; i<numVectors; i++) {
	    if(a.getNodeLength() != b.getNodeLength()) {
		System.out.println(i+" a:"+a.getNodeLength()+", b:"+b.getNodeLength());
	    }
	    app[i] = new double [a.getNodeLength()+i];
	    bpp[i] = new double [b.getNodeLength()+i];
	    cpp[i] = new double [c.getNodeLength()+i];
	}

	Jpetra.MultiVector a1 = new Jpetra.MultiVector("View", localMap, app, numVectors);
	Jpetra.MultiVector b1 = new Jpetra.MultiVector("View", localMap, bpp, numVectors);
	Jpetra.MultiVector c1 = new Jpetra.MultiVector("View", map2d, cpp, numVectors);

	for(i=0; i<4; i++) {
	    ierr = 0;
	    String transa = "N"; if (i>1) transa = "T";
	    String transb = "N"; if (i%2!=0) transb = "T";
	    double alpha = (double)i+1;
	    double beta = (double) i/2.0;
	    ierr += c.random(); // Fill c with random numbers
	    ierr += buildMatrixTests(c, transa, transb, alpha, a, b, beta, cGemm);
	    if(ierr != 0) {
		System.out.println("*** build: "+ierr+" ***");
		System.exit(0);
	    }
	    
	    ap = a;
	    bp = b;
	    cp = c;

	    ierr += cp.multiply(transa, transb, alpha, ap, bp, beta);
	    ierr += cp.update(-1.0, cGemm, 1.0);
	    ierr += cp.norm2(residual);

	    //if(verbose && ierr == 0) {
	    if(verbose) {
		System.out.println(i);
		System.out.println("XXXXX Replicated Local MultiVector GEMM tests XXXXX");
		System.out.println("      alpha = "+alpha+", beta = "+beta+", transa = "
				   +transa+", transb = "+transb);
		for(int j=0; j<numVectors; j++) System.out.println("      Residual["+j+"] = "+residual[j]);
	    }
	    if(ierr == 0 && badResidual(residual, numVectors)) return -1;
	}

	return 0;
    }

    public static int buildMatrixTests(Jpetra.MultiVector c, String transa, String transb, double alpha,
			 Jpetra.MultiVector a, Jpetra.MultiVector b, double beta,
			 Jpetra.MultiVector cGemm) {
	int i, j;
	double fi, fj;

	double [][] ap = new double [a.getNodeLength()][a.getNumVectors()];
	double [][] bp = new double [b.getNodeLength()][b.getNumVectors()];
	double [][] cp = new double [c.getNodeLength()][c.getNumVectors()];
	double [][] cGemmp = new double [cGemm.getNodeLength()][cGemm.getNumVectors()];

	int aRows = a.getNodeLength();
	int aCols = a.getNumVectors();
	int bRows = b.getNodeLength();
	int bCols = b.getNumVectors();
	int cRows = c.getNodeLength();
	int cCols = c.getNumVectors();
	int cGemmRows = cGemm.getNodeLength();
	int cGemmCols = cGemm.getNumVectors();
	ap = a.extractView();
	bp = b.extractView();
	cp = c.extractView();
	cGemmp = cGemm.extractView();

	/* DEBUG
	int opACols = (transa.equals("N")) ? a.getNumVectors() : a.getNodeLength();
	int opBRows = (transb.equals("N")) ? b.getNodeLength() : b.getNumVectors();
	int cGlobalInnerDim = (transa.equals("N")) ? a.getNumVectors() : a.getGlobalLength();
	*/
	int opACols = (transa.equals("N")) ? a.getNodeLength() : a.getNumVectors();
	int opBRows = (transb.equals("N")) ? b.getNumVectors() : b.getNodeLength();
	int cGlobalInnerDim = (transa.equals("N")) ? a.getNodeLength() : a.getGlobalLength();

	System.out.println("transa: "+transa+", transb: "+transb);
	System.out.println("ACols: "+opACols+", BRows: "+opBRows);

	boolean aIsLocal = !a.isDistributedGlobal();
	boolean bIsLocal = !b.isDistributedGlobal();
	boolean cIsLocal = !c.isDistributedGlobal();
	int aIndexBase = a.getMap().getIndexBase();
	int bIndexBase = b.getMap().getIndexBase();

	Jpetra.Map aMap = new Jpetra.Map(-1, aRows, aIndexBase, a.getMap().getComm());
	Jpetra.Map bMap = new Jpetra.Map(-1, bRows, bIndexBase, b.getMap().getComm());
	
	int [] aGlobalElements = new int [aRows];
	aMap.getGlobalElements(aGlobalElements);
	int [] bGlobalElements = new int [bRows];
	bMap.getGlobalElements(bGlobalElements);

	// Check for compatible dimensions
	if(c.getNodeLength() != cRows ||
	   opACols           != opBRows ||
	   c.getNumVectors() != cCols ||
	   c.getNodeLength() != cGemm.getNodeLength() ||
	   c.getNumVectors() != cGemm.getNumVectors()) {

	    if(c.getNodeLength() != cRows) System.out.println(1);
	    if(opACols != opBRows) System.out.println(2);
	    if(c.getNumVectors() != cCols) System.out.println(3);
	    if(c.getNodeLength() != cGemm.getNodeLength()) System.out.println(4);
	    if(c.getNumVectors() != cGemm.getNumVectors()) System.out.println(5);
	    return -2;
	}
	boolean case1 =  aIsLocal &&  bIsLocal &&  cIsLocal;
	boolean case2 = !aIsLocal && !bIsLocal &&  cIsLocal;
	boolean case3 = !aIsLocal &&  bIsLocal && !cIsLocal;

	if(!(case1 || case2 || case3)) {
	    return -3;
	}
	
	double sf = cGlobalInnerDim;
	double sfinv = 1.0/sf;

	// Define a
	if(aIsLocal) {
	    for(j=0; j<aRows; j++) {
		for(i=0; i<aCols; i++) {
		    fi = i+1;
		    fj = j+1;
		    ap[j][i] = (fi*sfinv)*fj;
		}
	    }
	}
	else {
	    for(j=0; j<aCols; j++) {
		for(i=0; i<aRows; i++) {
		    fi = aGlobalElements[i]+1;
		    fj = j+1;
		    ap[j][i] = (fi*sfinv)*fj;
		}
	    }
	}

	// Define b
	if(bIsLocal) {
	    for(j=0; j<bRows; j++) {
		for(i=0; i<bCols; i++) {
		    fi = i+1;
		    fj = j+1;
		    bp[j][i] = 1.0/((fi*sfinv)*fj);
		}
	    }
	}
	else {
	    for(j=0; j<bCols; j++) {
		for(i=0; i<bRows; i++) {
		    fi = bGlobalElements[i]+1;
		    fj = j+1;
		    bp[j][i] = 1.0/((fi*sfinv)*fj);
		}
	    }
	}

	// Define cGemm
	if(case1) {
	    for(j=0; j<cCols; j++) {
		for(i=0; i<cRows; i++) {
		    fi = (i+1)*cGlobalInnerDim;
		    fj = j+1;
		    cGemmp[j][i] = alpha * (fi/fj) + beta * cp[j][i];
		}
	    }
	}
	else if(case2) {
	    for(j=0; j<cCols; j++) {
		for(i=0; i<cRows; i++) {
		    fi = (i+1)*cGlobalInnerDim;
		    fj = j+1;
		    cGemmp[j][i] = alpha * (fi/fj) + beta * cp[j][i];
		}
	    }
	}
	else {
	    for(j=0; j<cCols; j++) {
		for(i=0; i<cRows; i++) {
		    fi = (aGlobalElements[i]+1) * cGlobalInnerDim;
		    fj = j+1;
		    cGemmp[j][i] = alpha * (fi/fj) + beta * cp[j][i];
		}
	    }
	}
	return 0;
    }

    public static boolean badResidual(double [] residual, int numVectors) {
	double threshold = 5.0E-6;
	for(int i=0; i<numVectors; i++) {
	    if(residual[i] > threshold) return true;
	}
	return false;
    }
}
