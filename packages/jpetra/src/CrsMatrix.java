package Jpetra;

/**
 *
 * @author Michael William Boldt
 * @version 
 */
public class CrsMatrix extends JpetraObject {

    protected Jpetra.Graph graph = null;
    protected boolean isAllocated = false;
    protected boolean hasStaticGraph = false;

    protected double [][] values = null;
    protected double normInf = -1;
    protected double normOne = -1;
    
    protected int numProcessRows = 0;
    protected int [] numEntriesPerRow = null;
    protected int [] numAllocatedEntriesPerRow = null;
    protected int [][] indices = null;
    Jpetra.MultiVector importVector = null;
    Jpetra.MultiVector exportData = null;
    int lenImports = 0;
    int lenExports= 0;
    double [] imports = null;
    double [] exports = null;
    int [] intImports = null;
    int [] intExports = null;
    String cv = null;

    /**
     * Creates a CrsMatrix object and allocates storage.
     *
     * @param cv               Copy or View
     *                         
     * @param rowMap           Map to use for this matrix
     *
     * @param numEntriesPerRow The ith entry indicates the approximate number of entries in the ith row
     */
    public CrsMatrix(String cv, final Jpetra.Map rowMap, int [] numEntriesPerRow) {
	this.numProcessRows = rowMap.getNumVnodeEquations();
	this.cv = cv;
	this.graph = new Jpetra.Graph(cv, rowMap, numEntriesPerRow);
	int ierr = allocate();
    }

    /**
     * Creates a CrsMatrix object and allocates storage.
     *
     * @param cv               Copy or View
     *
     * @param rowMap           Map to use for this matrix
     *
     * @param numEntriesPerRow Approximate number of entries in each row, 0 will result in fill during the insertion phase
     */
    public CrsMatrix(String cv, final Jpetra.Map rowMap, int numEntriesPerRow) {
	this.numProcessRows = rowMap.getNumVnodeEquations();
	this.cv = cv;
       // This constructor is not in Graph!!!	this.graph = new Jpetra.Graph(cv, rowMap, numEntriesPerRow);
	int ierr = allocate();
    }

    /**
     * Creates a CrsMatrix object and allocates storage.
     *
     * @param cv               Copy or View
     *
     * @param rowMap
     *
     * @param colMap
     *
     * @param numEntriesPerRow The ith entry indicates the approximate number of entries in the ith row
     */
    public CrsMatrix(String cv, final Jpetra.Map rowMap, final Jpetra.Map colMap, int [] numEntriesPerRow) {
	this.numProcessRows = rowMap.getNumVnodeEquations();
	this.cv = cv;
	this.graph = new Jpetra.Graph(cv, rowMap, colMap, numEntriesPerRow);
	int ierr = allocate();
    }

    /**
     * Creates a CrsMatrix object and allocates storage.
     *
     * @param cv               Copy or View
     *
     * @param rowMap
     *
     * @param colMap
     *
     * @param numEntriesPerRow Approximate number of entries in each row
     */
    public CrsMatrix(String cv, final Jpetra.Map rowMap, final Jpetra.Map colMap, int numEntriesPerRow) {
	this.numProcessRows = rowMap.getNumVnodeEquations();
	this.cv = cv;
	int [] tmp = { numEntriesPerRow };
	this.graph = new Jpetra.Graph(cv, rowMap, colMap, tmp);
	int ierr = allocate();
    }

    /**
     * Creates a CrsMatrix from a Graph.
     *
     * @param cv    Copy or View
     *
     * @param graph
     */
    public CrsMatrix(String cv, final Jpetra.Graph graph) {
	this.numProcessRows = graph.getNumProcessRows();
	this.cv = cv;
	int ierr = allocate();
    }

    /**
     * Copy constructor
     *
     * @param matrix CrsMatrix from which to make a copy
     */
    public CrsMatrix(CrsMatrix matrix) {
	this.isAllocated = matrix.isAllocated;
	this.numProcessRows = matrix.numProcessRows;
	this.cv = "Copy";
	this.graph = new Jpetra.Graph(matrix.graph);
	int ierr = allocate();
	for(int i=0; i<numProcessRows; i++) {
	    int numEntries = this.numEntriesPerRow[i];
	    for(int j=0; j<numEntries; j++) this.values[i][j] = matrix.values[i][j];
	}
    }

    protected int allocate() {

	int i, j;

	// Get direct references to graph info
	this.numEntriesPerRow = graph.getNumIndicesPerRow();
	this.numAllocatedEntriesPerRow = graph.getNumAllocatedIndicesPerRow();
	this.indices = graph.getIndices();

	// Allocate values array
	this.values = new double [this.numProcessRows][];

	if(this.cv.equals("Copy")) {
	    for(i=0; i<this.numProcessRows; i++) {
		int numAllocatedEntries = this.numAllocatedEntriesPerRow[i];
		
		if(numAllocatedEntries > 0) this.values[i] = new double[numAllocatedEntries];
		else this.values[i] = null;
		
		for(j=0; j<numAllocatedEntries; j++) this.values[i][j] = 0.0;
	    }
	}
	else {
	    for(i=0; i<this.numProcessRows; i++) this.values[i] = null;
	}

	setAllocated(true);
	return 0;
    }

    /**
     * Initialize all values in the graph of the matrix with <code>scalar</code>.
     */
    public int putScalar(double scalar) {
	for(int i=0; i<numProcessRows; i++) {
	    int numEntries = numEntriesPerRow[i];
	    for(int j=0; j<numEntries; j++) values[i][j] = scalar;
	}
	return 0;
    }

    public int insertGlobalValues(int row, int numEntries, double [] values, int [] indices) {
	if(indicesAreLocal()) return -2; // Cannot insert global values into local graph
	graph.setIndicesAreGlobal(true);
	row = graph.LRID(row); // Find local row number for this global row index

	return insertValues(row, numEntries, values, indices);
    }

    public int insertProcessValues(int row, int numEntries, double [] values, int [] indices) {
	if(indicesAreGlobal()) return -2;
	graph.setIndicesAreLocal(true);

	return insertValues(row, numEntries, values, indices);
    }

    public int insertValues(int row, int numEntries, double [] values, int [] indices) {
	if(hasStaticGraph()) return -2; // If the matrix graph is fully constructed, we can not insert new values

	int j;
	double [] tmpValues;
	int ierr = 0;

	if(row < 0 || row >= this.numProcessRows) return -1; // Not in row range
	
	if(this.cv.equals("View")) {
	    if(this.values[row] != null) ierr = 2; // This row has already been defined.  Issue warning.
	    this.values[row] = values;
	}
	else {
	    int start = this.numEntriesPerRow[row];
	    int stop = start + numEntries;
	    int numAllocatedEntries = this.numAllocatedEntriesPerRow[row];
	    if(stop > numAllocatedEntries) {
		if(numAllocatedEntries == 0) this.values[row] = new double[numEntries]; // Row was not allcoated, so do it.
		else {
		    ierr = 1; // Out of room.  Must delete and allocate more space.
		    tmpValues = new double[stop];
		    for(j=0; j<start; j++) tmpValues[j] = this.values[row][j]; // Copy existing entries
		    this.values[row] = tmpValues; // Set to new storage.
		}
	    }
	    for(j=start; j<stop; j++) this.values[row][j] = values[j-start];
	}
	
	graph.insertIndices(row, numEntries, indices); // Update graph
	
	return ierr;
    }

    public int transformToLocal() {
	return (transformToLocal(getRowMap(), getColMap()));
    }

    public int transformToLocal(Jpetra.BlockMap domainMap, Jpetra.BlockMap rangeMap) {
	if(!hasStaticGraph) graph.makeIndicesLocal(domainMap, rangeMap);
	sortEntries(); // Sort column entries from smalles to largest
	mergeRedundantEntries(); // Get rid of any redundant index values
	if(!hasStaticGraph) graph.transformToLocal(domainMap, rangeMap);

	return 0;
    }

    /**
     * Returns a copy of the specified local global in user-provided arrays.
     */
    public int extractGlobalRowCopy(int row, int length, int [] numEntries, double [] values, int [] indices) {
	int ierr = graph.extractGlobalRowCopy(row, length, numEntries, indices);
	if(ierr != 0) return ierr;
	return extractGlobalRowCopy(row, length, numEntries, values);
    }

    /**
     * Returns a copy of the specified global row values in user-provided array.
     */
    public int extractGlobalRowCopy(int row, int length, int [] numEntries, double [] values) {
	int row0 = graph.getRowMap().getLID(row); // Normalize row range
	return extractProcessRowCopy(row0, length, numEntries, values);
    }

    /**
     * Returns a copy of the specified local row in user-provided arrays.
     */
    public int extractProcessRowCopy(int row, int length, int [] numEntries, double [] values, int [] indices) {
	int ierr = graph.extractProcessRowCopy(row, length, numEntries, indices);
	if(ierr != 0) return ierr;
	return extractProcessRowCopy(row, length, numEntries, values);
    }

    /**
     * Returns a copy of the specified local row values in user-provided array.
     */
    public int extractProcessRowCopy(int row, int length, int [] numEntries, double [] values) {
	int j;

	if(row < 0 || row >= this.numProcessRows) return -1;
	
	numEntries[0] = this.numEntriesPerRow[row];
	if(length < numEntries[0]) return -2; // Not enough space for copy. Needed size is passed back in numEntries

	System.arraycopy(this.values[row], 0, values, 0, numEntries[0]);

	return 0;
    }

    public int sortEntries() {
	if(!indicesAreLocal()) return -1;
	if(isSorted()) return 0;

	// For each row, sort column entries from smallest to largest
	// Use shell sort. Stable sort so it is fast if indices are already sorted.

	for(int i=0; i<numProcessRows; i++) {
	    double [] values = this.values[i];
	    int numEntries = this.numEntriesPerRow[i];
	    int [] indices = this.indices[i];

	    int n = numEntries;
	    int m = n / 2;
	    
	    while(m > 0) {
		int max = n - m;
		for(int j=0; j<max; j++) {
		    for(int k=j; k>=0; k-=m) {
			if(indices[k+m] >= indices[k]) break;
			double dtemp = values[k+m];
			values[k+m] = values[k];
			values[k] = dtemp;
			int itemp = indices[k+m];
			indices[k+m] = indices[k];
			indices[k] = itemp;
		    }
		}
		m = m / 2;
	    }
	}
	graph.setSorted(true); // This also sorted the graph
	return 0;
    }

    public int mergeRedundantEntries() {
	int i, j, k;
	
	if(hasNoRedundancies()) return 0;
	if(!isSorted()) return -1; // Must have sorted entries

	// For each row, remove column indices that are repeated.
	// Also determine if matrix is upper or lower triangular or has no diagonal.
	// Note: Assumes that sortEntries was already called.

	for(i=0; i<numProcessRows; i++) {
	    int numEntries = this.numEntriesPerRow[i];
	    if(numEntries > 0) {
		final double [] values = this.values[i];
		final int [] indices = this.indices[i];
		int j0 = 0;
		int jj0 = indices[j0];
		for(j=0; j<numEntries; j++) {
		    int jj = indices[j];
		    if(jj == jj0) { // Check if index is repeated
			values[j0] += values[j];
			for(k=j; k<numEntries-1; k++) values[k] = values[k+1];
			numEntries--;
		    }
		    else {
			j0 = j; // Redefine comparison index values
			jj0 = indices[j0];
		    }
		}
	    }
	}

	return graph.removeRedundantIndices();
    }

    public int multiply(final Jpetra.Vector domainVector, Jpetra.Vector rangeVector) {

	final double [] x = domainVector.extractView()[0];
	//final double [] x = xmv[0];
	double [] y = rangeVector.extractView()[0];
	//double [] y = ymv[0];
	// should check if x.getNumVectors()==numVectors
	for(int i=0; i<numProcessRows; i++) {
	    int numEntries = numEntriesPerRow[i];
	    double sum = 0.0;
	    if(numEntries > 0) {
		final double [] vals = this.values[i];
		final int [] inds = this.indices[i];
		for(int j=0; j<numEntries; j++) {
		    sum += vals[j]*x[inds[j]];
		}
		y[i] = sum;
	    }
	}
	return(0);
    }

    public int [] getNumEntriesPerRow() { return graph.getNumIndicesPerRow(); } 

    public Jpetra.Map getRowMap() { return (Jpetra.Map)graph.getRowMap(); }

    public Jpetra.Map getColMap() { return (Jpetra.Map)graph.getColMap(); }

    protected boolean hasStaticGraph() { return hasStaticGraph; }

    protected int setAllocated(boolean flag) { isAllocated = flag; return 0; }

    /**
     * Returns the number of nonzero entries in this global matrix.
     */
    public int getNumGlobalNonzeros() { return graph.getNumGlobalNonzeros(); }

    /**
     * Returns the number of global matrix rows.
     */
    public int getNumGlobalRows() { return graph.getNumGlobalRows(); }

    /**
     * Returns the number of global matrix columns.
     */
    public int getNumGlobalCols() { return graph.getNumGlobalCols(); }

    /**
     * Returns the number of global nonzero diagonal entries.
     */
    public int getNumGlobalDiagonals() { return graph.getNumGlobalDiagonals(); }

    /**
     * Returns the number of nonzero entries in the calling process's portion of the matrix.
     */
    public int getNumProcessNonzeros() { return graph.getNumProcessNonzeros(); }

    /**
     * Returns the number of matrix rows owned by the calling process.
     */
    public int getNumProcessRows() { return graph.getNumProcessRows(); }

    /**
     * Returns the number of matrix columns owned by the calling process.
     */
    public int getNumProcessCols() { return graph.getNumProcessCols(); }

    /**
     * Returns the number of local nonzero diagonal entries.
     */
    public int getNumProcessDiagonals() { return getNumProcessDiagonals(); }

    /**
     * Returns the current number of nonzero entries in specified global row on this process.
     */
    public int getNumGloalEntries(int row) { return graph.getNumGlobalIndices(row); }

    /**
     * Returns the allocated number of nonzero entries in specified global row on this process.
     */
    public int getNumAllocatedGlobalEntries(int row) { return graph.getNumAllocatedGlobalIndices(row); }

    /**
     * Returns the maximum number of nonzero entries across all rows on this process.
     */
    public int getMaxNumEntries() { return graph.getMaxNumIndices(); }
    
    /**
     * If sortEntries() has been called, this query returns true, otherwise it returns false.
     */
    public boolean isSorted() { return graph.isSorted(); }

    /**
     * If mergeRedundantEntries() has been called, this query returns true, otherwise it returns false.
     */
    public boolean hasNoRedundancies() { return graph.hasNoRedundancies(); }

    /**
     * If matrix indices has been transformed to local, this query returns true, otherwise it returns false.
     */
    public boolean indicesAreLocal() { return graph.indicesAreLocal(); }

    /**
     * If matrix indices has not been transformed to local, this query returns true, otherwise it returns false.
     */
    public boolean indicesAreGlobal() { return graph.indicesAreGlobal(); }

    /**
     * If matrix is upper triangular, this query returns true, otherwise it returns false.
     */
    public boolean isUpperTriangular() { return graph.isUpperTriangular(); }

    /**
     * If matrix is lower triangular, this query returns true, otherwise it returns false.
     */
    public boolean isLowerTriangular() { return graph.isLowerTriangular(); }

    /**
     * If matrix is lower triangular, this query returns true, otherwise it returns false.
     */
    public boolean hasNoDiagonal() { return graph.hasNoDiagonal(); }

    /**
     * Returns the local row index for given global row index, returns -1 if no local row for this global row.
     */
    public int LRID(int GRID) { return graph.LRID(GRID); }

    /**
     * Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
     */
    public int GRID(int LRID) { return graph.GRID(LRID); }

    public boolean isFilled() { return graph.isFilled(); }
    
    protected boolean isAllocated() { return isAllocated; }

    public String toMatrixMarketString() {
     StringBuffer buffer  = new StringBuffer(0);
	buffer.append(getNumGlobalRows() + " " +getNumGlobalCols() + " "+ getNumGlobalNonzeros() + "\n");

	int numProcessRows = getNumProcessRows();
	int maxNumIndices = getMaxNumEntries();
	int [] indices = new int [maxNumIndices];
	double [] values = new double [maxNumIndices];
	int []numIndices = new int [1];
	for(int i = 0; i<numProcessRows; i++) {
	    int row = GRID(i);
	    extractGlobalRowCopy(row, maxNumIndices, numIndices, values, indices);
	    
	    for(int j=0; j<numIndices[0]; j++) {
		buffer.append(row+" "+indices[j]+" "+values[j] + "\n");
	    }
	}
     String matrixMarketString = buffer.toString();
	return(matrixMarketString);
    }
    public void print() {
	System.out.println("Number of Global Rows:      "+getNumGlobalRows());
	System.out.println("Number of Global Cols:      "+getNumGlobalCols());
	System.out.println("Number of Global Diagonals: "+getNumGlobalNonzeros());
	if(isLowerTriangular()) System.out.println(" ** Matrix is Lower Triangular **");
	if(isUpperTriangular()) System.out.println(" ** Matrix is Upper Triangular **");
	if(hasNoDiagonal())     System.out.println(" ** Matrix has No Diagonal     **");

	int numProcessRows = getNumProcessRows();
	int maxNumIndices = getMaxNumEntries();
	int [] indices = new int [maxNumIndices];
	double [] values = new double [maxNumIndices];
	int []numIndices = new int [1];

	System.out.println("   Process    Row Index    Col Index    Value");
	for(int i = 0; i<numProcessRows; i++) {
	    int row = GRID(i);
	    extractGlobalRowCopy(row, maxNumIndices, numIndices, values, indices);
	    
	    for(int j=0; j<numIndices[0]; j++) {
		System.out.println("   1           "+row+"        "+indices[j]+"        "+values[j]);
	    }
	}
    }
}
