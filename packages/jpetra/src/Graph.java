/*
 * Graph.java
 *
 * Created on June 28, 2001, 11:01 AM
 */

package Jpetra;

/**
 *
 * @author  Michael William Boldt
 * @version 
 */
public class Graph extends JpetraObject {

    private final BlockMap rowMap;
    private final BlockMap colMap;
    
    private BlockMap domainMap = null;
    private BlockMap rangeMap = null;
    private BlockMap importMap = null;
    private Importer importer = null;
    private BlockMap exportMap = null;
    private Exporter exporter = null;
    
    private boolean filled = false;
    private boolean allocated = false;
    private boolean sorted = false;
    private boolean storageOptimized = false;
    private boolean noRedundancies = false;
    private boolean indicesAreGlobal = false;
    private boolean indicesAreLocal = false;
    private boolean indicesAreContiguous = false;
    private boolean lowerTriangular = true;
    private boolean upperTriangular = true;
    private boolean noDiagonal = true;
    private boolean globalConstantsComputed = false;
    
    private int indexBase = 0;
    private int numGlobalEntries = 0;
    private int numGlobalBlockRows;
    private int numGlobalBlockCols;
    private int numGlobalBlockDiagonals = 0;
    private int numProcessEntries = 0;
    private int numProcessBlockRows;
    private int numProcessBlockCols;
    private int numProcessBlockDiagonals = 0;
    
    private int maxRowDim;
    private int maxColDim;
    private int globalMaxRowDim;
    private int globalMaxColDim;
    private int maxNumNonzeros = 0;
    private int globalMaxNumNonzeros = 0;
    
    private int numGlobalNonzeros = 0;
    private int numGlobalRows;
    private int numGlobalCols;
    private int numGlobalDiagonals = 0;
    private int numProcessNonzeros = 0;
    private int numProcessRows;
    private int numProcessCols;
    private int numProcessDiagonals = 0;
    
    private int [][] indices = null;
    private int [] numAllocatedIndicesPerRow = null;
    private int [] numIndicesPerRow = null;
    private int maxNumIndices = 0;
    private int globalMaxNumIndices = 0;
    private int [] allIndices = null;
    private int lenImporters = 0;
    private int lenExporters = 0;
    private int [] imports = null;
    private int [] exports = null;
    private String copyView = null;
    
    // #ifdef PETRA_LEVELSCHEDULING
    
    /** Creates new Graph */
    public Graph(String copyView, BlockMap rowMap, int [] numIndicesPerRow)
    // throws JpetraException
    {
        if(copyView != "Copy" && copyView != "View") {
            System.out.println("copyView must be either \"Copy\" or \"View\"");
            System.exit(1);
        }
        
        this.rowMap = rowMap;
        colMap = rowMap;
        this.copyView = copyView;
        
        initializeDefaults();

        if(allocate(numIndicesPerRow, 1) != 0) {
            System.out.println("allocate returned non-zero value");
            System.exit(1);
        }
    }
    
    public Graph(String copyView, BlockMap rowMap, BlockMap colMap, int [] numIndicesPerRow)
    //throws JpetraException
    {
        //if(copyView != "Copy" && copyView != "View")
        //    throw new JpetraException("copyView must be either \"Copy\" or \"View\"");
        if(copyView != "Copy" && copyView != "View") {
            System.out.println("copyView must be either \"Copy\" or \"View\"");
            System.exit(1);
        }
        
        this.rowMap = rowMap;
        this.colMap = colMap;
        this.copyView = copyView;
        
        initializeDefaults();
        //if (allocate(numIndicesPerRow, 1) != 0)
        //    throw new JpetraException("allocate returned non-zero value");
        if(allocate(numIndicesPerRow, 1) != 0) {
            System.out.println("allocate returned non-zero value");
            System.exit(1);
        }
    }
    
    public Graph(Graph graph)
    //throws JpetraException
    {
        rowMap = graph.rowMap;
        colMap = graph.colMap;
        domainMap = graph.domainMap;
        rangeMap = graph.rangeMap;
        importMap = graph.importMap;
        importer = graph.importer;
        exportMap = graph.exportMap;
        exporter = graph.exporter;
        filled = graph.filled;
        sorted = graph.sorted;
        noRedundancies = graph.noRedundancies;
        indicesAreGlobal = graph.indicesAreGlobal;
        indicesAreLocal = graph.indicesAreLocal;
        indicesAreContiguous = graph.indicesAreContiguous;
        lowerTriangular = graph.lowerTriangular;
        upperTriangular = graph.upperTriangular;
        noDiagonal = graph.noDiagonal;
        indexBase = graph.indexBase;
        numGlobalEntries = graph.numGlobalEntries;
        numGlobalBlockRows = graph.numGlobalBlockRows;
        numGlobalBlockCols = graph.numGlobalBlockCols;
        numGlobalBlockDiagonals = graph.numGlobalBlockDiagonals;
        numProcessEntries = graph.numProcessEntries;
        numProcessBlockRows = graph.numProcessBlockRows;
        numProcessBlockCols = graph.numProcessBlockCols;
        numProcessBlockCols = graph.numProcessBlockDiagonals;
        maxRowDim = graph.maxRowDim;
        maxColDim = graph.maxColDim;
        globalMaxRowDim = graph.globalMaxRowDim;
        globalMaxColDim = graph.globalMaxColDim;
        maxNumNonzeros = graph.maxNumNonzeros;
        globalMaxNumNonzeros = graph.globalMaxNumNonzeros;
        numGlobalNonzeros = graph.numGlobalNonzeros;
        numGlobalRows = graph.numGlobalRows;
        numGlobalCols = graph.numGlobalCols;
        numGlobalDiagonals = graph.numGlobalDiagonals;
        numProcessNonzeros = graph.numProcessNonzeros;
        numProcessRows = graph.numProcessRows;
        numProcessCols = graph.numProcessCols;
        numProcessDiagonals = graph.numProcessDiagonals;
        copyView = "Copy";
        
        //if(allocate(graph.numIndicesPerRow, 1) != 0)
        //    throw new JpetraException("allocate returned non-zero value");
        if(allocate(graph.numIndicesPerRow, 1) != 0) {
            System.out.println("allocate returned non-zero value");
            System.exit(1);
        }
        for(int i=0; i<numProcessBlockRows; i++) {
            numIndicesPerRow[i] = numAllocatedIndicesPerRow[i];
	    indices[i] = new int [numIndicesPerRow[i]];
            System.arraycopy(graph.indices[i], 0, indices[i], 0, numIndicesPerRow[i]);
        }
        
        maxNumIndices = graph.maxNumIndices;
        if(importMap != null) importMap = new BlockMap(graph.importMap);
        if(importer != null) importer = new Importer(graph.importer);
        if(exportMap != null) exportMap = new BlockMap(graph.exportMap);
        if(exporter != null) exporter = new Exporter(graph.exporter);
    }
    
    private void initializeDefaults() {
        numGlobalBlockRows = rowMap.getNumGlobalElements();
        numGlobalBlockCols = colMap.getNumGlobalElements();
        numProcessBlockRows = rowMap.getNumProcessElements();
        numProcessBlockCols = colMap.getNumProcessElements();
        
        numGlobalRows = rowMap.getNumGlobalEquations();
        numGlobalCols = colMap.getNumGlobalEquations();
        numProcessRows = rowMap.getNumProcessEquations();
        numProcessCols = colMap.getNumProcessEquations();
        
        globalMaxRowDim = rowMap.getMaxElementSize();
        maxRowDim = rowMap.getMaxElementSize();
        globalMaxColDim = colMap.getMaxElementSize();
        maxColDim = colMap.getMaxElementSize();
        
    }
    
    private int allocate(int [] numIndicesPerRow, int inc) {
        int i, j;
        
        indices = new int [numProcessBlockRows] [];
        this.numIndicesPerRow = new int [numProcessBlockRows];
        numAllocatedIndicesPerRow = new int [numProcessBlockRows];
        
        if(copyView == "Copy")
            for(i=0; i<numProcessBlockRows; i++) {
                final int numIndices = numIndicesPerRow[i*inc];
                if(numIndices > 0) indices[i] = new int [numIndices];
                else indices[i] = null;
                
                numAllocatedIndicesPerRow[i] = numIndices;
                this.numIndicesPerRow[i] = 0;
                for(j=0; j<numIndices; j++) indices[i][j] = indexBase - 1;
            }
        else
            for(i=0; i<numProcessBlockRows; i++) {
                indices[i] = null;
                numAllocatedIndicesPerRow[i] = 0;
                this.numIndicesPerRow[i] = 0;
            }
        // #ifdef PETRA_LEVELSCHEDULING
        
        allocated = true;
        return 0;
    }
    
    private int reAllocate() {
        numAllocatedIndicesPerRow = new int [numProcessBlockRows];
        
        for(int i=0; i<numProcessBlockRows; i++) numAllocatedIndicesPerRow[i] = numIndicesPerRow[i];
        storageOptimized = false;
        return 0;
    }
    
    public int insertGlobalIndices(int row, int numIndices, int [] indices) {
        if(indicesAreLocal) return -2;
        if(indicesAreContiguous) return -3;
        indicesAreGlobal = true;
        row = LRID(row);
        
        return(insertIndices(row, numIndices, indices));
    }
    
    public int insertProcessIndices(int row, int numIndices, int [] indices) {
        if(indicesAreGlobal) return -2;
        if(indicesAreContiguous) return -3;
        indicesAreLocal = true;
        return(insertIndices(row, numIndices, indices));
    }
    
    protected int insertIndices(int row, int numIndices, int [] indices) {
        sorted = false;
        
        int j;
        int [] tmpIndices;
        int ierr = 0;
        
        if(row < 0 || row >= numProcessBlockRows) return -1;
        if(copyView == "View") {
            if(this.indices[row] != null) ierr = 2;
            this.indices[row] = indices;
            numAllocatedIndicesPerRow[row] = numIndices;
            numIndicesPerRow[row] = numIndices;
        }
        else {
            int start = numIndicesPerRow[row];
            int stop = start + numIndices;
            int numAllocatedIndices = numAllocatedIndicesPerRow[row];
            if(stop > numAllocatedIndices) {
                if(numAllocatedIndices == 0) this.indices[row] = new int [numIndices];
                else {
                    ierr = 1;
                    tmpIndices = new int [stop];
                    System.arraycopy(this.indices[row], 0, tmpIndices, 0, start);
                    this.indices[row] = tmpIndices;
                }
                numAllocatedIndicesPerRow[row] = stop;
            }
            
            numIndicesPerRow[row] = stop;
            for(j=start; j<stop; j++) this.indices[row][j] = indices[j-start];
        }
        maxNumIndices = Math.max(maxNumIndices, numIndicesPerRow[row]);
        return ierr;
    }

    public int removeGlobalIndices(int row, int numIndices, int [] indices) {
	int j, k;
	int ierr = 0;
	int [] loc = new int [1];

	if(indicesAreLocal) return -2;
	if(copyView == "View") return -3;

	row = LRID(row);
	if(row < 0 || row >= numProcessBlockRows) return -1;

	int numCurrentIndices = numIndicesPerRow[row];
	
	for(j=0; j<numIndices; j++){
	    int index = indices[j];
	    if(findGlobalIndexLoc(row, index, j, loc)) {
		for(k=loc[0]+1; k<numCurrentIndices; k++) this.indices[row][k-1] = this.indices[row][k];
		numCurrentIndices--;
		numIndicesPerRow[row]--;
	    }
	}
	return ierr;
    }

    public int removeProcessIndices(int row, int numIndices, int [] indices) {
	if(indicesAreGlobal) return -2;

	int j, k;
	int ierr = 0;
	int [] loc = new int [1];

	if(copyView == "View") return -3;

	if(row < 0 || row >= numProcessBlockRows) return -1;

	int numCurrentIndices = numIndicesPerRow[row];
	for(j=0; j<numIndices; j++) {
	    int index = indices[j];
	    if(findProcessIndexLoc(row, index, j, loc)) {
		for(k=loc[0]+1; k<numCurrentIndices; k++) this.indices[row][k-1] = this.indices[row][k];
		numCurrentIndices--;
		numIndicesPerRow[row]--;
	    }
	}
	return ierr;
    }

    public int removeGlobalIndices(int row) {
	int j;
	int ierr = 0;

	if(indicesAreLocal) return -2;
	if(copyView == "View") return -3;

	row = LRID(row);

	if(row < 0 || row >= numProcessBlockRows) return -1;

	int numIndices = numIndicesPerRow[row];
	numIndicesPerRow[row] = 0;

	for(j=0; j<numIndices; j++) indices[row][j] = indexBase - 1;

	return ierr;
    }

    public int removeProcessIndices(int row) {
	int j;
	int ierr = 0;
       
	if(indicesAreGlobal) return -2;
	if(copyView == "View") return -3;

	int numIndices = numIndicesPerRow[row];
	numIndicesPerRow[row] = 0;
	
	for(j=0; j<numIndices; j++) indices[row][j] = -1;

	return ierr;
    }

    

    protected boolean findGlobalIndexLoc(int localRow, int index, int start, int [] loc) {
	int j;
	int numIndices = numIndicesPerRow[localRow];
	int [] indices = this.indices[localRow];
	
	if(indicesAreLocal) index = LCID(index);

	for(j=0; j<numIndices; j++) {
	    int j0 = (j + start) % numIndices;
	    if(indices[j0] == index) {
		loc[0] = j0;
		return true;
	    }
	}
	return false;
    }

    protected boolean findProcessIndexLoc(int localRow, int index, int start, int [] loc) {
	int j;
	int numIndices = numIndicesPerRow[localRow];

	if(indicesAreGlobal) return false;
	for(j=0;j<numIndices; j++) {
	    int j0 = (j + start) % numIndices;
	    if(indices[localRow][j0] == index) {
		loc[0] = j0;
		return true;
	    }
	}
	return false;
    }

    public int transformToLocal() {
	return(transformToLocal(rowMap, colMap));
    }

    public int transformToLocal(BlockMap domain, BlockMap range) {
	domainMap = domain;
	rangeMap = range;

	makeIndicesLocal(domainMap, rangeMap);

	sortIndices();

	removeRedundantIndices();

	computeGlobalConstants();

	filled = true;

	return 0;
    }

    private int computeGlobalConstants() {
	int i, j, k;

	if(globalConstantsComputed) return 0;

	int [] tempVec = new int [3];
	int [] tempRes = new int [3];
	int [] tempMaxNumI = new int [1];

	numProcessEntries = 0;
	maxNumIndices = 0;
	for(i=0; i<numProcessBlockRows; i++) {
	    numProcessEntries += numIndicesPerRow[i];
	    maxNumIndices = Math.max(maxNumIndices, numIndicesPerRow[i]);


	}

	// Case 1: Constant block size (including blocksize = 1)
	if(rowMap.hasConstantElementSize() && colMap.hasConstantElementSize()) {
	    tempVec[0] = numProcessEntries;
	    tempVec[1] = numProcessBlockDiagonals;
	    tempMaxNumI[0] = maxNumIndices;

	    tempRes = getComm().sumAll(tempVec);
	    numGlobalEntries = tempRes[0];
	    globalMaxNumIndices = tempRes[1];
	    tempRes = getComm().maxAll(tempMaxNumI);
	    globalMaxNumIndices = tempRes[0];

	    int rowElementSize = rowMap.getMaxElementSize();
	    int colElementSize = colMap.getMaxElementSize();
	    numProcessNonzeros = numProcessEntries * rowElementSize * colElementSize;
	    numGlobalNonzeros = numGlobalEntries * rowElementSize * colElementSize;
	    maxNumNonzeros = numGlobalEntries * rowElementSize * colElementSize;
	    globalMaxNumNonzeros = globalMaxNumIndices * rowElementSize * colElementSize;
	}
	// Case 2: Variable block size (more work)
	else {
	    numProcessNonzeros = 0;
	    maxNumNonzeros = 0;
	    int [] rowElementSizeList = rowMap.getElementSizeList();
	    int [] colElementSizeList = rowElementSizeList;
	    if(importer != null) colElementSizeList = importMap.getElementSizeList();
	    for(i=0; i<numProcessBlockRows; i++) {
		int numEntries = numIndicesPerRow[i];
		int [] indices = this.indices[i];
		if(numEntries > 0) {
		    int curNumNonzeros = 0;
		    int rowDim = rowElementSizeList[i];
		    for(j=0; j<numEntries; j++) {
			int colDim = colElementSizeList[indices[j]];
			curNumNonzeros += rowDim * colDim;
			maxColDim = Math.max(maxColDim, colDim);
		    }
		    maxNumNonzeros = Math.max(maxNumNonzeros, curNumNonzeros);
		    numProcessNonzeros += curNumNonzeros;
		}
	    }

	    // Sum up all nonzeros
	    
	    tempVec[0] = numProcessEntries;
	    tempVec[1] = numProcessBlockDiagonals;
	    tempVec[2] = numProcessNonzeros;
	    
	    tempRes = getComm().sumAll(tempVec);

	    numGlobalEntries = tempRes[0];
	    numGlobalBlockDiagonals = tempRes[1];
	    numGlobalNonzeros = tempRes[2];

	    tempVec[0] = maxNumIndices;
	    tempVec[1] = maxNumNonzeros;
	    
	    tempRes = getComm().sumAll(tempVec);

	    globalMaxNumIndices = tempRes[0];
	    globalMaxNumNonzeros = tempRes[1];
	}
	globalConstantsComputed = true;

	return 0;
    }

    public int sortIndices() {
	if(!indicesAreLocal) return -1;
	if(sorted) return 0;

	// For each row, sort column entries from smallest to largest.
	// Use shell sort, which is fast if indices are already sorted.

	for(int i=0; i<numProcessBlockRows; i++) {
	    int n = numIndicesPerRow[i];
	    final int [] list = indices[i];
	    int m = n / 2;

	    while(m > 0) {
		int max = n - m;
		for(int j=0; j<max; j++)
		    for(int k=j; k>=0; k-=m) {
			if(list[k+m] >= list[k])
			    break;
			int iTemp = list[k+m];
			list[k+m] = list[k];
			list[k] = iTemp;
		    }
		m = m / 2;
	    }
	}
	sorted = true;
	return 0;
    }

    public int removeRedundantIndices() {
	int i, j, k, jj;

	if(noRedundancies) return 0;
	if(!sorted) return -1; // Must have sorted index set
	if(indicesAreGlobal) return -2; // Indices must be local

	// For each row, remove column indices that are repeated.
	// Also, determine if graph is upper or lower triangular or has no diagonal.
	// Note: This function assums that sortIndices was already called.

	this.numProcessDiagonals = 0;
	this.numProcessBlockDiagonals = 0;
	for(i=0; i<this.numProcessBlockRows; i++) {
	    boolean diagFound = false;
	    int numIndices = this.numIndicesPerRow[i];
	    if(numIndices > 0) {
		final int [] indices = this.indices[i];
		int j0 = 0;
		jj = indices[0];
		if(jj > i) this.lowerTriangular = false;
		if(jj < i) this.upperTriangular = false;
		if(jj ==i) diagFound = true;
		for(j=i; j<numIndices; j++) {
		    jj = indices[j];
		    if(jj > i) this.lowerTriangular = false;
		    if(jj < i) this.upperTriangular = false;
		    if(jj ==i) diagFound = true;
		    if(jj == indices[j0]) { // Check if index is repeated
			this.numIndicesPerRow[i]--;
			this.numProcessNonzeros--;
			for(k=j; k<numIndices-1; k++) indices[k] = indices[k+1];
			numIndices--;
		    }
		    else j0=j;
		}
		if(diagFound) {
		    this.numProcessBlockDiagonals++;
		    this.numProcessDiagonals += rowMap.getElementSize(i);
		}
	    }
	}
	this.noDiagonal = (numProcessBlockDiagonals == 0);
	this.noRedundancies = true;
	return 0;
    }

    protected int makeIndicesLocal(BlockMap domain, BlockMap range) {
	int i, j, k;

	if(!indicesAreLocal) {
	    domainMap = domain;
	    rangeMap = range;

	    // For each row, check if column indices are owned or not.
	    // If owned, transform global index to local index.
	    // If not owned, add to importMap for later use.

	    numProcessBlockCols = domainMap.getNumProcessElements();

	    int incBlockCols = Math.max(Math.min(numProcessBlockCols/4, 100), 10);
	    int maxBlockCols = 0;
	    int [] colIndices = null;
	    if(domainMap.isDistributedGlobal()) {
		maxBlockCols = numProcessBlockCols;
		colIndices = new int [maxBlockCols];
	    }
	    int [] newColIndices = null;

	    int newNumProcessBlockCols = numProcessBlockCols;
	    for(i=0; i<numProcessBlockRows; i++) {
		final int numIndices = numIndicesPerRow[i];
		for(j=0; j<numIndices; j++) {
		    int GID = indices[i][j];
		    int LID = LCID(GID);
		    if(LID == -1) {
			boolean extNew = true;
			for(k=numProcessBlockCols; k<newNumProcessBlockCols; k++)
			    if(colIndices[k] == GID) {
				extNew = false;
				break;
			    }
			if(extNew) {
			    if(newNumProcessBlockCols >= maxBlockCols) { // Need to expand...
				maxBlockCols = Math.max(numProcessBlockCols+incBlockCols, maxBlockCols+incBlockCols);
				newColIndices = new int [maxBlockCols];
				System.arraycopy(colIndices, numProcessBlockCols, newColIndices, numProcessBlockCols, newNumProcessBlockCols-numProcessBlockCols);
				colIndices = newColIndices;
			    }
			    colIndices[newNumProcessBlockCols] = GID;
			    newNumProcessBlockCols++;
			}
		    }
		}
	    }
	    // Create importMap.  This map will be used to facilitate communication in matrix classes.
	    if(domainMap.isDistributedGlobal()) {
		// Find processes that own the off-procssor GIDs
		int numRemote = newNumProcessBlockCols - numProcessBlockCols;
		// int [] remoteColIndices = colIndices + numProcessBlockCols;
		int [] remoteColIndices = new int [colIndices.length-numProcessBlockCols];
		System.arraycopy(colIndices, numProcessBlockCols, remoteColIndices, 0, colIndices.length-numProcessBlockCols);
		int nLists = 1;
		int [] PIDList = null;
		int [] sizeList = null;
 		int [] remoteSizeList = null;
		boolean doSizes = !domainMap.hasConstantElementSize();

		if(numRemote > 0) PIDList = new int [numRemote];

		if(doSizes) {
		    if(newNumProcessBlockCols > 0) sizeList = new int [newNumProcessBlockCols];
		    // remoteSizeList = sizeList + numProcessBlockCols;
		    remoteSizeList = new int [sizeList.length-numProcessBlockCols];
		    System.arraycopy(sizeList, numProcessBlockCols, remoteSizeList, 0, remoteSizeList.length);
		    nLists++;
		}
		domainMap.getRemoteIDList(numRemote, remoteColIndices, PIDList, null, remoteSizeList);
		System.arraycopy(remoteSizeList, 0, sizeList, numProcessBlockCols, remoteSizeList.length);

		// Sort external column indices so that all columns coming from a given remote process are contiguous

		Util util = new Util();
		int [][] sortLists = new int [2][];
		sortLists[0] = remoteColIndices;
		sortLists[1] = remoteSizeList;
		util.sort(true, numRemote, PIDList, 0, null, nLists, sortLists);
		System.arraycopy(remoteSizeList, 0, sizeList, numProcessBlockCols, remoteSizeList.length);

		domainMap.getGlobalElements(colIndices);
		if(doSizes) domainMap.getElementSizeList(sizeList);

		numProcessBlockCols = newNumProcessBlockCols;

		// Make import map with same element sizes as domain map
		if(domainMap.hasConstantElementSize())
		    importMap = new BlockMap(-1, newNumProcessBlockCols, colIndices,
		    domainMap.getMaxElementSize(), domainMap.getIndexBase(), domainMap.getComm());

		else // Most general case where block size is variable
		    importMap = new BlockMap(-1, newNumProcessBlockCols, colIndices, sizeList, 
                    domainMap.getIndexBase(), domainMap.getComm());

		// Define an export map only if rowMap and rangeMap are different

		if(!rowMap.sameAs(rangeMap)) {
		    exporter = new Exporter(rowMap, rangeMap);
		    exportMap = rangeMap;
		}

		indicesAreLocal = true;
		indicesAreGlobal = false;

		// Transform indices to local index space

		for(i=0; i<numProcessBlockRows; i++) {
		    final int numIndices = numIndicesPerRow[i];
		    for(j=0; j<numIndices; j++) {
			int GID = indices[i][j];
			int LID = LCID(GID);
			if(LID != -1) indices[i][j] = LID;
			else {
			    System.out.println("Internal error in transformToLocal");
			    System.exit(1);
			}
		    }
		}
	    }
	    // Compute number of local columns

	    // Easy to do if constant element size
	    if(rangeMap.hasConstantElementSize())
		numProcessCols = numProcessBlockCols * rangeMap.getMaxProcessElementSize();
	    else {
		numProcessCols = 0;
		for(i=0; i<numProcessBlockRows; i++) {
		    final int numIndices = numIndicesPerRow[i];
		    for(j=0; j<numIndices; j++) {
			int elementSize = rangeMap.getElementSize(indices[i][j]);
			numProcessCols += elementSize;
		    }
		}
	    }
	}
	else {
	    // Do a sanity check on column indices.  They must all be in the range 0 to numProcessBlockCols.
	    // Indices will be sorted so we only need to check the last onw.

	    if(!sorted) sortIndices();

	    for(i=0; i<numProcessBlockRows; i++) {
		int numIndices = numIndicesPerRow[i];
		if(numIndices > 0)
		    if(indices[i][numIndices-1] >= numProcessBlockCols) return -1;
	    }
	}
	return 0;
    }

    // public int optimizeStorage() { } IS THIS POSSIBLE W/JAVA?

    public int extractGlobalRowCopy(int row, int lenOfIndices, int [] numIndices, int [] indices) {
	int j;

	row = LRID(row);
	if(row < 0 || row >= numProcessBlockRows) return -1;

	numIndices[0] = numIndicesPerRow[row];
	if(lenOfIndices < numIndices[0]) return -2;

	if(indicesAreLocal) for(j=0; j<numIndices[0]; j++) indices[j] = GCID(this.indices[row][j]);
	else System.arraycopy(this.indices[row], 0, indices, 0, numIndices[0]);

	return 0;
    }

    public int extractProcessRowCopy(int row, int lenOfIndices, int [] numIndices, int [] indices) {
	int j;

	if(row < 0 || row >= numProcessBlockRows) return -1;

	numIndices[0] = numIndicesPerRow[row];
	if(lenOfIndices < numIndices[0]) return -2;

	if(indicesAreGlobal) return -3;

	System.arraycopy(this.indices[row], 0, indices, 0, numIndices[0]);

	return 0;
    }

    public int extractGlobalRowView(int row, int [] numIndices, int [] view) {
	row = LRID(row);

	if(row < 0 || row >= numProcessBlockRows) return -1;
	if(indicesAreLocal) return -2;

	numIndices[0] = numIndicesPerRow[row];

	view = indices[row];

	return 0;
    }

    public int extractProcessRowView(int row, int [] numIndices, int [] view) {
	if(row < 0 || row >= numProcessBlockRows) return -1;
	if(indicesAreGlobal) return -2;

	numIndices[0]= numIndicesPerRow[row];

	view = indices[row];

	return 0;
    }

    public boolean isFilled() { return filled; }

    public boolean isSorted() { return sorted; }

    public boolean hasOptimizedStorage() { return storageOptimized; }

    public boolean indicesAreGlobal() { return indicesAreGlobal; }

    public boolean indicesAreLocal() { return indicesAreLocal; }

    public boolean isLowerTriangular() { return lowerTriangular; }

    public boolean isUpperTriangular() { return upperTriangular; }

    public boolean hasNoDiagonal() { return noDiagonal; }

    public boolean isProcessGlobalRow(int GID) { return rowMap.isMyGID(GID); }

    public int getNumProcessRows() { return numProcessRows; }
    
    public int getNumGlobalRows() { return numGlobalRows; }

    public int getNumProcessCols() { return numProcessCols; }

    public int getNumGlobalCols() { return numGlobalCols; }

    public int getNumGlobalNonzeros() { return numGlobalNonzeros; }

    public int getNumGlobalDiagonals() { return numGlobalDiagonals; }

    public int getNumProcessDiagonals() { return numProcessDiagonals; }

    public int getNumProcessBlockRows() { return numProcessBlockRows; }

    public int getNumGlobalBlockRows() { return numGlobalBlockRows; }

    public int getNumProcessBlockCols() { return numProcessBlockCols; }

    public int getNumProcessBlockDiagonals() { return numProcessBlockDiagonals; }

    public int getNumGlobalBlockDiagonals() { return numGlobalBlockDiagonals; }

    public int getNumGlobalEntries() { return numGlobalEntries; }

    public int getNumProcessEntries() { return numProcessEntries; }

    public int getMaxRowDim() { return maxRowDim; }

    public int getGlobalMaxRowDim() { return globalMaxRowDim; }

    public int getMaxColDim() { return maxColDim; }

    public int getGlobalMaxColDim() { return globalMaxColDim; }

    /**
     * Gets the local row index for a given global row index. Returns -1
     * if no local row for this global row.
     */
    public int LRID(int GRID) { return rowMap.getLID(GRID); }

    /**
     * Gets the global row index for given local row index. Returns indexBase-1
     * if we don't hae this local row.
     */
    public int GRID(int LRID) { return rowMap.getGID(LRID); }

    /**
     * Get the local column index for a given global column index. 
     * Returns -1 if no local row for this global row.
     */
    public int LCID(int GCID) {
	int index = domainMap.getLID(GCID);
	if(index != -1) return index;
	if(importMap == null) return -1;
	return importMap.getLID(GCID);
    }
    
    /**
     * Get the global column  index for a given local column index.
     * Returns indexBase-1 if no global row for this local row.
     */
    public int GCID(int LCID) {
	int index = domainMap.getGID(LCID); // Check row map first
	if(index != indexBase-1) return index;
	if(importMap == null) return indexBase-1;
	return importMap.getGID(LCID); // Check col map
    }

    /**
     * Returns true if given GRID belongs to the calling process in this map,
     * otherwise returns false.
     */
    public boolean isMyGRID(int GRID) { return (LRID(GRID) != -1); }

    /**
     * Returns true if given LRID belongs to the calling process in this map,
     * otherwise returns false.
     */
    public boolean isMyLRID(int LRID) { return (GRID(LRID) != indexBase-1); }

    /**
     * Returns true if given GCID belongs to the calling process in this map,
     * otherwise returns false.
     */
    public boolean isMyGCID(int GCID) { return (LCID(GCID) != -1); }

    /**
     * Returns true if given LCID belongs to the calling process in this map,
     * otherwise returns false.
     */
    public boolean isMyLCID(int LCID) { return (GCID(LCID) != indexBase-1); }

    public int getNumProcessNonzeros() { return numProcessNonzeros; }

    public int getNumGlobalIndices(int row) {
	row = LRID(row);
	if(row != -1) return numIndicesPerRow[row];
	else return 0;
    }

    public int getNumAllocatedGlobalIndices(int row) {
	row = LRID(row);
	if(row != -1) return numAllocatedIndicesPerRow[row];
	else return 0;
    }

    public boolean hasNoRedundancies() { return noRedundancies; }

    public int getMaxNumIndices() { return maxNumIndices; }

    public int getGlobalMaxNumIndices() { return globalMaxNumIndices; }

    public int getMaxNumNonzeros() { return maxNumNonzeros; }

    public int getGlobalMaxNumNonzeros() { return globalMaxNumNonzeros; }

    public int getNumProcessIndices(int row) { return numIndicesPerRow[row]; }

    public int getNumAllocatedProcessIndices(int row) { return numAllocatedIndicesPerRow[row]; }

    public int getIndexBase() { return indexBase; }

    public BlockMap getRowMap() { return rowMap; }

    public BlockMap getColMap() { return colMap; }

    public BlockMap getDomainMap() { return domainMap; }

    public BlockMap getRangeMap() { return rangeMap; }

    public Comm getComm() { return rowMap.getComm(); }

    public BlockMap getImporterMap() { 
	if(importMap == null) return rowMap; 
	else return importMap;
    }

    public Importer getImporter() { return importer; }

    public BlockMap getExporterMap() { return exportMap; }

    public Exporter getExporter() { return exporter; }

    public int doImporter(Graph sourceGraph, Importer importer, String combineMode) {
	if(!rowMap.sameAs(importer.getTargetMap())) return -2;
	if(!sourceGraph.rowMap.sameAs(importer.getSourceMap())) return -3;

	int numSameIDs = importer.getNumSameIDs();
	int numPermuteIDs = importer.getNumPermuteIDs();
	int numRemoteIDs = importer.getNumRemoteIDs();
	int numExporterIDs = importer.getNumExporterIDs();
	int [] exportLIDs = importer.getExporterLIDs();
	int [] remoteLIDs = importer.getRemoteLIDs();
	int [] permuteToLIDs = importer.getPermuteToLIDs();
	int [] permuteFromLIDs = importer.getPermuteFromLIDs();
	int sizeOfPacket = sourceGraph.getGlobalMaxNumIndices() + 2;
	int nSend = sizeOfPacket * importer.getNumSend();
	int nRecv = sizeOfPacket * importer.getNumRecv();

	return doTransfer(sourceGraph, combineMode, numSameIDs, numPermuteIDs, numRemoteIDs, numExporterIDs,
			  permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs, nSend, nRecv,
			  sizeOfPacket, lenExporters, exports, lenImporters, imports, false);
    }

    // public int doExporter(Graph source, Exporter exporter, String combineMode) { }

    // public int doImporter(Graph source, Exporter exporter, String combineMode) { }

    // public int doExporter(Graph source, Importer importer, String combineMode) { }

    private int doTransfer(Graph sourceGraph, String combineMode, int numSameIDs, int numPermuteIDs,
			   int numRemoteIDs, int numExporterIDs, int [] permuteToLIDs, int [] permuteFromLIDs,
			   int [] remoteLIDs, int [] exportLIDs, int nSend, int nRecv, int sizeOfPacket,
			   int lenExporters, int [] exports, int lenImporters, int [] imports, boolean doReverse) {
	int ierr = 0;
	ierr = copyAndPermute(this, sourceGraph, numSameIDs, numPermuteIDs, permuteToLIDs, permuteFromLIDs);

	if(ierr != 0) return ierr;
	if(combineMode == "Zero") return 0;

	// Need to send global row ID, numIndices for the row, then the indices.
	if(nSend > lenExporters) {
	    exports = new int [nSend];
	    lenExporters = nSend;
	}

	ierr = pack(sourceGraph, numExporterIDs, exportLIDs, exports);

	if(ierr != 0) return ierr;

	if(nRecv > lenImporters) {
	    imports = new int [nRecv];
	    lenImporters = nRecv;
	}

	// #ifdef PETRA_MPI...

	return ierr;
    }

    private int copyAndPermute(Graph target, Graph source, int numSameIDs, int numPermuteIDs,
			       int [] permuteToLIDs, int [] permuteFromLIDs) {
	int i;
	int ierr;
	int row;
	int [] numIndices = new int[1];
	int [] indices = new int [source.globalMaxNumIndices];
	int fromRow, toRow;

	// Do copy first
	if(numSameIDs > 0) {
	    if(source.indicesAreLocal) {
		int maxNumIndices = source.maxNumIndices;
		indices = new int [maxNumIndices]; // Need some temporary space

		for(i=0; i<numSameIDs; i++) {
		    row = target.GRID(i);
		    ierr = source.extractGlobalRowCopy(row, maxNumIndices, numIndices, indices);
		    if(ierr < 0) {
			System.out.println("error in extractGlobalRowCopy (copyAndPermute) "+ierr);
			System.exit(1);
		    }
		    ierr = target.insertGlobalIndices(row, numIndices[0], indices);
		    if(ierr < 0) {
			System.out.println("error in insertGlobalIndices (copyAndPermute)< "+ierr);
			System.exit(2);
		    }
		}
	    }
		else { // source.indicesAreGlobal
		    for(i=0; i<numSameIDs; i++) {
			row = target.GRID(i);
			if(source.extractGlobalRowView(row, numIndices, indices) != 0) {
			    System.out.println("error in Graph.copyAndPermute(3)");
			    System.exit(3);
			}
			if(target.insertGlobalIndices(row, numIndices[0], indices) != 0) {
			    System.out.println("error in Graph.copyAndPermute(4)");
			    System.exit(4);
			}
		    }
		}
	    }

	    // Do local permutation next
	    if(numPermuteIDs > 0) {
		if(source.indicesAreLocal) {
		    int maxNumIndices = source.maxNumIndices;
		    indices = new int [maxNumIndices]; // Need some temporary space

		    for(i=0; i<numPermuteIDs; i++) {
			fromRow = source.GRID(permuteFromLIDs[i]);
			toRow = target.GRID(permuteToLIDs[i]);
			if(source.extractGlobalRowCopy(fromRow, maxNumIndices, numIndices, indices) != 0) {
			     System.out.println("error in Graph.copyAndPermute(5)");
			     System.exit(5);
			}
			if(target.insertGlobalIndices(toRow, numIndices[0], indices) != 0) {
			     System.out.println("error in Graph.copyAndPermute(6)");
			     System.exit(6);
			}
		    }
		}
		else { // source.indicesAreGlobal
		    for(i=0; i<numPermuteIDs; i++) {
			fromRow = source.GRID(permuteFromLIDs[i]);
			toRow = target.GRID(permuteToLIDs[i]);
			if(source.extractGlobalRowView(fromRow, numIndices, indices) != 0) {
			     System.out.println("error in Graph.copyAndPermute(7)");
			     System.exit(7);
			}
			if(target.insertGlobalIndices(toRow, numIndices[0], indices) != 0) {
			     System.out.println("error in Graph.copyAndPermute(8)");
			     System.exit(8);
			}
		    }
		}
	    }
	return 0;
    }

    private int pack(Graph source, int numSendIDs, int [] sendLIDs, int [] sends) {
	int i;
	int [] numIndices = new int [1];
	int [] indices;
	int fromRow;

	if(numSendIDs > 0) {
	    
	    // Each segment of sends will be filled by a packed row of information for each row as follows:
	    // 1st int: GRID of row where GRID is the global row ID for the source graph.
	    // next int: numIndices, number of indices in row.
	    // next numIndices: the actual indices for the row.
	    // any remaining space (of length globalMaxNumIndices - numIndices ints) wil be wasted but
	    // we need fixed sized segments for concurrent communication routines.

	    int globalMaxNumIndices = source.globalMaxNumIndices;
	    indices = new int [globalMaxNumIndices];
	    int inc = globalMaxNumIndices + 2;

	    for(i=0; i<numSendIDs; i++) {
		fromRow = source.GRID(sendLIDs[i]);

		if(source.extractGlobalRowCopy(fromRow, globalMaxNumIndices, numIndices, indices) != 0) {
		    System.out.println("Error in Graph.pack() call to extractGlobalRowCopy returned error");
		    System.exit(1);
		}

		sends[i*inc] = fromRow;
		sends[i*inc+1] = numIndices[0];
		System.arraycopy(indices, 0, sends, i*inc+2, globalMaxNumIndices);
	    }
	}
	return 0;
    }

    protected void setIndicesAreGlobal(boolean flag) {  indicesAreGlobal = flag; }
    
    protected void setIndicesAreLocal(boolean flag) { indicesAreLocal = flag; }
	   
    protected int [][] getIndices() { return indices; }

    protected int [] getNumIndicesPerRow() { return numIndicesPerRow; }

    protected int [] getNumAllocatedIndicesPerRow() { return numAllocatedIndicesPerRow; }

    protected void setSorted(boolean flag) { sorted = flag; }

    protected void setNoRedundancies(boolean flag) { noRedundancies = flag; }

    public int [][] GETINDICES() {return indices;}

    public void DISPLAY() {
	for(int i=0; i<numGlobalRows; i++) {
	    System.out.print("R"+i+": ");
	    for(int j=0; j<numIndicesPerRow[i]; j++)
		System.out.print(indices[i][j]+" ");
	    System.out.println();
	}
    }
}
