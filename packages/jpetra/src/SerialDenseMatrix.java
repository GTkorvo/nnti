package Jpetra;

// Math provides abs(), max(), and min()
import java.lang.Math;
// For matrix multiplication
import org.netlib.blas.Dgemm;
import org.netlib.blas.Dsymm;

/**
 * The <code>SerialDenseMatrix</code> class is intended to provide very basic support for dense rectangular matrices.
 * <p>
 * <b>Constructing <code>SerialDenseMatrix</code> Objects</b>
 * <p>
 * There are four <code>SerialDenseMatrix</code> constructors.  The first constructs a zero-sized object which should be made
 * to appropriate length using the <code>shape()</code> or <code>reshape()</code> methods and then filled with <code>setElement()</code>. 
 * The second constructs an object sized to the dimensions specified, which should be filled with <code>setElement()</code>.
 * The third is a constructor that accepts user data as a 2D array (<code>int[][]</code>), and the fourth is a copy constructor.
 * <p>
 * <b>Extracting Data from <code>SerialDenseMatrix</code> Objects</b>
 * <p>
 * Once an <code>SerialDenseMatrix</code> is constructed, it is possible to view the data via access functions.
 * <p>
 * <b>Vector and Utility Functions</b>
 * <p>
 * Once an SerialDenseMatrix is constructed, several mathematical functions can be applied to
 * the object.  Specifically:
 * <ul>
 *   <li> Multiplication.
 *   <li> Norms.
 * </ul>
 *
 * The SerialDenseMatrix code originated as C++ code in ePetra and was ported to Java.
 * @author Jason Cross
 */ 
 public class SerialDenseMatrix extends JpetraObject {
    /**
     * number of rows in matrixA
     */
    int numRows;
    
    /**
     * number of columns in matrixA
     */
    int numCols;
   
    /**
     * the "Leading Dimension", or stride between vectors in memory
     */
    int lda;

    /**
     * set to <code>true</code> if <code>matrixA</code> is a copy of the source data<br>
     * set to <code>false</code> if <code>matrixA</code> is a reference to user data
     */
    boolean copied;

    /**
     * the data object of the dense matrix,<br>
     * it is 1-dimensional but access is 2-dimensinal
     */    
    double[] matrixA;
    
    
    /**
     * construct a <code>null</code> <code>SerialDenseMatrix</code><br>
     * needs to be shaped with <code>shape()</code> and<br>
     * needs to be filled with <code>setElement()</code>
     */
    public SerialDenseMatrix() {
        this.numRows = 0;
        this.numCols = 0;
        this.lda = 0;
        this.copied = false;
        // initilize matrixA?
    }
    
    /**
     * construct a zero filled <code>numRows</code> X <code>numCols</code> <code>SerialDenseMatrix</code><br>
     * needs to be filled with <code>setElement()</code>
     *
     * @param numRows number of rows in <code>this</code> <code>SerialDenseMatrix</code>
     * @param numCols number of columns in <code>this</code> <code>SerialDenseMatrix</code>
     */    
    public SerialDenseMatrix(int numRows, int numCols) {
        this.numRows = 0;
        this.numCols = 0;
        this.lda = 0;
        this.copied = false;
        
        // setup matrixA
        int errorcode = shape(numRows, numCols);
        if(errorcode != 0) {
            System.err.println("Shape returned non-zero value: " + errorcode);
        }
    
    }
    

    /**
     * construct a <code>numRows</code> X <code>numCols</code> <code>SerialDenseMatrix</code>
     * filled with the elements of a 2-dimensional array<br>
     *
     * @param source 2-dimensional array of ints to be copied for the elements of <code>matrixA</code>
     * @param numRows number of rows in <code>this</code> <code>SerialDenseMatrix</code>
     * @param numCols number of columns in <code>this</code> <code>SerialDenseMatrix</code>
     */      
    public SerialDenseMatrix(double[][] source, int numRows, int numCols) {
        this.numRows = numRows;
        this.numCols = numCols;
        this.lda = numRows;
        this.copied = true;
        
        // setup matrixA
        matrixA = new double[this.lda * this.numCols];
        // do copy
        int i,j;
        int offSet;
        for (j=0; j < numCols; j++) {
            offSet = j*lda;
            
            for (i=0; i < numRows; i++) {
                matrixA[i + offSet] = source[i][j];
            }
        }
        
    
    }
    
    
    /**
     * construct a copy of the <code>source</code> <code>SerialDenseMatrix</code>
     *
     * @param source <code>SerialDenseMatrix</code> to be copied
     */
    public SerialDenseMatrix(SerialDenseMatrix source) {
        /*this.numRows = source.getNumRows();
        this.numCols = source.getNumCols();
        this.lda = source.getLda();
        this.copied = true;
        */
        
        // setup matrixA
        shape(source.getNumRows(), source.getNumCols());
        // do the copy
        copyMatrix(source.getMatrixA(), source.getLda(), source.getNumRows(), source.getNumCols(), matrixA, this.lda);
    }
    
    /**
     * Set dimensions of a <code>SerialDenseMatrix</code> object; initilalize values to zero.
     * <p>
     * Allows user to define the dimensions of a <code>SerialDenseMatrix</code> at any point. This method can
     * be called at any point after construction.  Any values that were previously in this object are
     * destroyed and the resized matrix starts off with all zero values.
     * 
     *  @param numRows number of rows in object
     *  @param numCols number of columns in object 
     */
    public int shape(int numRows, int numCols) {
        this.numRows = numRows;
        this.numCols = numCols;
        this.lda = numRows;
        this.copied = true;
        
        this.matrixA = new double[this.lda * this.numCols];
        // in C++ matrixA is then filled with 0s, Java does this by default
        
        return 0;
    }
    
    /**
     * Reshape a <code>SerialDenseMatrix</code> object.
     *
     * Allows user to define the dimensions of a <code>SerialDenseMatrix</code> at any point. This method can
     * be called at any point after construction.  Any values that were previously in this object are
     * copied into the new shape.  If the new shape is smaller than the original, the upper left portion
     * of the original matrix (the principal submatrix) is copied to the new matrix.
     *
     *  @param numRows number of rows in object
     *  @param numCols number of columns in object
     */
    public int reshape(int numRows, int numCols) {
        // Allocate space for new matrix
        double[] tmpMatrixA = new double[numRows*numCols];
        // values 0 by defualt
        int tmpRow = Math.min(this.numRows, numRows);
        int tmpCol = Math.min(this.numCols, numCols);
        if (matrixA != null) {
            copyMatrix(matrixA, lda, tmpRow, tmpCol, tmpMatrixA, numRows); // Copy principal submatrix of A to new A
        }
        
        /*DeleteArrays(); // Get rid of anything that might be already allocated  -- For Java add garbage collection explicitly?*/
        this.numRows = numRows;
        this.numCols = numCols;
        lda = this.numRows;
        matrixA = tmpMatrixA; // Set reference to new A
        copied = true;
        
        return 0;    
    }
    
    /**
     * Copy <code>myMatrixA</code> to <code>myMatrixB</code>.
     *
     * @param myMatrixA source matrix
     * @param myLdaA the LDA of myMatrixA
     * @param myNumRows the number of rows of myMatrixA
     * @param myNumcols the number of columns of myMatrixA
     * @param myMatrixB the destination matrix
     * @param myLdaB the LDA of myMatrixB
     */
    void copyMatrix(double[] myMatrixA, int myLdaA, int myNumRows, int myNumCols, double[] myMatrixB, int myLdaB) {
        int i, j;
        int offSetA, offSetB;
        for(j=0; j < myNumCols; j++) {
            offSetA = j*myLdaA;
            offSetB = j*myLdaB;
            
            for (i=0; i< myNumRows; i++) {
                myMatrixB[i + offSetB] = myMatrixA[i + offSetA];
            }
        }
    }
    
    /**
     * Computes the 1-Norm of the <code>this</code> matrix.
     *
     * @return 1-Norm of the <code>this</code> matrix
     */
    public double getOneNorm() { 
        double normA = 0;
        
        int i,j;
        double sum;
        int offSet;
        for (j=0; j < numCols; j++) {
            sum = 0;
            offSet = j*lda;
            
            for (i=0; i < numRows; i++) {
                sum += Math.abs(matrixA[i + offSet]);
            }
            
            normA = Math.max(normA, sum);
        }
        
        return normA;  
    }
    
    /**
     * Computes the Infinity-Norm of the <code>this</code> matrix.
     *
     * @return Infinity-Norm of the <code>this</code> matrix
     */
    public double getInfNorm() {
    
        double normA = 0;

        // Comment is from C++ ePetra:
        // Loop across columns in inner loop.  Most expensive memory access, but 
        // requires no extra storage.
        int i, j;
        double sum;
        int index;
        for (i=0; i < numRows; i++) {
            sum = 0;
            index = i;
            
            for (j=0; j < numCols; j++) {
                sum += Math.abs(matrixA[index]);
                index += lda;
            }
            
            normA = Math.max(normA, sum);
        }
        return normA;   
    }
    
    /**
     * accessor for the row dimension of system
     *
     * @return row dimension of system
     */
    public int getNumRows() {
        return this.numRows;
    }
    
    /**
     * accessor for the column dimension of system
     *
     * @return column dimension of system
     */
    public int getNumCols() {
        return this.numCols;
    }

    /**
     * accessor for the <code>matrixA</code>
     *
     * @return <code>matrixA</code>
     */  
    public double[] getMatrixA() {
        return matrixA;
    }
    
    /**
     * accessor for the leading dimension of the <code>this</code> matrix
     *
     * @return the leading dimension of the <code>this</code> matrix
     */
    public int getLda() {
        return this.lda;
    }
    
    /**
     * accessor for an element of the <code>this</code> matrix
     *
     * @param rowIndex row index for the element
     * @param colIndex column index for the element
     * @return element at <code>rowIndex</code>, <code>colIndex</code>
     */
    public double getElement(int rowIndex, int colIndex) {
    
        if (rowIndex >= numRows) {
            System.err.println("Row index = " + rowIndex + " Out of Range 0 - " + (numRows-1));
        }
        
        if (colIndex >= numCols) {
            System.err.println("Column index = " + colIndex + " Out of Range 0 - " + (numCols-1));
        }
        
        return matrixA[colIndex*lda + rowIndex];
    }
    
    /**
     * accessor for setting an element of the <code>this</code> matrix
     *
     * @param rowIndex row index for the element
     * @param colIndex column index for the element
     */
    public void setElement(int rowIndex, int colIndex, double value) {
        
        if (rowIndex >= numRows) {
            System.err.println("Row index = " + rowIndex + " Out of Range 0 - " + (numRows-1));
        }
        
        if (colIndex >= numCols) {
            System.err.println("Column index = " + colIndex + " Out of Range 0 - " + (numCols-1));
        }
        
        matrixA[colIndex*lda + rowIndex] = value;
    }
    
    public int  multiply (String transA, String transB, double scalarAB, SerialDenseMatrix myA, SerialDenseMatrix myB, double scalarThis ) {
        // Check for compatible dimensions
        if (transA.charAt(0) != 'T' && transA.charAt(0) != 'N') System.err.println("-2"); // Return erro
        if (transB.charAt(0) != 'T' && transB.charAt(0) != 'N') System.err.println("-3");
        
        int nrowsMyA = (transA.charAt(0) == 'T') ? myA.getNumCols() : myA.getNumRows();
        int ncolsmyA = (transA.charAt(0) == 'T') ? myA.getNumRows() : myA.getNumCols();
        int nrowsMyB = (transB.charAt(0) == 'T') ? myB.getNumCols() : myB.getNumRows();
        int ncolsMyB = (transB.charAt(0) == 'T') ? myB.getNumRows() : myB.getNumCols();
        
        if (this.numRows        != myA.getNumRows()  ||
            myA.getNumCols()   != myB.getNumRows()   ||
            this.numCols        != myB.getNumCols()    ) System.err.println("-1"); // Return error
        
        
        // Call BLAS GEMM function
        //GEMM(transA, transB, M_, N_, A_ncols, ScalarAB, A.A(), A.LDA(), 
        //   B.A(), B.LDA(), ScalarThis, A_, LDA_);
        
        // offSets are set to 0... don't think we need them
        Dgemm.dgemm(transA, transB, this.numRows, this.numCols, myA.getNumCols(), scalarAB, myA.getMatrixA(), myA.getLda(), 0,
           myB.getMatrixA(), myB.getLda(), 0, scalarThis, this.matrixA, this.lda, 0);
        
        /*
         * C++ from ePetra, need to ask about it
         * 
         * long int nflops = 2*M_;
         * nflops *= N_;
         * nflops *= A_ncols;
         * if (ScalarAB != 1.0) nflops += M_*N_;
         * if (ScalarThis != 0.0) nflops += M_*N_;
         * UpdateFlops(nflops);
         */
        
        return(0);
    }
    
    public int multiply (String sideA, double scalarAB, SerialSymDenseMatrix myA, SerialDenseMatrix myB, double scalarThis ) {
        // Check for compatible dimensions
        if (sideA.charAt(0) == 'R') {
            if (this.numRows != myB.getNumRows() || 
                this.numCols != myA.getNumCols() ||
                myB.getNumCols() != myA.getNumRows() ) {System.err.println("-1");} // Return error
        }
        else if (sideA.charAt(0) == 'L') {
            if (this.numRows != myA.getNumRows() || 
                numCols != myB.getNumCols() ||
                myA.getNumCols() != myB.getNumRows() ) {System.err.println("-1");}// Return error
        }
        else {
            System.err.println("-2"); // Return error, incorrect value for SideA
        }
        
        // Call BLAS SYMM function, offsets are set to 0
        Dsymm.dsymm(sideA, myA.getUpLo(), this.numRows, this.numCols, scalarAB, myA.getMatrixA(), 0, myA.getLda(), 
           myB.getMatrixA(), 0, myB.getLda(), scalarThis, this.matrixA, 0, this.lda);
        
        /*
         * C++ from ePetra, need to ask about it
         *
         * long int nflops = 2*M_;
         * nflops *= N_;
         * nflops *= A.N();
         * if (ScalarAB != 1.0) nflops += M_*N_;
         * if (ScalarThis != 0.0) nflops += M_*N_;
         * UpdateFlops(nflops);
         */
         
        return(0);
    }
    
    /**
     * Prints out the matrix.  For developer debugging purposes only.
     * May be removed or modified at any time without notice.
     */
    public void print() {
        int i, j;
        int index;
        for (i=0; i < numRows; i++) {
            index = i;
            
            System.out.print("[");
            for (j=0; j < numCols; j++) {
                System.out.print(" " + matrixA[index]);
                index += lda;
            }
            System.out.println(" ]");
        }
        
    }
}