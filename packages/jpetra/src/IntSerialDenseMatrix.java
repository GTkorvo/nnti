package Jpetra;

// Math provides abs(), max(), and min()
import java.lang.Math;
 
/**
 * The <code>IntSerialDenseMatrix</code> class is intended to provide very basic support for dense rectangular matrices.
 * <p>
 * <b>Constructing <code>IntSerialDenseMatrix</code> Objects</b>
 * <p>
 * There are four <code>IntSerialDenseMatrix</code> constructors.  The first constructs a zero-sized object which should be made
 * to appropriate length using the <code>shape()</code> or <code>reshape()</code> methods and then filled with <code>setElement()</code>. 
 * The second constructs an object sized to the dimensions specified, which should be filled with <code>setElement()</code>.
 * The third is a constructor that accepts user data as a 2D array (<code>int[][]</code>), and the fourth is a copy constructor.
 * <p>
 * <b>Extracting Data from <code>IntSerialDenseMatrix</code> Objects</b>
 * <p>
 * Once an <code>IntSerialDenseMatrix</code> is constructed, it is possible to view the data via access functions.
 * <p>
 * <b>Vector and Utility Functions</b>
 * <p>
 * Once an IntSerialDenseMatrix is constructed, several mathematical functions can be applied to
 * the object.  Specifically:
 * <ul>
 *   <li> Multiplication.
 *   <li> Norms.
 * </ul>
 *
 * The IntSerialDenseMatrix code originated as C++ code in ePetra and was ported to Java.
 * @author Jason Cross
 */ 
 public class IntSerialDenseMatrix extends JpetraObject {
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
    int[] matrixA;
    
    
    /**
     * construct a <code>null</code> <code>IntSerialDenseMatrix</code><br>
     * needs to be shaped with <code>shape()</code> and<br>
     * needs to be filled with <code>setElement()</code>
     */
    public IntSerialDenseMatrix() {
        this.numRows = 0;
        this.numCols = 0;
        this.lda = 0;
        this.copied = false;
        // initilize matrixA?
    }
    
    /**
     * construct a zero filled <code>numRows</code> X <code>numCols</code> <code>IntSerialDenseMatrix</code><br>
     * needs to be filled with <code>setElement()</code>
     *
     * @param numRows number of rows in <code>this</code> <code>IntSerialDenseMatrix</code>
     * @param numCols number of columns in <code>this</code> <code>IntSerialDenseMatrix</code>
     */    
    public IntSerialDenseMatrix(int numRows, int numCols) {
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
     * construct a <code>numRows</code> X <code>numCols</code> <code>IntSerialDenseMatrix</code> filled with the elements of a 2-dimensional array<br>
     *
     * @param source 2-dimensional array of ints to be copied for the elements of <code>matrixA</code>
     * @param numRows number of rows in <code>this</code> <code>IntSerialDenseMatrix</code>
     * @param numCols number of columns in <code>this</code> <code>IntSerialDenseMatrix</code>
     */      
    public IntSerialDenseMatrix(int[][] source, int numRows, int numCols) {
        this.numRows = numRows;
        this.numCols = numCols;
        this.lda = numRows;
        this.copied = true;
        
        // setup matrixA
        matrixA = new int[this.lda * this.numCols];
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
     * construct a copy of the <code>source</code> <code>IntSerialDenseMatrix</code>
     *
     * @param source <code>IntSerialDenseMatrix</code> to be copied
     */
    public IntSerialDenseMatrix(IntSerialDenseMatrix source) {
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
     * Set dimensions of a <code>IntSerialDenseMatrix</code> object; initilalize values to zero.
     * <p>
     * Allows user to define the dimensions of a <code>IntSerialDenseMatrix</code> at any point. This method can
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
        
        this.matrixA = new int[this.lda * this.numCols];
        // in C++ matrixA is then filled with 0s, Java does this by default
        
        return 0;
    }
    
    /**
     * Reshape a <code>IntSerialDenseMatrix</code> object.
     *
     * Allows user to define the dimensions of a <code>IntSerialDenseMatrix</code> at any point. This method can
     * be called at any point after construction.  Any values that were previously in this object are
     * copied into the new shape.  If the new shape is smaller than the original, the upper left portion
     * of the original matrix (the principal submatrix) is copied to the new matrix.
     *
     *  @param numRows number of rows in object
     *  @param numCols number of columns in object
     */
    public int reshape(int numRows, int numCols) {
        // Allocate space for new matrix
        int[] tmpMatrixA = new int[numRows*numCols];
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
    void copyMatrix(int[] myMatrixA, int myLdaA, int myNumRows, int myNumCols, int[] myMatrixB, int myLdaB) {
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
    public int getOneNorm() { 
        int normA = 0;
        
        int i,j;
        int sum, offSet;
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
    public int getInfNorm() {
    
        int normA = 0;

        // Comment is from C++ ePetra:
        // Loop across columns in inner loop.  Most expensive memory access, but 
        // requires no extra storage.
        int i, j;
        int sum, index;
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
    public int[] getMatrixA() {
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
    public int getElement(int rowIndex, int colIndex) {
    
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
    public void setElement(int rowIndex, int colIndex, int value) {
        
        if (rowIndex >= numRows) {
            System.err.println("Row index = " + rowIndex + " Out of Range 0 - " + (numRows-1));
        }
        
        if (colIndex >= numCols) {
            System.err.println("Column index = " + colIndex + " Out of Range 0 - " + (numCols-1));
        }
        
        matrixA[colIndex*lda + rowIndex] = value;
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