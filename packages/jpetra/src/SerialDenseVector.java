package Jpetra;

/**
 * The <code>SerialDenseVector</code> class enables the construction and use of integer-valued, 
 * dense vectors.  It derives from the SerialDenseMatrix class.
 * <p>
 * The <code>SerialDenseVector</code> class is intended to provide convenient vector notation but derives all signficant 
 * functionality from SerialDenseMatrix.
 * <p>
 * <b>Constructing <code>SerialDenseVector</code> Objects</b>
 * <p>
 * There are four <code>SerialDenseVector</code> constructors.  The first constructs a zero-length object which should be made
 * to appropriate length using the <code>size()</code> or <code>resize()</code> methods and and then filled with <code>setElement()</code>. 
 * The second constructs an object sized to the dimension specified, which should be filled with <code>setElement()</code>. 
 * The third is a constructor that accepts user data as a 1D array, and the fourth is a copy constructor.
 * <p>
 * The third constructor has two data access modes (specified by the DataAccess argument):
 * <ol>
 *   <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
 *        user data is not needed after construction.
 *   <li> View mode - Creates a "view" of the user data. In this case, the
 *        user data is required to remain intact for the life of the object.
 * </ol>
 * 
 * <b>Warning:</b> View mode is extremely dangerous from a data hiding perspective.
 * Therefore, we strongly encourage users to develop code using Copy mode first and 
 * only use the View mode in a secondary optimization phase.
 * <p>
 * <b>Extracting Data from <code>SerialDenseVector</code> Objects</b>
 * <p>
 * Once a <code>SerialDenseVector</code> is constructed, it is possible to view the data via access functions.
 */
public class SerialDenseVector extends SerialDenseMatrix {
    /**
     * construct a <code>null</code> <code>SerialDenseVector</code><br>
     * needs to be sized with <code>zise()</code> and<br>
     * needs to be filled with <code>setElement()</code>
     */
    public SerialDenseVector() {
        super();
    }
    
    /**
     * construct a zero filled <code>legnth</code> long <code>SerialDenseVector</code><br>
     * needs to be filled with <code>setElement()</code>
     *
     * @param length the length of <code>this</code> <code>SerialDenseVector</code>
     */   
    public SerialDenseVector(int length) {    
        super(length, 1);
    }
    
    /**
     * construct a copy of the <code>source</code> <code>SerialDenseVector</code>
     *
     * @param source <code>SerialDenseVector</code> to be copied
     */
    public SerialDenseVector(SerialDenseVector source) {  
        super(source);
    }
    
    /**
     * construct a <code>length</code> long <code>SerialDenseVector</code> filled with the elements of an array<br>
     *
     * @param copied set to true if the array data should be copied
     * @param source array of ints to be copied for the elements of <code>this</code> vector
     */   
    public SerialDenseVector(boolean copied, double[] sourceData) {
        this.numRows = sourceData.length;
        this.numCols = 1;
        this.lda = numRows;
        this.copied = copied;
        
        if (copied == true) {
            System.arraycopy(sourceData, 0, matrixA, 0, this.numRows);
        }
        else {
            matrixA = sourceData;
        }    
    }
    
    /**
     * Set the length of a <code>SerialDenseVector</code> object; initilalize values to zero.
     * <p>
     * Allows user to define the length of a <code>SerialDenseVector/code> at any point. This method can
     * be called at any point after construction.  Any values that were previously in this object are
     * destroyed and the resized vector starts off with all zero values.
     * 
     *  @param length length of the vector
     */    
    public void size(int length) {
        shape(length, 1);
    }

    /**
     * Reshape a <code>IntSerialDense</code> object.
     *
     * Allows user to define the length of a <code>SerialDenseMatrix</code> at any point. This method can
     * be called at any point after construction.  Any values that were previously in this object are
     * copied into the new length.  If the new length is smaller than the original, the top of the original vector
     * is copied into the new vector.
     *
     *  @param numRows number of rows in object
     *  @param numCols number of columns in object
     */    
    public void resize(int length) {
        reshape(length, 1);
    }
    
    /**
     * accessor for the length of the <code>this</code> vector
     *
     * @return the length of the <code>this</code> vector
     */    
    public int getLength() {
        return numRows;
    }    

    /**
     * accessor for an element of the <code>this</code> vector
     *
     * @param index index for the element
     * @return element at <code>position</code>
     */    
    public double getElement(int index) {
        if (index >= numRows) {
            System.err.println("Row index = " + index + " Out of Range 0 - " + (numRows-1));
        }
        
        return matrixA[index];
    }

    /**
     * accessor for setting an element of the <code>this</code> vector
     *
     * @param position index for the element to set
     * @param value value to set the element to
     */      
    public void setElement(int index, double value) {
        if (index >= numRows) {
            System.err.println("Row index = " + index + " Out of Range 0 - " + (numRows-1));
        }
        
        matrixA[index] = value;
    }
    
}