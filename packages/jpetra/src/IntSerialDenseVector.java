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

package Jpetra;

/**
 * The <code>IntSerialDenseVector</code> class enables the construction and use of integer-valued, 
 * dense vectors.  It derives from the IntSerialDenseMatrix class.
 * <p>
 * The <code>IntSerialDenseVector</code> class is intended to provide convenient vector notation but derives all signficant 
 * functionality from IntSerialDenseMatrix.
 * <p>
 * <b>Constructing <code>IntSerialDenseVector</code> Objects</b>
 * <p>
 * There are four <code>IntSerialDenseVector</code> constructors.  The first constructs a zero-length object which should be made
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
 * <b>Extracting Data from <code>IntSerialDenseVector</code> Objects</b>
 * <p>
 * Once a <code>IntSerialDenseVector</code> is constructed, it is possible to view the data via access functions.
 */
public class IntSerialDenseVector extends IntSerialDenseMatrix {
    /**
     * construct a <code>null</code> <code>IntSerialDenseVector</code><br>
     * needs to be sized with <code>zise()</code> and<br>
     * needs to be filled with <code>setElement()</code>
     */
    public IntSerialDenseVector() {
        super();
    }
    
    /**
     * construct a zero filled <code>legnth</code> long <code>IntSerialDenseVector</code><br>
     * needs to be filled with <code>setElement()</code>
     *
     * @param length the length of <code>this</code> <code>IntSerialDenseVector</code>
     */   
    public IntSerialDenseVector(int length) {    
        super(length, 1);
    }
    
    /**
     * construct a copy of the <code>source</code> <code>IntSerialDenseVector</code>
     *
     * @param source <code>IntSerialDenseVector</code> to be copied
     */
    public IntSerialDenseVector(IntSerialDenseVector source) {  
        super(source);
    }
    
    /**
     * construct a <code>length</code> long <code>IntSerialDenseVector</code> filled with the elements of an array<br>
     *
     * @param copied set to true if the array data should be copied
     * @param source array of ints to be copied for the elements of <code>this</code> vector
     */   
    public IntSerialDenseVector(boolean copied, int[] sourceData) {
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
     * Set the length of a <code>IntSerialDenseVector</code> object; initilalize values to zero.
     * <p>
     * Allows user to define the length of a <code>IntSerialDenseVector/code> at any point. This method can
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
     * Allows user to define the length of a <code>IntSerialDenseMatrix</code> at any point. This method can
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
    public int getElement(int index) {
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
    public void setElement(int index, int value) {
        if (index >= numRows) {
            System.err.println("Row index = " + index + " Out of Range 0 - " + (numRows-1));
        }
        
        matrixA[index] = value;
    }
    
}
