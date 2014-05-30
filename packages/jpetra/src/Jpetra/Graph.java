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

package Jpetra;

import java.util.TreeMap;

/**
 *
 * @author  Jason Cross
 */
public class Graph {
    private int[] nonZeroEntries;  // contigous 2-d array of non-zero ID entries
    private int[] numEntries;  // number of entries per array
    VectorSpace dVectorSpace; // the domain vector space
    /*private boolean filled;*/
    //boolean rowOriented;
    /*private TreeMap OuterMap;*/
    
    public Graph(VectorSpace vectorSpace, boolean rowOriented) {
        //this.rowOriented = rowOriented;
        //this.filled = false;
        this.dVectorSpace = vectorSpace;
    }
    
    public Graph(VectorSpace vectorSpace, int[] nonZeroEntries, int[] numEntries) {
        //this.rowOriented = rowOriented;
        //this.filled = true;
        this.dVectorSpace = vectorSpace;
        this.nonZeroEntries = nonZeroEntries;
        this.numEntries = numEntries;
        
    }
    
    //public boolean isRowOriented() {
    //    return this.rowOriented;
    //}
    
    
    
    public int getNumNonZeros() {
        return this.nonZeroEntries.length;
    }
    
    public int[] getNumEntriesArray() {
        return this.numEntries;
    }
    
    public int getIndex(int entryId) {
        return nonZeroEntries[entryId];
    }
    
    public int[] getNonZeroEntriesArray() {
        return this.nonZeroEntries;
    }
}
