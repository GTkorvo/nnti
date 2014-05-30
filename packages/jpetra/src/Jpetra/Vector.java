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

/**
 * <code>Vector</code> is a simple wrapper class to <code>MultiVector</code>
 * which creates a <code>MultiVector</code> with only one column.
 *
 * @author  Jason Cross
 */
public class Vector extends MultiVector {
    public Vector(VectorSpace vectorSpace) {
        super(vectorSpace);
    }
    
    public Vector(VectorSpace vectorSpace, double[] vectorValues) {
        super(vectorSpace);
        this.values = new double[1][];
        this.values[0] = vectorValues;
    }
}
