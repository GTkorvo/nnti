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

// this class from MVIO was modified for Jpetra so that
// it would compile under Java 1.4 specifications instead of 1.5
// MVIO can be found at http://www.math.uib.no/~bjornoh/mtj/mvio/
// original MVIO copyright message:
/*
 * Copyright (C) 2003, 2004 Bjørn-Ove Heimsund
 * 
 * This file is part of MVIO.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 */

package Jpetra.MatrixMarketIO;

//import java.lang.StringBuffer;

/**
 * Contains information on a matrix in the
 * <a href="http://math.nist.gov/MatrixMarket">Matrix Market</a> exchange
 * format. Supports all valid matrices.
 */
public class MatrixInfo {

	/**
	 * What kind of numbers are stored
	 */
	public static final int REAL = 0;
        public static final int INTEGER= 1;
        public static final int COMPLEX = 2;
        public static final int PATTERN = 3;

	/**
	 * Symmetry structure of the matrix, if any
	 */
	public static final int GENERAL = 4;
        public static final int SYMMETRIC = 5;
        public static final int SKEWSYMMETRIC = 6;
        public static final int HERMITIAN = 7;

	/**
	 * True if the matrix is sparse, else false
	 */
	private boolean sparse;

	/**
	 * Type of data stored
	 */
	private int field;

	/**
	 * Matrix symmetry
	 */
	private int symmetry;

	/**
	 * Creates a specific type
	 * 
	 * @param sparse
	 *            True for sparse matrices, else false
	 * @param field
	 *            Type of data stored
	 * @param symmetry
	 *            Matrix symmetry
	 */
	public MatrixInfo(
		boolean sparse,
		int field,
		int symmetry) {
		this.sparse = sparse;
		this.field = field;
		this.symmetry = symmetry;

		validate();
	}

	/**
	 * Validates the representation
	 */
	private void validate() {
		if (isDense() && isPattern())
			throw new IllegalArgumentException("Matrix cannot be dense with pattern storage");
		if (isReal() && isHermitian())
			throw new IllegalArgumentException("Data cannot be real with hermitian symmetry");
		if (!isComplex() && isHermitian())
			throw new IllegalArgumentException("Data must be complex with hermitian symmetry");
		if (isPattern() && isSkewSymmetric())
			throw new IllegalArgumentException("Storage cannot be pattern and skew symmetrical");
	}

	/**
	 * Returns <code>true</code> if the matrix is in coordinate format, else
	 * <code>false</code>
	 */
	public boolean isSparse() {
		return sparse;
	}

	/**
	 * Returns <code>true</code> if the matrix is in coordinate format, else
	 * <code>false</code>
	 */
	public boolean isCoordinate() {
		return sparse;
	}

	/**
	 * Returns <code>true</code> if the matrix is in array format, else
	 * <code>false</code>
	 */
	public boolean isDense() {
		return !sparse;
	}

	/**
	 * Returns <code>true</code> if the matrix is in array format, else
	 * <code>false</code>
	 */
	public boolean isArray() {
		return !sparse;
	}

	/**
	 * Returns <code>true</code> if the matrix stores real numbers, else
	 * <code>false</code>
	 */
	public boolean isReal() {
		return field == REAL;
	}

	/**
	 * Returns <code>true</code> if the matrix stores integers, else <code>false</code>
	 */
	public boolean isInteger() {
		return field == INTEGER;
	}

	/**
	 * Returns <code>true</code> if the matrix stores complex numbers, else
	 * <code>false</code>
	 */
	public boolean isComplex() {
		return field == COMPLEX;
	}

	/**
	 * Returns <code>true</code> if the matrix does not store any numbers,
	 * else <code>false</code>
	 */
	public boolean isPattern() {
		return field == PATTERN;
	}

	/**
	 * Returns <code>true</code> if the matrix form is general, else <code>false</code>
	 */
	public boolean isGeneral() {
		return symmetry == GENERAL;
	}

	/**
	 * Returns <code>true</code> if the matrix is symmetrical, else <code>false</code>
	 */
	public boolean isSymmetric() {
		return symmetry == SYMMETRIC;
	}

	/**
	 * Returns <code>true</code> if the matrix is skew-symmetrical, else
	 * <code>false</code>
	 */
	public boolean isSkewSymmetric() {
		return symmetry == SKEWSYMMETRIC;
	}

	/**
	 * Returns <code>true</code> if the matrix is Hermitian, else <code>false</code>
	 */
	public boolean isHermitian() {
		return symmetry == HERMITIAN;
	}

	/**
	 * Returns a string representation of the specifier. Can be used to provide
	 * a header for writing to a file. It is a two-line output, which can look
	 * like this:
	 * 
	 * <pre>
	 *  %%MatrixMarket matrix coordinate real general
	 * </pre>
	 */
	public String toString() {
		//StringBuilder buf = new StringBuilder();
                StringBuffer buf = new StringBuffer();
		buf.append("%%MatrixMarket matrix ");

		if (isSparse())
			buf.append("coordinate ");
		else
			buf.append("array ");

		if (isReal())
			buf.append("real ");
		else if (isComplex())
			buf.append("complex ");
		else if (isPattern())
			buf.append("pattern ");
		else if (isInteger())
			buf.append("integer ");
		else // This should never happen
			throw new IllegalArgumentException("Unknown field specification");

		if (isGeneral())
			buf.append("general\n");
		else if (isSymmetric())
			buf.append("symmetric\n");
		else if (isSkewSymmetric())
			buf.append("skew-symmetric\n");
		else if (isHermitian())
			buf.append("Hermitian\n");
		else // This should never happen
			throw new IllegalArgumentException("Unknown symmetry specification");

		return buf.toString();
	}

}
