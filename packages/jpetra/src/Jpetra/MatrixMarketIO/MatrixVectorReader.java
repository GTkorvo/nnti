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

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.util.LinkedList;
import java.util.List;

import Jpetra.MatrixMarketIO.*;

/**
 * Reads matrices and vectors
 */
public class MatrixVectorReader extends BufferedReader {

	/**
	 * Reads the entries of the matrix or vector
	 */
	private StreamTokenizer st;

	/**
	 * Constructor for MatrixVectorReader
	 * 
	 * @param in
	 *            A Reader
	 */
	public MatrixVectorReader(Reader in) {
		super(in);
		setup();
	}

	/**
	 * Constructor for MatrixVectorReader
	 * 
	 * @param in
	 *            A Reader
	 * @param sz
	 *            Input buffer size
	 */
	public MatrixVectorReader(Reader in, int sz) {
		super(in, sz);
		setup();
	}

	/**
	 * Sets up the stream tokenizer
	 */
	private void setup() {
		st = new StreamTokenizer(this);
		st.resetSyntax();
		st.eolIsSignificant(false);
		st.lowerCaseMode(true);

		// Parse numbers as words
		st.wordChars('0', '9');
		st.wordChars('-', '.');

		// Characters as words
		st.wordChars('\u0000', '\u00FF');

		// Skip comments
		st.commentChar('%');

		// Skip whitespace and newlines
		st.whitespaceChars(' ', ' ');
		st.whitespaceChars('\u0009', '\u000e');
	}

	/**
	 * Shifts the indices. Useful for converting between 0- and 1-based
	 * indicing.
	 * 
	 * @param num
	 *            Added to every index
	 * @param indices
	 *            Indices to shift
	 */
	public void add(int num, int[] indices) {
		for (int i = 0; i < indices.length; ++i)
			indices[i] += num;
	}

	/**
	 * Reads the matrix info for the Matrix Market exchange format. The line
	 * must consist of exactly 5 space-separated entries, the first being
	 * "%%MatrixMarket"
	 */
	public MatrixInfo readMatrixInfo() throws IOException {
		String[] component = readLine().trim().split(" +");
		if (component.length != 5)
				throw new IOException(
						"Current line unparsable. It must consist of 5 tokens");

		// Read header
		if (!component[0].equalsIgnoreCase("%%MatrixMarket"))
				throw new IOException("Not in Matrix Market exchange format");

		// This will always be "matrix"
		if (!(component[1].equalsIgnoreCase("matrix") || component[1].equalsIgnoreCase("mtx")))
				throw new IOException("Expected \"matrix\", got "
						+ component[1]);

		// Sparse or dense?
		boolean sparse = false;
		if (component[2].equalsIgnoreCase("coordinate"))
			sparse = true;
		else if (component[2].equalsIgnoreCase("array"))
			sparse = false;
		else
			throw new IOException("Unknown layout " + component[2]);

		// Dataformat
		int field = -1;
		if (component[3].equalsIgnoreCase("real"))
			field = MatrixInfo.REAL;
		else if (component[3].equalsIgnoreCase("integer"))
			field = MatrixInfo.INTEGER;
		else if (component[3].equalsIgnoreCase("complex"))
			field = MatrixInfo.COMPLEX;
		else if (component[3].equalsIgnoreCase("pattern"))
			field = MatrixInfo.PATTERN;
		else
			throw new IOException("Unknown field specification " + component[3]);

		// Matrix pattern
		int symmetry = -1;
		if (component[4].equalsIgnoreCase("general"))
			symmetry = MatrixInfo.GENERAL;
		else if (component[4].equalsIgnoreCase("symmetric"))
			symmetry = MatrixInfo.SYMMETRIC;
		else if (component[4].equalsIgnoreCase("skew-symmetric"))
			symmetry = MatrixInfo.SKEWSYMMETRIC;
		else if (component[4].equalsIgnoreCase("Hermitian"))
			symmetry = MatrixInfo.HERMITIAN;
		else
			throw new IOException("Unknown symmetry specification "
					+ component[4]);

		// Pack together. This also verifies the format
		return new MatrixInfo(sparse, field, symmetry);
	}

	/**
	 * Reads the vector info for the Matrix Market exchange format. The line
	 * must consist of exactly 4 space-separated entries, the first being
	 * "%%MatrixMarket"
	 */
	public VectorInfo readVectorInfo() throws IOException {
		String[] component = readLine().trim().split(" +");
		if (component.length != 4)
				throw new IOException(
						"Current line unparsable. It must consist of 4 tokens");

		// Read header
		if (!component[0].equalsIgnoreCase("%%MatrixMarket"))
				throw new IOException("Not in Matrix Market exchange format");

		// This will always be "vector"
		if (!component[1].equalsIgnoreCase("vector"))
				throw new IOException("Expected \"vector\", got "
						+ component[1]);

		// Sparse or dense?
		boolean sparse = false;
		if (component[2].equalsIgnoreCase("coordinate"))
			sparse = true;
		else if (component[2].equalsIgnoreCase("array"))
			sparse = false;
		else
			throw new IOException("Unknown layout " + component[2]);

		// Dataformat
		int field = -1;
		if (component[3].equalsIgnoreCase("real"))
			field = VectorInfo.REAL;
		else if (component[3].equalsIgnoreCase("integer"))
			field = VectorInfo.INTEGER;
		else if (component[3].equalsIgnoreCase("complex"))
			field = VectorInfo.COMPLEX;
		else if (component[3].equalsIgnoreCase("pattern"))
			field = VectorInfo.PATTERN;
		else
			throw new IOException("Unknown field specification " + component[3]);

		// Pack together. This also verifies the format
		return new VectorInfo(sparse, field);
	}

	/**
	 * Checks if a Matrix Market header is present ("%%MatrixMarket")
	 * 
	 * @return True if a header was found, else false
	 * @throws IOException
	 */
	public boolean hasInfo() throws IOException {
		// Read a line, then skip back
		mark(1024);
		String[] component = readLine().trim().split(" +");
		reset();

		return component[0].equalsIgnoreCase("%%MatrixMarket");
	}

	/**
	 * Reads all the comments (lines starting with '%'). Positions the reader
	 * at the first non-comment line. Can only be called after reading the
	 * matrix or vector info. The comments read does not include '%' or the
	 * newline
	 */
	public String[] readComments() throws IOException {
		/*
                 * we don't need comments so this is disabled
                List list = new LinkedList();
		while (true) {
			mark(1024); // Line length equal 1024 at most
			String line = readLine().trim();
			if (line.length() > 0)
				if (line.charAt(0) != '%') {
					reset();
					break;
				} else
					list.add(line.substring(1));
		}

		//return list.toArray(new String[list.size()]);
                return (String[]) list.toArray();  //this doesn't cast properly.. danget
                 */
            return null;
	}

	/**
	 * Reads in the size of a matrix. Skips initial comments
	 */
	public MatrixSize readMatrixSize(MatrixInfo info) throws IOException {
		// Always read the matrix size
		int numRows = getInt(), numColumns = getInt();

		// For coordinate matrices we also read the number of entries
		if (info.isDense())
			return new MatrixSize(numRows, numColumns, info);
		else {
			int numEntries = getInt();
			return new MatrixSize(numRows, numColumns, numEntries);
		}
	}

	/**
	 * Reads in the size of an array matrix. Skips initial comments
	 */
	public MatrixSize readArraySize() throws IOException {
		int numRows = getInt(), numColumns = getInt();

		return new MatrixSize(numRows, numColumns, numRows * numColumns);
	}

	/**
	 * Reads in the size of a coordinate matrix. Skips initial comments
	 */
	public MatrixSize readCoordinateSize() throws IOException {
		int numRows = getInt(), numColumns = getInt(), numEntries = getInt();

		return new MatrixSize(numRows, numColumns, numEntries);
	}

	/**
	 * Reads in the size of a vector. Skips initial comments
	 */
	public VectorSize readVectorSize(VectorInfo info) throws IOException {
		// Always read the vector size
		int size = getInt();

		// For coordinate vectors we also read the number of entries
		if (info.isDense())
			return new VectorSize(size);
		else {
			int numEntries = getInt();
			return new VectorSize(size, numEntries);
		}
	}

	/**
	 * Reads in the size of a dense vector. Skips initial comments
	 */
	public VectorSize readVectorArraySize() throws IOException {
		int size = getInt();

		return new VectorSize(size);
	}

	/**
	 * Reads in the size of a coordinate vector. Skips initial comments
	 */
	public VectorSize readVectorCoordinateSize() throws IOException {
		int size = getInt(), numEntries = getInt();

		return new VectorSize(size, numEntries);
	}

	/**
	 * Reads the array data
	 */
	public void readArray(double[] data) throws IOException {
		int size = data.length;
		for (int i = 0; i < size; ++i)
			data[i] = getDouble();
	}

	/**
	 * Reads the array data
	 */
	public void readArray(float[] data) throws IOException {
		int size = data.length;
		for (int i = 0; i < size; ++i)
			data[i] = getFloat();
	}

	/**
	 * Reads the array data
	 */
	public void readArray(int[] data) throws IOException {
		int size = data.length;
		for (int i = 0; i < size; ++i)
			data[i] = getInt();
	}

	/**
	 * Reads the array data
	 */
	public void readArray(long[] data) throws IOException {
		int size = data.length;
		for (int i = 0; i < size; ++i)
			data[i] = getLong();
	}

	/**
	 * Reads the array data. The first array will contain real entries, while
	 * the second contain imaginary entries
	 */
	public void readArray(double[] dataR, double[] dataI) throws IOException {
		int size = dataR.length;
		if (size != dataI.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			dataR[i] = getDouble();
			dataI[i] = getDouble();
		}
	}

	/**
	 * Reads the array data. The first array will contain real entries, while
	 * the second contain imaginary entries
	 */
	public void readArray(float[] dataR, float[] dataI) throws IOException {
		int size = dataR.length;
		if (size != dataI.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			dataR[i] = getFloat();
			dataI[i] = getFloat();
		}
	}

	/**
	 * Reads a coordinate vector
	 */
	public void readCoordinate(int[] index, double[] data) throws IOException {
		int size = index.length;
		if (size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			index[i] = getInt();
			data[i] = getDouble();
		}
	}

	/**
	 * Reads a coordinate vector
	 */
	public void readCoordinate(int[] index, float[] data) throws IOException {
		int size = index.length;
		if (size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			index[i] = getInt();
			data[i] = getFloat();
		}
	}

	/**
	 * Reads a coordinate vector
	 */
	public void readCoordinate(int[] index, int[] data) throws IOException {
		int size = index.length;
		if (size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			index[i] = getInt();
			data[i] = getInt();
		}
	}

	/**
	 * Reads a coordinate vector
	 */
	public void readCoordinate(int[] index, long[] data) throws IOException {
		int size = index.length;
		if (size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			index[i] = getInt();
			data[i] = getLong();
		}
	}

	/**
	 * Reads a coordinate vector. First data array contains real entries, and
	 * the second contains imaginary entries
	 */
	public void readCoordinate(int[] index, float[] dataR, float[] dataI)
			throws IOException {
		int size = index.length;
		if (size != dataR.length || size != dataI.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			index[i] = getInt();
			dataR[i] = getFloat();
			dataI[i] = getFloat();
		}
	}

	/**
	 * Reads a coordinate vector. First data array contains real entries, and
	 * the second contains imaginary entries
	 */
	public void readCoordinate(int[] index, double[] dataR, double[] dataI)
			throws IOException {
		int size = index.length;
		if (size != dataR.length || size != dataI.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			index[i] = getInt();
			dataR[i] = getDouble();
			dataI[i] = getDouble();
		}
	}

	/**
	 * Reads a pattern vector
	 */
	public void readPattern(int[] index) throws IOException {
		int size = index.length;
		for (int i = 0; i < size; ++i)
			index[i] = getInt();
	}

	/**
	 * Reads a coordinate matrix
	 */
	public void readCoordinate(int[] row, int[] column, double[] data)
			throws IOException {
		int size = row.length;
		if (size != column.length || size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
			data[i] = getDouble();
		}
	}

	/**
	 * Reads a coordinate matrix
	 */
	public void readCoordinate(int[] row, int[] column, float[] data)
			throws IOException {
		int size = row.length;
		if (size != column.length || size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
			data[i] = getFloat();
		}
	}

	/**
	 * Reads a coordinate matrix
	 */
	public void readCoordinate(int[] row, int[] column, int[] data)
			throws IOException {
		int size = row.length;
		if (size != column.length || size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
			data[i] = getInt();
		}
	}

	/**
	 * Reads a coordinate matrix
	 */
	public void readCoordinate(int[] row, int[] column, long[] data)
			throws IOException {
		int size = row.length;
		if (size != column.length || size != data.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
			data[i] = getLong();
		}
	}

	/**
	 * Reads a pattern matrix
	 */
	public void readPattern(int[] row, int[] column) throws IOException {
		int size = row.length;
		if (size != column.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
		}
	}

	/**
	 * Reads a coordinate matrix. First data array contains real entries, and
	 * the second contains imaginary entries
	 */
	public void readCoordinate(int[] row, int[] column, double[] dataR,
			double[] dataI) throws IOException {
		int size = row.length;
		if (size != column.length || size != dataR.length || size != dataI.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
			dataR[i] = getDouble();
			dataI[i] = getDouble();
		}
	}

	/**
	 * Reads a coordinate matrix. First data array contains real entries, and
	 * the second contains imaginary entries
	 */
	public void readCoordinate(int[] row, int[] column, float[] dataR,
			float[] dataI) throws IOException {
		int size = row.length;
		if (size != column.length || size != dataR.length || size != dataI.length)
				throw new IllegalArgumentException(
						"All arrays must be of the same size");
		for (int i = 0; i < size; ++i) {
			row[i] = getInt();
			column[i] = getInt();
			dataR[i] = getFloat();
			dataI[i] = getFloat();
		}
	}

	/**
	 * Reads an integer
	 */
	public int getInt() throws IOException {
		st.nextToken();
		if (st.ttype == StreamTokenizer.TT_WORD)
			return Integer.parseInt(st.sval);
		else if (st.ttype == StreamTokenizer.TT_EOF)
			throw new EOFException("End-of-File encountered during parsing");
		else
			throw new IOException("Unknown token found during parsing");
	}

	/**
	 * Reads a long
	 */
	public long getLong() throws IOException {
		st.nextToken();
		if (st.ttype == StreamTokenizer.TT_WORD)
			return Long.parseLong(st.sval);
		else if (st.ttype == StreamTokenizer.TT_EOF)
			throw new EOFException("End-of-File encountered during parsing");
		else
			throw new IOException("Unknown token found during parsing");
	}

	/**
	 * Reads a double
	 */
	public double getDouble() throws IOException {
		st.nextToken();
		if (st.ttype == StreamTokenizer.TT_WORD)
			return Double.parseDouble(st.sval);
		else if (st.ttype == StreamTokenizer.TT_EOF)
			throw new EOFException("End-of-File encountered during parsing");
		else
			throw new IOException("Unknown token found during parsing");
	}

	/**
	 * Reads a float
	 */
	public float getFloat() throws IOException {
		st.nextToken();
		if (st.ttype == StreamTokenizer.TT_WORD)
			return Float.parseFloat(st.sval);
		else if (st.ttype == StreamTokenizer.TT_EOF)
			throw new EOFException("End-of-File encountered during parsing");
		else
			throw new IOException("Unknown token found during parsing");
	}

}
