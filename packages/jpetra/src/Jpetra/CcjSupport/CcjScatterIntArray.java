/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package Jpetra.CcjSupport;

import CCJ.*;
import java.io.Serializable;

class CcjScatterIntArray implements Partitionable {

    public int [] array;
    public int n;

    public CcjScatterIntArray(int[] array) {
	//array = new int[length];
	//n = length;
        this.array = array;
        n = this.array.length;
    }


    public int size() {
	return array.length;
    }


    public void setElementAt(int index, int groupSize, Serializable object) {

	int [] other = (int []) object;

	int d = n / groupSize;
	int m = n % groupSize;

	int offset = index * d;
	int size = d;

	if (index >= (groupSize - m)) {
	    size += 1;
	    offset += index - (groupSize - m);
	}

	System.arraycopy(other, 0, array, offset, size);
    }


    public Serializable elementAt(int index, int groupSize) {
	int size = array.length / groupSize;
	if ((array.length % groupSize) != 0) {
	    throw new RuntimeException("Unbalanced distribution not supported.");
	}
	int [] other = new int[size];
	int myIndex = index * size;

	System.arraycopy(array, myIndex, other, 0, size);
	return other;
    }


    public String toString() {
	return "IntArray([" + array.length + "])";
    }

}


