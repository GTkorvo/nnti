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

/* based on IntArray.java which is Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package Jpetra.CcjSupport;
import CCJ.*;

import java.io.Serializable;


class CcjScatter2DArray implements Partitionable {
    public int [][] array;
    
    public CcjScatter2DArray(int[][] in) {
        array = in;
    }
    
    public int size() {
        return array.length;
    }
    
    public void setElementAt(int index, int groupSize, Serializable object) {
        array[index]=(int []) object;
    }
    
    public Serializable elementAt(int index, int groupSize) {
        return array[index];
    }
    
    public String toString() {
        return "IntArray([" + array.length + "][])";
    }
}
