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
 *
 * @author  Jason Cross
 */
public class SerialDirectory extends JpetraObject implements Directory {
    private ElementSpace elementSpace;
    
    public SerialDirectory(ElementSpace elementSpace) {
        this.elementSpace = elementSpace;
    }
    
    public int[][] getDirectoryEntries(int[] globalElements) {
        int[][] vnodeIdsLocalElementSpaceIds = new int[globalElements.length][];
        
        for (int i=0; i < globalElements.length; i++) {
            vnodeIdsLocalElementSpaceIds[i] = new int[]{0, globalElements[i]};
        }
        return vnodeIdsLocalElementSpaceIds; 
    }
    
}
