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

import Jpetra.*;

public class IntSDMTest {
    public static void main(String[] args) {
        new IntSDMTest();
    }

    public IntSDMTest() {
        int numRows = 5;
        int numCols = 5;
        IntSerialDenseMatrix isdm = new IntSerialDenseMatrix(numRows, numCols);
        
        // fill up isdm
        int i, j;
        int k = 0;
        for(i=0; i < numRows; i++) {
            for(j=0; j < numCols; j++) {
                isdm.setElement(i, j, k++);
            }
        }
        
        isdm.print();
        
        IntSerialDenseMatrix isdm1 = new IntSerialDenseMatrix(isdm);
        
        isdm1.print();
        
        System.out.println("isdm OneNorm: " + isdm.getOneNorm());
        System.out.println("isdm1 OneNorm: " + isdm1.getOneNorm());
        
        System.out.println("isdm InfoNorm: " + isdm.getInfNorm());
        System.out.println("isdm1 InfoNorm: " + isdm1.getInfNorm());
        
        isdm1.reshape(numRows+2, numCols+2);
        
        isdm.print();
        isdm1.print();
        
        System.out.println("isdm OneNorm: " + isdm.getOneNorm());
        System.out.println("isdm1 OneNorm: " + isdm1.getOneNorm());
        
        System.out.println("isdm InfoNorm: " + isdm.getInfNorm());
        System.out.println("isdm1 InfoNorm: " + isdm1.getInfNorm());
        
        isdm1.reshape(numRows-2, numCols-2);
        isdm.shape(numRows-3, numCols-3);
        
        isdm.print();
        isdm1.print();
        
        System.out.println("isdm OneNorm: " + isdm.getOneNorm());
        System.out.println("isdm1 OneNorm: " + isdm1.getOneNorm());
        
        System.out.println("isdm InfoNorm: " + isdm.getInfNorm());
        System.out.println("isdm1 InfoNorm: " + isdm1.getInfNorm());
    }
}
