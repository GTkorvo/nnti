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

public class IntSDVTest {
    public static void main(String[] args) {
        new IntSDVTest();
    }

    public IntSDVTest() {
        int length = 5;
        IntSerialDenseVector isdv = new IntSerialDenseVector(length);
        
        // fill up isdv
        int i;
        int k = 0;
        for(i=0; i < length; i++) {
            isdv.setElement(i, k++);
        }
        
        System.out.println("isdv:");
        isdv.print();
        
        IntSerialDenseVector isdv1 = new IntSerialDenseVector(isdv);

        System.out.println("isdv1:");        
        isdv1.print();
        
        System.out.println("isdm OneNorm: " + isdv.getOneNorm());
        System.out.println("isdv1 OneNorm: " + isdv1.getOneNorm());
        
        System.out.println("isdv InfoNorm: " + isdv.getInfNorm());
        System.out.println("isdv1 InfoNorm: " + isdv1.getInfNorm());
        
        isdv1.resize(length+2);

        System.out.println("isdv:");        
        isdv.print();
        System.out.println("isdv1:");        
        isdv1.print();
        
        System.out.println("isdv OneNorm: " + isdv.getOneNorm());
        System.out.println("isdv1 OneNorm: " + isdv1.getOneNorm());
        
        System.out.println("isdv InfoNorm: " + isdv.getInfNorm());
        System.out.println("isdv1 InfoNorm: " + isdv1.getInfNorm());
        
        isdv1.resize(length-2);
        isdv.size(length-3);
        
        System.out.println("isdv:");        
        isdv.print();
        System.out.println("isdv1:");        
        isdv1.print();
        
        System.out.println("isdv OneNorm: " + isdv.getOneNorm());
        System.out.println("isdv1 OneNorm: " + isdv1.getOneNorm());
        
        System.out.println("isdv InfoNorm: " + isdv.getInfNorm());
        System.out.println("isdv1 InfoNorm: " + isdv1.getInfNorm());
    }
}
