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

package test;

import Jpetra.*;
import java.util.Arrays;

public class CcjCommSumAllTest2 {
    boolean verbose = false;
    
    public static void main(String[] args) {
        new CcjCommSumAllTest2(args);
    }
    
    public CcjCommSumAllTest2(String[] args) {
        JpetraObject output = new JpetraObject();
        output.initializeOutput();
        Jpetra.CcjComm comm = new Jpetra.CcjComm("test/ccjhosts.txt");
        
        if(args.length > 0 && args[0].equals("-v")) verbose = true;
        
        if(verbose) System.out.println("Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
        
        int offSet = comm.getVnodeID()*4;
        int[] myIntArray = new int[]{offSet, offSet+1, offSet+2, offSet+3};
        double[] myDoubleArray = new double[]{offSet, offSet+1, offSet+2, offSet+3};
        
        int[] finalIntArray = new int[4];
        double[] finalDoubleArray = new double[4];
        int j=0;
        for (int i = 0; i < comm.getNumVnodes(); i++) {
            j=i*4;
            
            finalIntArray[0] += j;
            finalDoubleArray[0] += j;
            finalIntArray[1] += j+1;
            finalDoubleArray[1] += j+1;
            finalIntArray[2] += j+2;
            finalDoubleArray[2] += j+2;
            finalIntArray[3] += j+3;
            finalDoubleArray[3] += j+3;
        }
        
        int[] globalIntSums = comm.sumAll(myIntArray);
        double[] globalDoubleSums = comm.sumAll(myDoubleArray);
        
        if (!java.util.Arrays.equals(globalIntSums, finalIntArray)) {
            System.err.println("intSumAll Failed!");
        }
        if (!java.util.Arrays.equals(globalDoubleSums, finalDoubleArray)) {
            System.err.println("doubleSumAll Failed!");
        }
        
        if (verbose) {
            for(int i=0; i < globalIntSums.length; i++) {
                System.out.println(globalIntSums[i] + " " + globalDoubleSums[i]);
            }
            System.out.println("------");
            for(int i=0; i < finalIntArray.length; i++) {
                System.out.println(finalIntArray[i] + " " + finalDoubleArray[i]);
            }
        }
        
        System.exit(1);
    }
}
