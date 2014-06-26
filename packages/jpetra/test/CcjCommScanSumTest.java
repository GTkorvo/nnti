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

public class CcjCommScanSumTest {
    boolean verbose = false;
    
    public static void main(String[] args) {
        new CcjCommScanSumTest(args);
    }
    
    public CcjCommScanSumTest(String[] args) {
        JpetraObject output = new JpetraObject();
        output.initializeOutput();
        
        Jpetra.CcjComm comm = new Jpetra.CcjComm("test/ccjhosts.txt");
        if (comm.getNumVnodes() !=4 ) {
            System.out.println("Fatal Error: This test must be run with 4 virutal nodes!");
        }
        if(args.length > 0 && args[0].equals("-v")) verbose = true;
        
        if(verbose) System.out.println("Processor "+comm.getVnodeID()+" of "+comm.getNumVnodes());
        
        int[] myIntArray;
        double[] myDoubleArray;
        if (comm.getVnodeID()==0) {
            myIntArray = new int[]{5,7,2,1000};
            myDoubleArray = new double[]{5,7,2,1000};
        }
        else if (comm.getVnodeID()==1) {
            myIntArray = new int[]{2,1,98,3};
            myDoubleArray = new double[]{2,1,98,3};
        }
        else if (comm.getVnodeID()==2) {
            myIntArray = new int[]{-2,4,102,9};
            myDoubleArray = new double[]{-2,4,102,9};
        }
        else {
            myIntArray = new int[]{5,6,7,8};
            myDoubleArray = new double[]{5,6,7,8};
        }
        
        int[] myIntScanSum = comm.scanSums(myIntArray);
        double[] myDoubleScanSum = comm.scanSums(myDoubleArray);
        
        // need to finish fixing this
        if (comm.getVnodeID()==0) {
            if (myIntScanSum[0] != 5 || myIntScanSum[1] != 7 || myIntScanSum[2] != 2 || myIntScanSum[3] != 1000) {
                System.out.println("intScanSum failed for vnode0");
            }
            if (!java.util.Arrays.equals(myDoubleScanSum, new double[]{5,7,2,1000})) {
                System.out.println("doubleScanSum failed for vnode0");
            }
        }
        else if (comm.getVnodeID()==1) {
            if (!java.util.Arrays.equals(myIntScanSum, new int[]{7,8,100,1003})) {
                System.out.println("intScanSum failed for vnode1");
            }
            if (!java.util.Arrays.equals(myDoubleScanSum, new double[]{7,8,100,1003})) {
                System.out.println("doubleScanSum failed for vnode1");
            }
        }
        else if (comm.getVnodeID()==2) {
            if (!java.util.Arrays.equals(myIntScanSum, new int[]{5,12,202,1012})) {
                System.out.println("intScanSum failed for vnode2");
            }
            if (!java.util.Arrays.equals(myDoubleScanSum, new double[]{5,12,202,1012})) {
                System.out.println("doubleScanSum failed for vnode2");
            }
        }
        else {
            if (!java.util.Arrays.equals(myIntScanSum, new int[]{10,18,209,1020})) {
                System.out.println("intScanSum failed for vnode3");
            }
            if (!java.util.Arrays.equals(myDoubleScanSum, new double[]{10,18,209,1020})) {
                System.out.println("doubleScanSum failed for vnode3");
            }
        }
        
        if (verbose) {
            for(int i=0; i < myIntScanSum.length; i++) {
                System.out.println(myIntScanSum[i] + " " + myDoubleScanSum[i]);
            }
        }
        
        System.exit(1);
    }
}
