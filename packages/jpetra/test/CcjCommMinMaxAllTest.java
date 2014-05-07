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

public class CcjCommMinMaxAllTest {
    boolean verbose = false;
    
    public static void main (String[] args) {
        new CcjCommMinMaxAllTest(args);
     }
     
     public CcjCommMinMaxAllTest(String[] args) {
         JpetraObject output = new JpetraObject();
         output.initializeOutput();
         
         Jpetra.CcjComm comm = new Jpetra.CcjComm("test/ccjhosts.txt");
         if (comm.getNumVnodes() !=3 ) {
             System.out.println("Fatal Error: This test must be run with 3 virutal nodes!");
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
         else {
             myIntArray = new int[]{-2,4,102,9};
             myDoubleArray = new double[]{-2,4,102,9};
         }
         
         int[] realIntMax;
         int[] realIntMin;
         double[] realDoubleMax;
         double[] realDoubleMin;
         
         if (comm.getNumVnodes() == 1) {
             realIntMax = new int[]{5,7,2,1000};
             realDoubleMax = new double[]{5,7,2,1000};
             
             realIntMin = new int[]{5,7,2,1000};
             realDoubleMin = new double[]{5,7,2,1000};
         }
         else if (comm.getNumVnodes() == 2) {
             realIntMax = new int[]{5,7,98,1000};
             realDoubleMax = new double[]{5,7,98,1000};
             
             realIntMin = new int[]{2,1,2,3};
             realDoubleMin = new double[]{2,1,2,3};                  
         }
         else { 
             realIntMax = new int[]{5,7,102,1000};
             realDoubleMax = new double[]{5,7,102,1000};
             
             realIntMin = new int[]{-2,1,2,3};
             realDoubleMin = new double[]{-2,1,2,3};               
         }
         //int[] myIntMax = new int[myIntArray.length];
         //System.arraycopy(myIntArray, 0, myIntMax, 0, myIntArray.length);
         int[] myIntMax = comm.maxAll(myIntArray);
         
         //double[] myDoubleMax = new double[myDoubleArray.length];
         //System.arraycopy(myDoubleArray, 0, myDoubleMax, 0, myDoubleArray.length);
         double[] myDoubleMax = comm.maxAll(myDoubleArray);
         
         //int[] myIntMin = new int[myIntArray.length];
         //System.arraycopy(myIntArray, 0, myIntMin, 0, myIntArray.length);
         int[] myIntMin = comm.minAll(myIntArray);
         
         //double[] myDoubleMin = new double[myDoubleArray.length];
         //System.arraycopy(myDoubleArray, 0, myDoubleMin, 0, myDoubleArray.length);
         double[] myDoubleMin = comm.minAll(myDoubleArray);        
         
         if (!java.util.Arrays.equals(myIntMax, realIntMax)) {
             System.err.println("intMax Failed!");
         }
         if (!java.util.Arrays.equals(myIntMin, realIntMin)) {
             System.err.println("intMin Failed!");
         }
         if (!java.util.Arrays.equals(myDoubleMax, realDoubleMax)) {
             System.err.println("doubleMax Failed!");
         }
         if (!java.util.Arrays.equals(myDoubleMin, realDoubleMin)) {
             System.err.println("doubleMin Failed!");
         }
         if (verbose) {
             for(int i=0; i < myIntMax.length; i++) {
                 System.out.println(myIntMax[i] + " " + myDoubleMax[i] + " " + myIntMin[i] + " " + myDoubleMin[i]);   
             }
         }        
         
         System.exit(0);        
     }  
}
