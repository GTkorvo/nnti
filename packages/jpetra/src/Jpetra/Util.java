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
 * @author  Michael William Boldt
 * @author Jason Cross
 */
public class Util extends JpetraObject {
    public Util() {
    }
    
    public Util(Util util) {
    }
    
    public static int max(int x, int y) {
        return (( (x) > (y) ) ? x : y);
    }
    
    public static double max(double x, double y) {
        return (( (x) > (y) ) ? x : y);
    }
    
    public static double min(double x, double y) {
        return (( (x) < (y) ) ? x : y);
    }
    
    public static int min(int x, int y) {
        return (( (x) < (y) ) ? x : y);
    }
    
    public static void sort(boolean sortAscending, int[] keys, double[][] doubleCompanions, int[][] intCompanions) {
        int numKeys = keys.length;
        int numDoubleCompanions = doubleCompanions.length;
        int numIntCompanions = intCompanions.length;
        int i;
        int n = numKeys;
        final int [] list = keys;
        int m = n / 2;
        
        while(m > 0) {
            int max = n-m;
            for(int j=0; j<max; j++) {
                for(int k=j; k>=0; k-=m) {
                    if(list[k+m] >= list[k])
                        break;
                    int temp = list[k+m];
                    list[k+m] = list[k];
                    list[k] = temp;
                    
                    for(i=0; i<numDoubleCompanions; i++) {
                        double dTemp = doubleCompanions[i][k+m];
                        doubleCompanions[i][k+m] = doubleCompanions[i][k];
                        doubleCompanions[i][k] = dTemp;
                    }
                    for(i=0; i<numIntCompanions; i++) {
                        int iTemp = intCompanions[i][k+m];
                        intCompanions[i][k+m] = intCompanions[i][k];
                        intCompanions[i][k] = iTemp;
                    }
                }
            }
            m = m / 2;
        }
    }
}

