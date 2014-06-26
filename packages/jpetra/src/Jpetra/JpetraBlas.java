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
public class JpetraBlas implements Blas {
    public JpetraBlas() {
        // empty
    }
    
    public double asum(double[] x) {
        double result = 0;
        for(int i=0; i < x.length; i++) {
            result += Math.abs(x[i]);
        }
        
        return result;
    }
    
    public double dot(double[] x, double[] y) {
        double result = 0;
        for(int i=0; i < x.length; i++) {
            result += x[i] * y[i];
        }
        
        return result;
    }
    
    public int iamax(double[] x) {
        double max = Math.abs(x[0]);
        int maxIndex = 0;
        for(int i=1; i < x.length; i++) {
            if (max < Math.abs(x[i])) {
                max = Math.abs(x[i]);
                maxIndex = i;
            }
        }
        
        return maxIndex;
    }
    
    public double norm2(double[] x) {
        double result = 0;
        
        for(int i=0; i < x.length; i++) {
            result += x[i] * x[i];
        }
        
        return Math.sqrt(result);
    }
    
    public void scale(double scalar, double[] x) {
        for(int i=0; i < x.length; i++) {
            x[i] = x[i] * scalar;
        }
    }
    
}
