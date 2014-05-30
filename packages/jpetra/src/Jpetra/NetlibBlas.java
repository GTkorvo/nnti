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

import org.netlib.blas.Ddot;
import org.netlib.blas.Dnrm2;
import org.netlib.blas.Dscal;
import org.netlib.blas.Dasum;
import org.netlib.blas.Idamax;
/**
 *
 * @author  Jason Cross
 */
public class NetlibBlas extends JpetraObject implements Blas {
    
    public void NetlibBlas () {}
    
    public double dot(double[] x, double[] y) {
        return Ddot.ddot(x.length-1, x, 1, 1, y, 1, 1);
    }
    
    public double norm2(double[] x) {
        return Dnrm2.dnrm2(x.length-1,x,1,1);
    }
    
    public void scale(double scalar, double[] x) {
        Dscal.dscal(x.length,scalar,x,0,1);
    }
    
    public double asum(double[] x) {
        return Dasum.dasum(x.length-1,x,1,1);
    }
    
    public int iamax(double[] x) {
        return Idamax.idamax(x.length-1,x,1,1);
    }
    
}
