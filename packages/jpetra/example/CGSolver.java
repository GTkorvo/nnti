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

import Jpetra.*;

/**
 * @author Mike Heroux
 */

public class CGSolver {
    private Jpetra.CrsMatrix A;
    private Jpetra.Vector x;
    private Jpetra.Vector b;
    private Jpetra.Vector r;
    private Jpetra.Vector p;
    private Jpetra.Vector w;
    private int numIters;
    private double [] rho;
    private double alpha;
    private double beta;

    /**
     * Initialize member variables
     */
    public CGSolver(Jpetra.CrsMatrix A, Jpetra.Vector x, Jpetra.Vector b) 
	throws Jpetra.JpetraException {
	this.A = A;
	this.x = x;
	this.b = b;
	numIters = 0;
	rho = new double[3];
	p = new Jpetra.Vector(x.getMap());
	r = new Jpetra.Vector(x.getMap());
	w = new Jpetra.Vector(x.getMap());
    }

    /**
     * Perform maxIters number of iterations, or until 2-norm of residual < tolerance.  
     * Returns true if converged otherwise false.
     */
    public Jpetra.Vector iterate(int maxIters, double tolerance) {

	int numIters0 = 0; // number of iterations for this call to iterate
	double [] ptw = new double[1];

	if (numIters==0) {
	    A.multiply(x, p);
	    r.update(1.0, b, -1.0, p, 0.0); // r = b - A*x0
	    r.dot(r,rho); // rho[0] = ||r||^2
	    rho[1] = rho[0];
            System.out.println("Iteration = " + numIters +
                               "   Square of Norm of residual  " + rho[0]);
	}

	while (Math.sqrt(rho[0]) > tolerance && numIters0 < maxIters) {
	    rho[2] = rho[1];
	    rho[1] = rho[0];
	    numIters++;
	    numIters0++;
	    if (numIters==1) {
		p.update(1.0, r, 0.0); // p = r
	    }
	    else {
		beta = rho[1]/rho[2]; // beta_k = rho[k-1]/rho[k-2]
		p.update(1.0, r, beta); // p = r + beta_k*p
	    }
	    A.multiply(p,w); // w = A*p
	    w.dot(p,ptw); // ptw = p^T*w
	    alpha = rho[1]/ptw[0];
            //System.out.println("alpha = " + alpha);
	    //x.print("X before x update");
	    //p.print("P before x update");
	    x.update(alpha, p, 1.0); // x = x + alpha_k*p
	    //x.print("X after x update");
	    //r.print("R before r update");
	    //w.print("W before r update");
	    r.update(-alpha, w, 1.0); // r = r - alpha_k*w
	    //r.print("R after r update");
	    r.dot(r,rho); // rho[0] = ||r||^2	    
            //System.out.println("rho[2] = " + rho[2]);
            //System.out.println("rho[1] = " + rho[1]);
            //System.out.println("rho[0] = " + rho[0]);
            System.out.println("Iteration = " + numIters +
                               "   Square of Norm of residual  " + rho[0]);
	}
	return(x);
    }    

    /**
     * Extract reference to direction vector p.
     */
    public Jpetra.Vector getp() { return(p);}   

    /**
     * Extract reference to residual vector r.
     */
    public Jpetra.Vector getr() { return(r);}
}
