//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include "AztecOO.h"

#ifdef EPETRA_MPI
#define AZ_MPI
#define AZTEC_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_MsrMatrix.h"
#include "Epetra_LinearProblem.h"

int main(int argc, char *argv[])
{
  int    proc_config[AZ_PROC_SIZE];// Processor information.

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  AZ_set_proc_config(proc_config,MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
  AZ_set_proc_config(proc_config,0);
#endif

  cout << comm << endl;

  // int temp; if (comm.MyPID()==0) cin >> temp; comm.Barrier();

  if(argc != 2) cerr << "Error: enter name of data file on command line" << endl; 
  /* Set exact solution to NULL */

  int    *update;                  /* vector elements updated on this node. */
  int    *bindx;
  double *val;
  double *xguess, *b, *xexact, *xsolve;
  int    N_update;           /* # of block unknowns updated on this node    */
  int    numLocalEquations, numGlobalEquations, n_nonzeros;

  xexact = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

    Trilinos_Util_read_hb(argv[1], comm.MyPID(), &numGlobalEquations, &n_nonzeros,
             &val,  &bindx, &xguess, &b, &xexact);

  Trilinos_Util_distrib_msr_matrix(comm, &numGlobalEquations, &n_nonzeros, &N_update,
		  &update, &val, &bindx, &xguess, &b, &xexact);

#ifdef DEBUG
  for (i = 0; i<N_update; i++)
    if (val[i] == 0.0 ) cout << "Zero diagonal at row " << i << endl;
#endif

  /* convert matrix to a local distributed matrix */
  int * external = 0;
  int * update_index = 0;
  int * extern_index = 0;
  int * data_org = 0;

  AZ_transform(proc_config, &external, bindx, val, update,
	       &update_index, &extern_index, &data_org,
	       N_update, 0, 0, 0, 0,
               AZ_MSR_MATRIX);

  cout << comm << ": Completed AZ_transform." << endl;

  AZ_MATRIX * Amat = AZ_matrix_create(N_update);
  AZ_set_MSR(Amat, bindx, val, data_org, N_update, update, AZ_LOCAL);

  Epetra_MsrMatrix A(proc_config, Amat); // Create Epetra_MsrMatrix

  const Epetra_Map & map = A.RowMatrixRowMap();
  Epetra_Vector xx(Copy, map, xexact);

  Epetra_Vector bb(Copy, map, b);
  Epetra_Vector borig(Copy, map, b);


  // Create an initial guess of random numbers
  Epetra_Vector x(map); x.Random();
  
  // This next step is VERY important to get the right vector ordering
  // since AZ_transform potentially changes the matrix ordering
  AZ_reorder_vec(&(bb[0]), data_org, update_index, 0);
  AZ_reorder_vec(&(x[0]), data_org, update_index, 0);
  cout << "Building Epetra_LinearProblem" << endl;

  Epetra_LinearProblem problem(&A, &x, &bb);
  // Construct a solver object for this problem

  cout << "Building AztecOO solver" << endl;

  AztecOO solver(problem);


  // Assert symmetric
  // problem->AssertSymmetric();

  // Set Problem Difficulty Level
  //problem->SetPDL(easy);

  solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  //solver.SetAztecOption(AZ_precond, AZ_ls);
  //solver.SetAztecOption(AZ_scaling, 8);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
  //solver.SetAztecOption(AZ_output, 0);
  //solver.SetAztecOption(AZ_graph_fill, 2);
  //solver.SetAztecOption(AZ_overlap, 1);
  //solver.SetAztecOption(AZ_poly_ord, 9);
  //solver.SetAztecParam(AZ_ilut_fill, 1.0);
  //solver.SetAztecParam(AZ_drop, 0.0);
  //double rthresh = 1.4;
  //cout << "Rel threshold = " << rthresh << endl;
  // solver.SetAztecParam(AZ_rthresh, rthresh);
  //double athresh = 10.0;
  //cout << "Abs threshold = " << athresh << endl;
  //solver.SetAztecParam(AZ_athresh, athresh);
  //solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);


  

  //solver.SetAztecOption(AZ_reorder, 2);

  int Niters = 320;
  solver.SetAztecOption(AZ_kspace, 320);
   
  cout << "Computing norms" << endl;

  double normInf = A.NormInf();
  double normOne = A.NormOne();
  if (comm.MyPID()==0) 
    cout << "Inf-norm of A = " << normInf << endl
	 << "One-norm of A = " << normOne<< endl << endl;
  cout << "Solving problem" << endl;
  solver.Iterate(Niters, 1.0e-14);

  Epetra_Vector bcomp(map);
  assert(A.Multiply(false, x, bcomp)==0);
 
  Epetra_Vector resid(map);
 
  assert(resid.Update(1.0, bb, -1.0, bcomp, 0.0)==0);
  double residual;
  assert(resid.Norm2(&residual)==0);
  if (comm.MyPID()==0) cout << "Residual    = " << residual << endl;

  // Now we must apply inverse ordering to get the solution back
  Epetra_Vector xfinal(x);
  AZ_invorder_vec(&(x[0]), data_org, update_index, 0, &(xfinal[0]));
  assert(resid.Update(1.0, xx, -1.0, xfinal, 0.0)==0);

  assert(resid.Norm2(&residual)==0);
  if (comm.MyPID()==0)
    cout << "2-norm of difference between computed and exact solution  = " << residual << endl;

  if (residual>1.0e-5) {
    cout << "Difference between computed and exact solution is large..." << endl     
	 << "Computing norm of A times this difference.  "
	 << "If this norm is small, then matrix is singular"
	 << endl;
    assert(A.Multiply(false, resid, bcomp)==0);
    assert(bcomp.Norm2(&residual)==0);
    if (comm.MyPID()==0)
      cout << "2-norm of A times difference between computed and exact solution  = " 
	   << residual << endl;
    
  }

  free ((void *) xguess);
  free ((void *) b);
  free ((void *) xexact);
  free ((void *) val);
  free ((void *) bindx);
  free ((void *) update);

  // Must delete
  if (external!=0) AZ_free ((void *) external);
  if (update_index!=0) AZ_free ((void *) update_index);
  if (extern_index!=0) AZ_free ((void *) extern_index);
  if (data_org!=0) AZ_free ((void *) data_org);

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
