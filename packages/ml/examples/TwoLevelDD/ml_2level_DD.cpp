
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ************************************************************************
//@HEADER

// Goal of this example is to present how to define domain decomposition
// preconditioner using class ML_Epetra::MultiLevelPreconditioner.
//
// \author Marzio Sala, SNL 9214
// \date Last modified on 17-Nov-04

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), required Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-triutils (for the definition of the linear systems)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  // Create the linear problem using the class `Trilinos_Util::CrsMatrixGallery.'
  // The matrix is here symmetric; however, we will create a
  // non-symmetric preconditioner (symmetric preconditioner can be
  // created as well with minor minofications)

  CrsMatrixGallery Gallery("laplace_3d", Comm);
  Gallery.Set("problem_size", 8000);
  
  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  int options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE];
  
  // SetDefaults() will call AZ_defaults(options,params), and will also set the
  // preconditioner as `AZ_dom_decomp'. We will overwrite some values later in
  // this file.
  // NOTE THAT THE VECTORS ARE NOT COPIED! Only the pointer is copied, 
  // so do not destroy options and params
  // before the end of the linear system solution!
  //
  // You can also call SetDefaults() without passing `options' and `params.' This
  // way, the code will allocate a int and a double vector, that must be freed by
  // the user.

  ML_Epetra::SetDefaults("DD",MLList,options,params);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // Some parameters are reported here to better explain the process
  // even if they are as defaults. 
  // NOTE: To use `METIS' as aggregation scheme, you need to configure
  // ML with the option --with-ml_metis. Otherwise, the code will
  // creates aggregates containing all the local nodes (that is,
  // the dimension of the coarse problem will be equal to the
  // number of processors)
 
  MLList.set("aggregation: type", "METIS");
  MLList.set("smoother: type","Aztec");
  
  // put 64 nodes on each aggregate. This number can be too small
  // for large problems. In this case, either augment ir, or increase
  // the number of levels. Also, use only presmoothing, and KLU as
  // coarse solver (KLU is enabled by default with Amesos)

  MLList.set("aggregation: nodes per aggregate", 64);
  MLList.set("smoother: pre or post", "pre");
  MLList.set("coarse: type","Amesos-KLU");
  
  // now we need to define the solvers on each subdomain (== processor).
  // Here we use an incomplete Cholesky factorization, with no fill-in
  // and no overlap. To that aim, we use Aztec's preconditioning function.
  // Aztec requires two more vectors. Note: the following options and params
  // will be used ONLY for the smoother, and will NOT affect the Aztec solver

  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_icc;

  // create the preconditioning object. We suggest to use `new' and
  // `delete' because the destructor contains some calls to MPI (as
  // required by ML and possibly Amesos). This is an issue only if the
  // destructor is called **after** MPI_Finalize().
 
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 as tolerance on the
  // relative residual  
  solver.Iterate(500, 1e-12);

  // delete the preconditioner. Do it BEFORE MPI_Finalize
  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

  if (residual > 1e-5)
    exit(EXIT_FAILURE);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
  
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO) */
