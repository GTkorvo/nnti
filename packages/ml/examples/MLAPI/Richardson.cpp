
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

#include "ml_config.h"
#include "ml_common.h"

#if defined(HAVE_ML_MLAPI)

#include "MLAPI.h"

using namespace Teuchos;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // Initialize the workspace and set the output level
    
    Init();

    // global dimension of the problem
    
    int NumGlobalElements = 10000;

    // define the space for fine level vectors and operators.
    
    Space S(NumGlobalElements);

    // define the linear system matrix.
    
    Operator A = Gallery("laplace_2d", S);

    // set parameters for aggregation and smoothers
    // NOTE: only a limited subset of the parameters accepted by
    // class ML_Epetra::MultiLevelPreconditioner is supported
    // by MLAPI::MultiLevelSA
    
    Teuchos::ParameterList MLList;
    MLList.set("max levels",3);
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("aggregation: damping factor", 0.0);
    MLList.set("smoother: type","symmetric Gauss-Seidel");
    MLList.set("smoother: sweeps",1);
    MLList.set("smoother: damping factor",1.0);
    MLList.set("coarse: max size",3);
    MLList.set("coarse: type","Amesos-KLU");

    MultiLevelSA P(A, MLList);

    // Here we define a simple Richardson method for the
    // solution of A x = b. The preconditioner is P,
    // the exact solution (x_ex) is a random vector, the
    // starting solution (x) is the zero vector. 
    
    MultiVector x_ex(S);
    MultiVector x(S);
    MultiVector b(S);
    MultiVector r(S);
    MultiVector z(S);

    x_ex.Random();
    b = A * x_ex;
    x = 0.0;
    
    double OldNorm   = 1.0;
    double Tolerance = 1e-13;
    int    MaxIters  = 30;

    // ================ //
    // Richardson cycle //
    // ================ //

    for (int i = 0 ; i < MaxIters ; ++i) {

      r = b - A * x; // new residual
      z = P * r;     // apply preconditioner with zero initial guess
      x = x + z;     // update solution

      // compute the A-norm of the error

      double NewNorm = sqrt((x - x_ex) * (A * (x - x_ex)));

      if (GetMyPID() == 0 && i) {
        cout << "||x - x_ex||_A = ";
        cout.width(15);
        cout << NewNorm << ", ";
        cout << "reduction = ";
        cout.width(15);
        cout << NewNorm / OldNorm << endl;
      }

      if (NewNorm < Tolerance)
        break;

      OldNorm = NewNorm;

    }

    // finalize the MLAPI workspace
    
    Finalize();

  }
  catch (const int e) {
    cout << "Caught integer exception, code = " << e << endl;
  } 
  catch (...) {
    cout << "problems here..." << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);

}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("The ML API requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // if defined(HAVE_ML_MLAPI)
