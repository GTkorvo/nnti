/*
// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

// Copy this file over any example that is using Triutils

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  CrsMatrixGallery Gallery("laplace_1d", Comm);
  Gallery.Set("problem_size",1000000 * Comm.NumProc());

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  assert (Problem != 0);
  // retrive pointers to solution (lhs), right-hand side (rhs)
  // and matrix itself (A)
  Epetra_MultiVector* lhs = Problem->GetLHS();
  Epetra_MultiVector* rhs = Problem->GetRHS();
  Epetra_RowMatrix* A = Problem->GetMatrix();
    
  Epetra_Time Time(Comm);

  for (int i = 0 ; i < 100 ; ++i)
    A->Multiply(false, *lhs, *rhs);

  cout << Time.ElapsedTime() << endl;

  return(EXIT_SUCCESS);
} // end of main()
