// @HEADER
// ***********************************************************************
// 
//                IFPACK
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_TEUCHOS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_Preconditioner.h"
#include "Ifpack.h"
#include "AztecOO.h"
#include "Ifpack_CrsIct.h"
#include "Ifpack_DropFilter.h"
#ifdef HAVE_IFPACK_AMESOS
#include "Amesos_TestRowMatrix.h"
#endif

using namespace Trilinos_Util;

// =======================================================================
// GOAL: test that the names in the factory do not change. This test
//       will not solve any linear system.
//
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // size of the global matrix. 
  const int NumPoints = 9;

  CrsMatrixGallery Gallery("minij", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "interlaced");
  Epetra_RowMatrix* A = Gallery.GetMatrix();

  Ifpack Factory;
  Ifpack_Preconditioner* Prec;

  Prec = Factory.Create("point relaxation", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("point relaxation stand-alone", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("block relaxation", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("block relaxation stand-alone", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("IC", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("ICT", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("ILU", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("ILUT", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

#ifdef HAVE_IFPACK_AMESOS
  Prec = Factory.Create("Amesos", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;

  Prec = Factory.Create("Amesos stand-alone", A);
  assert (Prec != 0);
  Prec->Initialize();
  Prec->Compute();
  cout << *Prec;
  delete Prec;
#endif
  
#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  if (Comm.MyPID() == 0)
    cout << "Test `PrecondititonerFactory.exe' passed!" << endl;
  return(EXIT_SUCCESS);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --enable-teuchos to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
