/*@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
//@HEADER
*/

#include "Epetra_ConfigDefs.h"
#include "EpetraExt_Version.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_DataAccess.h"

#include "EpetraExt_MatlabEngine.h"

#define BUFSIZE 200
#define MATLABBUF 1024 * 16

int main(int argc, char *argv[]) {
cout << "going to setup MPI...\n";

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  cout << "mpit setup complete\n";

  int MyPID = comm.MyPID();

  char s [BUFSIZE] ;
  char matlabBuffer [MATLABBUF];
  cout << "going to init matlab\n";
  EpetraExt::MatlabEngine engine (comm);
  cout << "matlab starten";
  /* MultiVector test
  cout << MyPID << " going to do multivector test...\n";
  int numGlobalElements = 100;
  int M = numGlobalElements/comm.NumProc();
  int N = 3;
  int numMyElements = M * N;
  double* A = new double[numMyElements];
  double* Aptr = A;
  int startValue = 0;

  cout << MyPID << " allocated space for A, now filling A\n";
  for(int col=0; col < N; col++) {
	startValue = (col * numGlobalElements) + (M * MyPID);
	for(int row=0; row < M; row++) {
          *Aptr++ = row+startValue;
      }
  }
  cout << MyPID << " A filled\n";

  Epetra_Map map (numGlobalElements, 0, comm);
  Epetra_MultiVector multiVector (Copy, map, A, M, N);
  //cout << multiVector;
  engine.PutMultiVector(multiVector, "TEST");
  */
  
  /*SerialDenseMatrix test
  cout << MyPID << " going to do SerialDenseMatrix test...\n";
  double* A = new double[30];
  cout << MyPID << " allocated space for A, now filling A\n";
  double* Aptr = A;
  int M = 5;
  int N = 6;
  int startValue = M*N*comm.MyPID();
  for(int i=0; i < M*N; i++) {
      *Aptr++ = i + startValue;
  }
  cout << MyPID << " A filled\n";

  Epetra_SerialDenseMatrix sdMatrix (View, A, M, M, N);
  engine.PutSerialDenseMatrix(sdMatrix, "TEST", 0);
  cout << sdMatrix;
  */

  /* SerialDenseVector test
  double* A = new double[30];
  double* Aptr = A;
  int length = 30;
  for(int i=0; i < length; i++) {
      *Aptr++ = i;
  }

  Epetra_SerialDenseVector sdVector (Copy, A, length);
  engine.PutSerialDenseMatrix(sdVector, "SDVECTOR");
  cout << sdVector;
  */

  ///*IntSerialDenseMatrix test
  cout << MyPID << " going to do IntSerialDenseMatrix test...\n";
  int* A = new int[30];
  cout << MyPID << " allocated space for A, now filling A\n";
  int* Aptr = A;
  int M = 5;
  int N = 6;
  for(int i=0; i < M*N; i++) {
      *Aptr++ = i;
  }
  cout << MyPID << " A filled\n";

  Epetra_IntSerialDenseMatrix isdMatrix (View, A, M, M, N);
  engine.PutIntSerialDenseMatrix(isdMatrix, "TEST", 0);
  cout << isdMatrix;
  //*/


  /* SerialDenseVector test
  int* A = new int[30];
  int* Aptr = A;
  int length = 30;
  for(int i=0; i < length; i++) {
      *Aptr++ = i;
  }

  Epetra_IntSerialDenseVector isdVector (Copy, A, length);
  engine.PutIntSerialDenseMatrix(isdVector, "ISDVECTOR");
  cout << isdVector;
  */

  /*while(1) {

	// do nothing
	}*/

  /*if (comm.NumProc() == 1) {
  int err;
  while(1) {
      // Prompt the user and get a string
      printf(">> ");
      if (fgets(s, BUFSIZE, stdin) == NULL) {
          printf("Bye\n");
          break ;
      }
      printf ("command :%s:\n", s) ;
      
      // Send the command to MATLAB
      // output goes to stdout
      err = engine.EvalString(s, matlabBuffer, MATLABBUF);
      if (err != 0) {
          printf("there was an error: %d", err);
		  err = 0;
      }
      else {
      	  printf("Matlab Output:\n%s", matlabBuffer);
      }
  }
  }*/

  //delete engine ;

  
  engine.EvalString("size(TEST)", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";
  engine.EvalString("TEST", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";

  cout << "\nall done\n";
  return(0);
  

  //Epetra_Map * map;
  //Epetra_CrsMatrix * A; 
  //Epetra_Vector * x; 
  //Epetra_Vector * b;
  //Epetra_Vector * xexact;

/*  int nx = 20*comm.NumProc();
  int ny = 30;
  int npoints = 7;
  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

   
  // Call routine to read in HB problem
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact);

  double residual;
  residual = A->NormInf(); double rAInf = residual;
  if (verbose) cout << "Inf Norm of A                                                     = " << residual << endl;
  residual = A->NormOne(); double rAOne = residual;
  if (verbose) cout << "One Norm of A                                                     = " << residual << endl;
  xexact->Norm2(&residual); double rxx = residual;
  if (verbose) cout << "Norm of xexact                                                    = " << residual << endl;
  Epetra_Vector tmp1(*map);
  A->Multiply(false, *xexact, tmp1);
  tmp1.Norm2(&residual); double rAx = residual;
  if (verbose) cout << "Norm of Ax                                                        = " << residual << endl;
  b->Norm2(&residual); double rb = residual;
  if (verbose) cout << "Norm of b (should equal norm of Ax)                               = " << residual << endl;
  tmp1.Update(1.0, *b, -1.0);
  tmp1.Norm2(&residual);
  if (verbose) cout << "Norm of difference between compute Ax and Ax from file            = " << residual << endl;
  comm.Barrier();

  EpetraExt::BlockMapToMatrixMarketFile("Test_map.mm", *map, "Official EpetraExt test map", 
				 "This is the official EpetraExt test map generated by the EpetraExt regression tests");

  EpetraExt::RowMatrixToMatrixMarketFile("Test_A.mm", *A, "Official EpetraExt test matrix", 
				 "This is the official EpetraExt test matrix generated by the EpetraExt regression tests");

  EpetraExt::VectorToMatrixMarketFile("Test_x.mm", *x, "Official EpetraExt test initial guess", 
				 "This is the official EpetraExt test initial guess generated by the EpetraExt regression tests");
				       
  EpetraExt::VectorToMatrixMarketFile("Test_xexact.mm", *xexact, "Official EpetraExt test exact solution", 
				 "This is the official EpetraExt test exact solution generated by the EpetraExt regression tests");
				       
  EpetraExt::VectorToMatrixMarketFile("Test_b.mm", *b, "Official EpetraExt test right hand side", 
				 "This is the official EpetraExt test right hand side generated by the EpetraExt regression tests");
				       
  Epetra_Map * map1;
  Epetra_CrsMatrix * A1; 
  Epetra_Vector * x1; 
  Epetra_Vector * b1;
  Epetra_Vector * xexact1;

  EpetraExt::MatrixMarketFileToMap("Test_map.mm", comm, map1);

  if (map->SameAs(*map1))
    if (verbose) cout << "Maps are equal.  In/Out works." << endl;
  else
    if (verbose) cout << "Maps are not equal.  In/Out fails." << endl;

  EpetraExt::MatrixMarketFileToCrsMatrix("Test_A.mm", *map1, A1);
  EpetraExt::MatrixMarketFileToVector("Test_x.mm", *map1, x1);
  EpetraExt::MatrixMarketFileToVector("Test_xexact.mm", *map1, xexact1);
  EpetraExt::MatrixMarketFileToVector("Test_b.mm", *map1, b1);


  residual = A1->NormInf(); double rA1Inf = residual;
  if (verbose) cout << "Inf Norm of A1                                                    = " << residual << endl;
  residual = A1->NormOne(); double rA1One = residual;
  if (verbose) cout << "One Norm of A1                                                    = " << residual << endl;
  xexact1->Norm2(&residual); double rxx1 = residual;
  if (verbose) cout << "Norm of xexact1                                                   = " << residual << endl;
  Epetra_Vector tmp11(*map1);
  A->Multiply(false, *xexact1, tmp11);
  tmp11.Norm2(&residual); double rAx1 = residual;
  if (verbose) cout << "Norm of A1*x1                                                     = " << residual << endl;
  b1->Norm2(&residual); double rb1 = residual;
  if (verbose) cout << "Norm of b1 (should equal norm of Ax)                              = " << residual << endl;
  tmp11.Update(1.0, *b1, -1.0);
  tmp11.Norm2(&residual);
  if (verbose) cout << "Norm of difference between compute A1x1 and A1x1 from file        = " << residual << endl;

  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

  delete A1;
  delete x1;
  delete b1;
  delete xexact1;
  delete map1;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
*/
}
