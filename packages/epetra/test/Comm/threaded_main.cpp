//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
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

// Petra_Comm Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#ifdef PETRA_MPI
#include <mpi.h>
#endif
#include "Petra_Comm.h"
#include "Petra_Time.h"

int main(int argc, char *argv[]) {

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // I'm alive !!!



  Petra_Comm & petracomm = *new Petra_Comm( MPI_COMM_WORLD );
  int MyPID =  petracomm.MyPID();
  int NumProc =  petracomm.NumProc();
  cout << "Processor "<<MyPID<<" of " << NumProc << " is alive."<<endl;

  if (NumProc!=2) {
    cout << " This special routine only works for 2 processors " << endl;
    abort();
  }


  Petra_Comm * other_comm;
  MPI_Status status;

  unsigned int icomm = &petracomm;

  if (MyPID==1) cout << "Address of Petra_Comm object on PE 1 = " << &petracomm << endl;

  if (MyPID==0) MPI_Recv((void *) &other_comm, 1, MPI_UNSIGNED, 1, 99, MPI_COMM_WORLD, &status);
  else MPI_Send ( (void *) &icomm, 1, MPI_UNSIGNED, 0, 99, MPI_COMM_WORLD);

  if (MyPID==0) cout << "Address of other Petra_Comm object on PE 0 = " << other_comm << endl;
  
  int otherPID = other_comm->MyPID();

  if (MyPID==0) cout << "Processor "<<MyPID<<" of " << NumProc
		     << " has a neighbor processor with ID "
		     << otherPID << " of " << other_comm->NumProc() <<endl;
 
  delete &petracomm;
  MPI_Finalize();

  return 0;
}

