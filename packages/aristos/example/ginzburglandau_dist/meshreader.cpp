#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include <iostream>

int meshreader(const Epetra_Comm & Comm,
               Epetra_IntSerialDenseVector & ipindx,
               Epetra_SerialDenseMatrix & ipcoords,
               Epetra_IntSerialDenseVector & pindx,
               Epetra_SerialDenseMatrix & pcoords,
               Epetra_IntSerialDenseMatrix & t,
               Epetra_IntSerialDenseMatrix & e,
               const char * geomfile)
{
  int MyPID = Comm.MyPID();
  int numip = 0, numcp = 0, nump = 0, numelems = 0, numedges = 0;
  char FileName[120];
  FILE* fid;

  sprintf(FileName, "%s.%03d", geomfile, MyPID);
  fid = fopen(FileName, "r");

  printf("Reading node info from: %s ...\n", FileName);
  fscanf(fid, "%d %d", &numip, &numcp);
  nump = numip + numcp;
  ipindx.Size(numip);
  ipcoords.Shape(numip, 2);
  pindx.Size(nump);
  pcoords.Shape(nump, 2);
  for (int i=0; i<numip; i++) {
    fscanf(fid, "%d %lf %lf %*d", &ipindx(i), &ipcoords(i,0),
           &ipcoords(i,1));
    pindx(i) = ipindx(i);
    pcoords(i,0) = ipcoords(i,0); pcoords(i,1) = ipcoords(i,1);
  }
  for (int i=numip; i<nump; i++) {
    fscanf(fid, "%d %lf %lf %*d", &pindx(i), &pcoords(i,0),
           &pcoords(i,1));
  }

  fscanf(fid, "%*[^\n]");   /* Skip to the End of the Line */
  fscanf(fid, "%*1[\n]");   /* Skip One Newline */

  fscanf(fid, "%*[^\n]");   /* Skip to the End of the Line */
  fscanf(fid, "%*1[\n]");   /* Skip One Newline */

  for (int i=0; i<nump; i++) {
    fscanf(fid, "%*[^\n]");   /* Skip to the End of the Line */
    fscanf(fid, "%*1[\n]");   /* Skip One Newline */
  }

  printf("Reading element info from: %s ...\n", FileName);
  fscanf(fid, "%d", &numelems);
  t.Shape(numelems, 3);
  for (int i=0; i<numelems; i++) {
    fscanf(fid, "%d %d %d", &t(i,0), &t(i,1), &t(i,2));
  }

  printf("Reading edge info from: %s ...\n", FileName);
  fscanf(fid, "%d", &numedges);
  e.Shape(numedges, 3);
  for (int i=0; i<numedges; i++) {
    fscanf(fid, "%d %d %d", &e(i,0), &e(i,1), &e(i,2));
  }

  fclose(fid);
  printf(" Done.\n");


  return(0);

}
