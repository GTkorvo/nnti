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
#include "EpetraExt_MultiVectorIn.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToMultiVector( const char *filename, const Epetra_BlockMap & map, Epetra_MultiVector * & A) {

  const int lineLength = 1025;
  const int tokenLength = 35;
  char line[lineLength];
  char token[tokenLength];
  char token1[tokenLength];
  char token2[tokenLength];
  char token3[tokenLength];
  char token4[tokenLength];
  char token5[tokenLength];
  int M, N;

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file

  // Check first line, which should be "%%MatrixMarket matrix coordinate real general" (without quotes)
  if(fgets(line, lineLength, handle)==0) return(-1);
  if(sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5 )==0) return(-1);
  if (strcmp(token1, "%%MatrixMarket") ||
      strcmp(token2, "matrix") ||
      strcmp(token3, "array") ||
      strcmp(token4, "real") ||
      strcmp(token5, "general")) return(-1);

  // Next, strip off header lines (which start with "%")
  do {
    if(fgets(line, lineLength, handle)==0) return(-1);
  } while (line[0] == '%');

  // Next get problem dimensions: M, N
  if(sscanf(line, "%d %d", &M, &N)==0) return(-1);

  // Compute the offset for each processor for when it should start storing values
  int numMyPoints = map.NumMyPoints();
  int offset;
  map.Comm().ScanSum(&numMyPoints, &offset, 1); // ScanSum will compute offsets for us
  offset -= numMyPoints; // readjust for my PE

  // Now construct vector/multivector
  if (N==1)
    A = new Epetra_Vector(map);
  else
    A = new Epetra_MultiVector(map, N);

  double ** Ap = A->Pointers();

  for (int j=0; j<N; j++) {
    double * v = Ap[j];

    // Now read in lines that we will discard
    for (int i=0; i<offset; i++)
      if(fgets(line, lineLength, handle)==0) return(-1);
    
    // Now read in each value and store to the local portion of the the  if the row is owned.
    double V;
    for (int i=0; i<numMyPoints; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
      if(sscanf(line, "%lg\n", &V)==0) return(-1);
      v[i] = V;
    }
  }

  return(0);
}
} // namespace EpetraExt
