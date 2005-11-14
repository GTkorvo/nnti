/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
 /* @HEADER@ */

#include "TSFMatrixLaplacian1D.hpp"
#include "TSFEpetraMatrix.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;

MatrixLaplacian1D::MatrixLaplacian1D(int nLocalRows, 
                                     const VectorType<double>& type)
  : OperatorBuilder<double>(nLocalRows, type), op_()
{
  int rank = MPISession::getRank();
  int nProc = MPISession::getNProc();
  RefCountPtr<MatrixFactory<double> > mFact 
    = vecType().createMatrixFactory(domain(), range());

  int lowestLocalRow = nLocalRows * rank;

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(mFact.get());
  for (int i=0; i<nLocalRows; i++)
    {
      int row = lowestLocalRow + i;
      Array<int> colIndices;
      if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
        {
          colIndices = tuple(row);
        }
      else
        {
          colIndices = tuple(row-1, row, row+1);
        }
      icmf->initializeNonzerosInRow(row, colIndices.size(),
                                    &(colIndices[0]));
    }
  icmf->finalize();
      
  op_ = mFact->createMatrix();
      
  RefCountPtr<LoadableMatrix<double> > mat = op_.matrix();

  /* fill in with the Laplacian operator */
  for (int i=0; i<nLocalRows; i++)
    {
      int row = lowestLocalRow + i;
      Array<int> colIndices;
      Array<double> colVals;
      if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
        {
          colIndices = tuple(row);
          colVals = tuple(1.0);
        }
      else
        {
          colIndices = tuple(row-1, row, row+1);
          colVals = tuple(-1.0, 2.0, -1.0);
        }
      mat->addToRow(row, colIndices.size(), 
                    &(colIndices[0]), &(colVals[0]));
    }
}
