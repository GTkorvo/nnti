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

#include "TSFEpetraMatrix.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFVectorSpace.hpp"  // changed from Impl
 //#include "TSFVectorImpl.hpp"
#include "TSFVector.hpp"
#include "TSFLinearOperator.hpp"  // changed from Impl
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFIfpackOperator.hpp"
#include "TSFGenericLeftPreconditioner.hpp"
#include "TSFGenericRightPreconditioner.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;

EpetraMatrix::EpetraMatrix(const Epetra_CrsGraph& graph,
                           const RefCountPtr<const EpetraVectorSpace>& domain,
                           const RefCountPtr<const EpetraVectorSpace>& range)
  : matrix_(rcp(new Epetra_CrsMatrix(Copy, graph))),
    range_(range),
    domain_(domain)
{}

void EpetraMatrix::generalApply(const Thyra::ETransp M_trans,
                                const Thyra::VectorBase<double>    &x,
                                Thyra::VectorBase<double>          *y,
                                const double            alpha,
                                const double            beta) const
{
  const EpetraVector* tx = dynamic_cast<const EpetraVector*>(&x);
  TEST_FOR_EXCEPTION(tx==0, runtime_error, 
                     "EpetraMatrix::generalApply() could not convert " 
                     << x.description() << " to an EpetraVector");

  EpetraVector* ty = dynamic_cast<EpetraVector*>(y);
  TEST_FOR_EXCEPTION(ty==0, runtime_error, 
                     "EpetraMatrix::generalApply() could not convert " 
                     << y->description() << " to an EpetraVector");

  const Epetra_Vector* epx = tx->epetraVec().get();
  Epetra_Vector* epy = ty->epetraVec().get();

  bool trans = M_trans==Thyra::TRANS;
  int ierr=0;
  if (beta==0.0)
    {
      ierr = matrix_->Multiply(trans, *epx, *epy);
      TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                         "EpetraMatrix::generalApply() detected ierr="
                         << ierr << " in matrix_->Multiply()");
      if (alpha != 1.0)
        {
          Thyra::Vt_S(y, alpha);
        }
    }
  else
    {
      Epetra_Vector tmp(M_trans == NOTRANS 
                        ? matrix_->OperatorRangeMap() 
                        : matrix_->OperatorDomainMap(), 
                        false);
      ierr = matrix_->Multiply(trans, *epx, tmp);
      TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                         "EpetraMatrix::generalApply() detected ierr="
                         << ierr << " in matrix_->Multiply()");
      epy->Update(alpha, tmp, beta);
      TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                         "EpetraMatrix::generalApply() detected ierr="
                         << ierr << " in epy->update()");
    }
}




void EpetraMatrix::addToRow(int globalRowIndex,
                            int nElemsToInsert,
                            const int* globalColumnIndices,
                            const double* elementValues)
{
  int ierr = crsMatrix()->SumIntoGlobalValues(globalRowIndex,
                                              nElemsToInsert,
                                              (double*) elementValues,
                                              (int*) globalColumnIndices);

  TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                     "failed to add to row " << globalRowIndex
                     << " in EpetraMatrix::addToRow() with nnz="
                     << nElemsToInsert 
                     << ". Error code was " << ierr);
}

void EpetraMatrix::addToElementBatch(int numRows, 
                                     int rowBlockSize,
                                     const int* globalRowIndices,
                                     int numColumnsPerRow,
                                     const int* globalColumnIndices,
                                     const double* values,
                                     const int* skipRow)
{
  Epetra_CrsMatrix* crs = crsMatrix();

  int numRowBlocks = numRows/rowBlockSize;
  int row = 0;

  for (int rb=0; rb<numRowBlocks; rb++)
    {
      const int* cols = globalColumnIndices + rb*numColumnsPerRow;
      for (int r=0; r<rowBlockSize; r++, row++)
        {
          if (skipRow[row]) continue;
          const double* rowVals = values + row*numColumnsPerRow;
          int ierr=crs->SumIntoGlobalValues(globalRowIndices[row], 
                                            numColumnsPerRow,
                                            (double*) rowVals,
                                            (int*) cols);
          TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                             "failed to add to row " << globalRowIndices[row]
                             << " in EpetraMatrix::addToRow() with nnz="
                             << numColumnsPerRow
                             << ". Error code was " << ierr);
        }
    }
}


void EpetraMatrix::zero()
{
  crsMatrix()->PutScalar(0.0);
}



void EpetraMatrix::getILUKPreconditioner(int fillLevels,
                                         int overlapFill,
                                         double relaxationValue,
                                         double relativeThreshold,
                                         double absoluteThreshold,
                                         LeftOrRight leftOrRight,
                                         Preconditioner<double>& rtn) const
{
  LinearOperator<double> ilu = new IfpackOperator(this, 
                                                  fillLevels,
                                                  overlapFill,
                                                  relaxationValue,
                                                  relativeThreshold,
                                                  absoluteThreshold);

  if (leftOrRight == Left)
    {
      rtn = new GenericLeftPreconditioner<double>(ilu);
    }
  else
    {
      rtn = new GenericRightPreconditioner<double>(ilu);
    }
}


void EpetraMatrix::print(ostream& os) const 
{
  int nProc = MPISession::getNProc();
  int rank = MPISession::getRank();
  for (int i=0; i<nProc; i++)
    {
      MPIComm::world().synchronize();
      if (i==rank) crsMatrix()->Print(os);
    }
}

string EpetraMatrix::description() const 
{
  string rtn = "EpetraMatrix[nRow=" 
    + Teuchos::toString(crsMatrix()->NumGlobalRows())
    + ", nCol=" + Teuchos::toString(crsMatrix()->NumGlobalCols())
    + "]";
  return rtn;
}

Epetra_CrsMatrix* EpetraMatrix::crsMatrix()
{
  return matrix_.get();
}

const Epetra_CrsMatrix* EpetraMatrix::crsMatrix() const 
{
  return matrix_.get();
}



Epetra_CrsMatrix& EpetraMatrix::getConcrete(const LinearOperator<double>& A)
{
  EpetraMatrix* ep 
    = dynamic_cast<EpetraMatrix*>(A.ptr().get());
  TEST_FOR_EXCEPTION(ep==0, std::runtime_error,
                     "EpetraMatrix::getConcrete called on a matrix that "
                     "could not be cast to an EpetraMatrix");
  return *(ep->crsMatrix());
}





void EpetraMatrix::getRow(const int& row, 
                          Teuchos::Array<int>& indices, 
                          Teuchos::Array<double>& values) const
{
  const Epetra_CrsMatrix* crs = crsMatrix();

  int numEntries;
  int* epIndices;
  double* epValues;

  int info = crs->ExtractGlobalRowView(row, numEntries, epValues, epIndices);
  TEST_FOR_EXCEPTION(info != 0, runtime_error,
                     "call to ExtractGlobalRowView not successful");

  indices.resize(numEntries);
  values.resize(numEntries);
  for (int i = 0; i < numEntries; i++)
    {
      indices[i] = *epIndices;
      values[i] = *epValues;
      epIndices++;
      epValues++;
    }
}
