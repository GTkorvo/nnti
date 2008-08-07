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


#ifndef TSFBLOCKOPERATORIMPL_HPP
#define TSFBLOCKOPERATORIMPL_HPP


#include "TSFProductVectorSpaceDecl.hpp"
#include "Thyra_DefaultBlockOperatorImpl.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFNonmemberOpHelpers.hpp"

using namespace TSFExtended;
using namespace TSFExtendedOps;
using namespace Teuchos;
using std::ostream;



namespace TSFExtended
{

/*==================================================================*/
template <class Scalar> inline
BlockOperator<Scalar>::BlockOperator(const VectorSpace<Scalar>& domain,
  const  VectorSpace<Scalar>& range)
  : Thyra::DefaultBlockedLinearOp<Scalar>(),
    domain_(domain),
    range_(range)
{
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > prodDomain
    = rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar> >(domain.ptr());
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > prodRange
    = rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar> >(range.ptr());
  beginBlockFill(prodRange, prodDomain);

  for (int i=0; i<prodRange->numBlocks(); i++)
  {
    VectorSpace<Scalar> rb = prodRange->getBlock(i);
    for (int j=0; j<prodDomain->numBlocks(); j++)
    {
      VectorSpace<Scalar> db = prodDomain->getBlock(j);
      LinearOperator<Scalar> Z = zeroOperator(db, rb);
      setNonconstBlock(i, j, Z.ptr());
    }
  }
}



/*==================================================================*/
template <class Scalar> inline
void BlockOperator<Scalar>::getRow(const int& row, 
  Teuchos::Array<int>& indices,
  Teuchos::Array<Scalar>& values) const
{
  /* find col block in which row exists and row within the block */
  int K = 0;
  int blockRow = 0;
  int rowInBlock = 0;
  for (int i = 0; i < this->numBlockRows(); i++)
  {
    int numR = -1;
    for (int j=0; j<this->numBlockCols(); j++)
    {
      if (this->getBlock(i,j).get()==0) continue;
      numR = this->getBlock(i, j)->range()->dim();
    }
    TEST_FOR_EXCEPTION(numR==-1, std::runtime_error,
      "Empty block row detected in "
      "BlockOperator<Scalar>::getRow()");
    K += numR;
    if (row < K)
    {
      blockRow = i;
      rowInBlock = row - (K - numR);
      break;
    }
  }

  /* get the row elements for each block in the row.  */
  int offset = 0;
  Teuchos::Array<int> localInd;
  Teuchos::Array<Scalar> localVal;
  for (int i = 0; i < this->numBlockCols(); i++)
  {
    if (this->getBlock(blockRow, i).get()==0) continue;
    RefCountPtr<const RowAccessibleOp<Scalar> > raOp 
      = rcp_dynamic_cast<const RowAccessibleOp<Scalar> >(this->getBlock(blockRow, i));
    raOp->getRow(rowInBlock,localInd, localVal);
    for (unsigned int j = 0; j < localInd.size(); j++)
    {
      indices.append(localInd[j] + offset);
      values.append(localVal[j]);
    }
    offset += this->getBlock(blockRow, i)->domain()->dim();
  }

}

/*==================================================================*/


template <class Scalar> inline
std::ostream& BlockOperator<Scalar>
::describe(
  std::ostream                         &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ,const std::string                   leadingIndent
  ,const std::string                   indentSpacer
  ) const
{
  return DefaultBlockOperator<Scalar, Scalar>::describe(out,
    verbLevel,
    leadingIndent,
    indentSpacer);
}

template <class Scalar>
inline void BlockOperator<Scalar>
::print(ostream& os) const
{
  os << this->description();
}


template <class Scalar>
inline LinearOperator<Scalar> makeLinearOperator(BlockOperator<Scalar>* b)
{
  RefCountPtr<LinearOpBase<Scalar> > p = rcp(b);
  return LinearOperator<Scalar>(p);
}
}

#endif
