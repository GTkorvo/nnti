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

#ifndef TSFBLOCKOPERATORDECL_HPP
#define TSFBLOCKOPERATORDECL_HPP

#include "TSFConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "TSFExplicitlyTransposeableOp.hpp"
#include "TSFRowAccessibleOp.hpp"
#include "TSFHandleable.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"




namespace TSFExtended
{
template <class Scalar>
class LinearOperator;

using Teuchos::RefCountPtr;
using Teuchos::Array;

/**
 * Class BlockOperator provides an abstract interface for configuration
 * and building block operators
 *
 * @author Paul T Boggs (ptboggs@sandia.gov)
 */
template <class Scalar>
class BlockOperator : public Thyra::DefaultBlockedLinearOp<Scalar>,
                      public RowAccessibleOp<Scalar>,
                      public Printable
{
public:

  /** ctor with domain and range specified.  The blocks must be
   *	specified later and all filled before use.
   */
  BlockOperator(const VectorSpace<Scalar> &domain, 
    const VectorSpace<Scalar> &range);

  /** */
  RefCountPtr<const VectorSpaceBase<Scalar> > range() const 
    {return range_.ptr();}

  /** */
  RefCountPtr<const VectorSpaceBase<Scalar> > domain() const 
    {return domain_.ptr();}



  /** */
  int numBlockRows() const {return this->productRange()->numBlocks();}

  /** */
  int numBlockCols() const {return this->productDomain()->numBlocks();}
    

 

  /** */
  void apply(
    const EConj                             conj
    ,const Thyra::MultiVectorBase<Scalar>    &X
    ,Thyra::MultiVectorBase<Scalar>           *Y
    ,const Scalar                      alpha = Teuchos::ScalarTraits<Scalar>::one()
    ,const Scalar                      beta  = Teuchos::ScalarTraits<Scalar>::zero()
    ) const ;

  /** */
  void applyTranspose(
    const EConj                            conj
    ,const Thyra::MultiVectorBase<Scalar>    &X
    ,Thyra::MultiVectorBase<Scalar>         *Y
    ,const Scalar                     alpha = Teuchos::ScalarTraits<Scalar>::one()
    ,const Scalar                     beta  = Teuchos::ScalarTraits<Scalar>::zero()
    ) const;

  /** Get entire row of the block matrix  */
  void getRow(const int& row, Array<int>& indices,
		Array<Scalar>& values) const;

    

  /** */
  std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const;

  /** */
  void print(std::ostream& os) const ;
private:
  VectorSpace<Scalar> domain_;
  VectorSpace<Scalar> range_;
}; 
}

#endif
