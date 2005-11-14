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
#include "TSFSingleScalarTypeOp.hpp"
#include "Thyra_DefaultBlockOperatorDecl.hpp"




namespace TSFExtended
{
  template <class Scalar>
  class LinearOperator;

  using namespace Teuchos;

  /**
   * Class BlockOperator provides an abstract interface for configuration
   * and building block operators
   *
   * @author Paul T Boggs (ptboggs@sandia.gov)
   */
  template <class Scalar>
  class BlockOperator : public SingleScalarTypeOpBase<Scalar>,
                        public Handleable<SingleScalarTypeOpBase<Scalar> >,
                        public RowAccessibleOp<Scalar>,
                        public Thyra::DefaultBlockOperator<Scalar, Scalar>,
                        public Printable
  {
  public:
    GET_RCP(SingleScalarTypeOpBase<Scalar>);

    /** ctor with domain and range specified.  The blocks must be
     *	specified later and all filled before use.
     */
    BlockOperator(const VectorSpace<Scalar> &domain, 
		  const VectorSpace<Scalar> &range);

     /** \brief Return a smart pointer for the range space for <tt>this</tt> operator.
   */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > range() const ;

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<Scalar> > domain() const ;

    /** 
     * Compute alpha*M*x + beta*y, where M=*this.
     * @param M_trans specifies whether the operator is transposed:
     *                op(M) = M, for M_trans == NOTRANS
     *                op(M) = M', for M_trans == TRANS
     * @param x       vector of length this->domain()->dim()
     * @param y       vector of length this->range()->dim()
     * @param alpha   scalar multiplying M*x (default is 1.0)
     * @param beta    scalar multiplying y (default is 0.0)
     */
    virtual void generalApply(
                       const Thyra::ETransp            M_trans
                       ,const Thyra::VectorBase<Scalar>    &x
                       ,Thyra::VectorBase<Scalar>          *y
                       ,const Scalar            //alpha = 1.0
                       ,const Scalar           // beta  = 0.0
                       ) const;

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
    void getRow(const int& row, Teuchos::Array<int>& indices,
		Teuchos::Array<Scalar>& values) const;

    

    /** */
    std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const;

    /** */
    void print(ostream& os) const ;
  private:

  }; 
}

#endif
