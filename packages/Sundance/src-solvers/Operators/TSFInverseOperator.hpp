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

#ifndef TSFINVERSEOPERATOR_HPP
#define TSFINVERSEOPERATOR_HPP

#include "SundanceDefs.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "Thyra_ZeroLinearOpBase.hpp"
#include "TSFOpWithBackwardsCompatibleApply.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "TSFSolverState.hpp"

namespace TSFExtended
{
  /** 
   * TSFInverseOperator represents the inverse of some other operator.  An
   * inverse operator object will contain an operator and a solver.  The 
   * operator data member is the operator whose inverse this represents.  The
   * solver data member is the solver that will be used in applying the
   * inverse.  If the solver is null, the operator is assumed to have
   * self-contained ability to solve systems, as for a dense matrix that 
   * does solves by factoring and backsolves.
   */
  template <class Scalar> 
  class InverseOperator : public OpWithBackwardsCompatibleApply<Scalar>
  {
  public:


    /**
     * Ctor with a linear operator and a solver specified.
     */
    InverseOperator(const LinearOperator<Scalar>& op, 
      const LinearSolver<Scalar>& solver)
      : op_(op), solver_(solver) {;}


    /** Virtual dtor */
    virtual ~InverseOperator(){;}


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
                              ,const Scalar            alpha = 1.0
                              ,const Scalar            beta  = 0.0
                              ) const 
    {
      TEST_FOR_EXCEPTION(dynamic_cast<Thyra::ZeroLinearOpBase<Scalar>* >(op_.ptr().get()) != 0, std::runtime_error,
                         "InverseOperator<Scalar>::apply() called on a ZeroOperator.");
      TEST_FOR_EXCEPTION(op_.domain().dim() != op_.range().dim(), std::runtime_error,
                         "InverseOperator<Scalar>::apply() called on a non-square operator.");
      LinearOperator<Scalar> applyOp;      
      if (M_trans == Thyra::NOTRANS)
        {
          applyOp = op_;
        }
      else
        {
          applyOp = op_.transpose();
        }




      if (alpha==Teuchos::ScalarTraits<Scalar>::zero())
        {
          Vt_S(y, beta);
        }
      else
        {
          Vector<Scalar> temp = createMember(*(x.space()));
          Vector<Scalar> result;
          assign(temp.ptr().get(), x);
          SolverState<Scalar> haveSoln = solver_.solve(applyOp, temp, result);
          TEST_FOR_EXCEPTION(haveSoln.finalState() != SolveConverged, 
                             std::runtime_error,
                             "InverseOperator<Scalar>::apply() " 
                             << haveSoln.stateDescription());
          Vt_S(result.ptr().get(), alpha);
          V_StVpV(y, beta, *y, *result.ptr().get());
        }      
    }


    /** 
     * Return the domain of the operator. 
     */
    virtual Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > domain() const 
    {return op_.ptr()->domain();}
    

    /** 
     * Return the range of the operator. 
     */
    virtual Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > range() const 
    {return op_.ptr()->range();}

  private:
    const LinearOperator<Scalar> op_;
    const LinearSolver<Scalar> solver_;  
    std::string msg_;
  };
}

#endif
