/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef PLAYA_INVERSEOPERATOR_DECL_HPP
#define PLAYA_INVERSEOPERATOR_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "Teuchos_RCP.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaSolverState.hpp"

namespace Playa
{
using Teuchos::RCP;

/** 
 * PlayaInverseOperator represents the inverse of some other operator.  An
 * inverse operator object will contain an operator and a solver.  The 
 * operator data member is the operator whose inverse this represents.  The
 * solver data member is the solver that will be used in applying the
 * inverse.  If the solver is null, the operator is assumed to have
 * self-contained ability to solve systems, as for a dense matrix that 
 * does solves by factoring and backsolves.
 */
template <class Scalar> 
class InverseOperator : public LinearOpWithSpaces<Scalar>,
                        public Printable
{
public:
  /**
   * Ctor with a linear operator and a solver specified.
   */
  InverseOperator(const LinearOperator<Scalar>& op, 
    const LinearSolver<Scalar>& solver);

  /** Virtual dtor */
  virtual ~InverseOperator(){;}

  /** 
   * Apply the operator. 
   * 
   * \param applyType Indicates whether to apply the operator, its transpose,
   * or its conjugate transpose. 
   * \param in The vector on which the operator is to act
   * \param out The vector into which the result of the operation 
   * is to be written. This vector should already be initialized by the
   * appropriate space.
   **/
  virtual void apply(
    Teuchos::ETransp applyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const ;


  
  /** */
  void print(std::ostream& os) const ;

  /** */
  LinearOperator<Scalar> op() const {return op_;}


private:
  const LinearOperator<Scalar> op_;
  const LinearSolver<Scalar> solver_;  
  std::string msg_;
};


/** \brief Implicit inverse operator. */
template <class Scalar> 
LinearOperator<Scalar> 
inverse(const LinearOperator<Scalar>& op, 
  const LinearSolver<Scalar>& solver);
  

}

#endif
