// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_H

#include "MoochoPack_ReducedHessianSecantUpdate_Strategy.hpp"
#include "MoochoPack_BFGSUpdate_Strategy.hpp"
#include "MoochoPack_quasi_newton_stats.hpp"
#include "MoochoPack_act_set_stats.hpp"
#include "ConstrainedOptPack_MatrixHessianSuperBasic.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Perform BFGS updates on only the free independent (super basic) variables.
 *
 * This method should be very efficient for few super basic variables.
 */
class ReducedHessianSecantUpdateBFGSProjected_Strategy : public ReducedHessianSecantUpdate_Strategy
{
public:
  
  /** \brief <<std comp>> members for the strategy object that will
   * perform the guts of the BFGS update.
   */
  STANDARD_COMPOSITION_MEMBERS( BFGSUpdate_Strategy, bfgs_update )

  /** \brief Set the ratio of the number of inequality constraints in the
    * active-set of the last two calls before a projected updating
    * for superbasic variables only is started.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, act_set_frac_proj_start )

  /** \brief Set the tolerance for determining if a projected BFGS update is
   * valid ???
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, project_error_tol )

  /** \brief Set the tolerance for Langrange multipliers for fixed variables
   * below which rows/cols from rHL_RR will not be dropped.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, super_basic_mult_drop_tol )
    
  /** \brief . */
    ReducedHessianSecantUpdateBFGSProjected_Strategy(
    const bfgs_update_ptr_t&      bfgs_update                = NULL
    ,value_type                   act_set_frac_proj_start    = 0.8
    ,value_type                   project_error_tol          = 1e-5
    ,value_type                   super_basic_mult_drop_tol  = 1e-5
    );      

  /** \brief . */
  bool perform_update(
    DVectorSlice* s_bfgs, DVectorSlice* y_bfgs, bool first_update
    ,std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,MatrixOp *rHL_k
    );
  /** \brief . */
  void print_step( std::ostream& out, const std::string& leading_str ) const;

private:
  
  // /////////////////////////////
  // Private data members

  quasi_newton_stats_iq_member    quasi_newton_stats_;
  act_set_stats_iq_member         act_set_stats_;

}; // end class ReducedHessianSecantUpdateBFGSProjected_Strategy

}  // end namespace MoochoPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_H
