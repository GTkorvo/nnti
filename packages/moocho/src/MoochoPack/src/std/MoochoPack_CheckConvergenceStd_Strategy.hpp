// ////////////////////////////////////////////////////////////////
// ReducedSpaceSQPPack/src/std/CheckConvergenceStd_Strategy.h
//
// Copyright (C) 2001
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//

#if !defined CHECK_CONVERGENCE_STD_STRATEGY_H
#define CHECK_CONVERGENCE_STD_STRATEGY_H

#include "ReducedSpaceSQPPack/src/std/CheckConvergence_Strategy.h"

namespace ReducedSpaceSQPPack {

///
/** Implementation of CheckConvergence_Strategy interface
 *
 * This object can not change the flow of control or do anything fancy.  It just
 *  checks convergence by calculating norm errors and comparing with tolerance
 *  It can update iteration quantities if desired.
 *
 * See the printed documentation generated by \c this->print_step().
 */
class CheckConvergenceStd_Strategy :
		public CheckConvergence_Strategy
	{
	public:
		CheckConvergenceStd_Strategy(
		  EOptErrorCheck opt_error_check = OPT_ERROR_REDUCED_GRADIENT_LAGR,
		  EScaleKKTErrorBy scale_opt_error_by = SCALE_BY_ONE,
		  EScaleKKTErrorBy scale_feas_error_by = SCALE_BY_ONE,
		  EScaleKKTErrorBy scale_comp_error_by = SCALE_BY_ONE,
		  bool scale_opt_error_by_Gf = true
		  );
		
		/** @name Overridden from CheckConvergence_Strategy */
		//@{
		///
		virtual bool Converged( Algorithm& _algo);

		///
		virtual void print_step( const Algorithm& _algo, std::ostream& out, const std::string& L ) const;

		//@}

	protected:

		value_type CalculateScalingFactor( rSQPState& state, EScaleKKTErrorBy scale_by ) const;

	}; // end class CheckConvergenceStd_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // CHECK_CONVERGENCE_STD_STRATEGY_H

