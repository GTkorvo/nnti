// ////////////////////////////////////////////////////////////////////////////
// CalcReducedGradLagrangianStd_AddedStep.h
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
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

#ifndef CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H
#define CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H

#include "../rSQPAlgo_Step.h"

namespace ReducedSpaceSQPPack {

///
/** Calculates the reduced gradient of the Lagrangian.
  * rGL = rGf + Z' * nu + V' * lambda(dep)
  */
class CalcReducedGradLagrangianStd_AddedStep : public rSQPAlgo_Step {
public:

	// ////////////////////
	// Overridden

	///
	bool do_step(Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CalcReducedGradLagrangianStd_AddedStep

}	// end namespace ReducedSpaceSQPPack 

#endif	// CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H
