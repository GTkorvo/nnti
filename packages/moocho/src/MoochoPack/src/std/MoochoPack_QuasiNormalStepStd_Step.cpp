// ////////////////////////////////////////////////////////////////////////////
// RangeSpaceStepStd_Step.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>

#include "ReducedSpaceSQPPack/include/std/RangeSpaceStepStd_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"

namespace ReducedSpaceSQPPack {

bool RangeSpaceStepStd_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	using BLAS_Cpp::no_trans;
	using AbstractLinAlgPack::V_InvMtV;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_MtV;

	rSQPAlgo         &algo        = rsqp_algo(_algo);
	rSQPState        &s           = algo.rsqp_state();
	const Range1D    con_decomp   = s.decomp_sys().con_decomp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// Get iteration quantities
	IterQuantityAccess<VectorWithOpMutable>
		&c_iq   = s.c(),
		&py_iq  = s.py(),
		&Ypy_iq = s.Ypy();
	IterQuantityAccess<MatrixWithOpNonsingular>
		&R_iq = s.R();
	IterQuantityAccess<MatrixWithOp>
		&Y_iq = s.Y();

	// Solve the system py = - inv(R) * c(con_decomp)
	VectorWithOpMutable &py_k = py_iq.set_k(0);
	V_InvMtV( &py_k, R_iq.get_k(0), no_trans, *c_iq.get_k(0).sub_view(con_decomp) );
	Vt_S( &py_k, -1.0 );

	// Ypy = Y * py
	V_MtV( &Ypy_iq.set_k(0), Y_iq.get_k(0), no_trans, py_k );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||py||   = " << py_iq.get_k(0).norm_inf() << std::endl
			<< "\n||Ypy||2 = " << Ypy_iq.get_k(0).norm_2()  << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out	<< "\npy_k =\n"  << py_iq.get_k(0);
		out	<< "\nYpy_k =\n" << py_iq.get_k(0);
	}

	return true;
}

void RangeSpaceStepStd_Step::print_step(
	const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Calculate the range space step\n"
		<< L << "py_k = - inv(R_k) * c_k(con_decomp)\n"
		<< L << "Ypy_k = Y_k * py_k\n";
}

} // end namespace ReducedSpaceSQPPack
