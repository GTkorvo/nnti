// //////////////////////////////////////////////////////////
// ExampleNLPFirstOrderDirectRun.cpp
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
//

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "ExampleNLPFirstOrderDirectRun.hpp"
#include "ExampleNLPFirstOrderDirect.hpp"
#include "ReducedSpaceSQPPack/Configurations/MamaJama/rSQPAlgo_ConfigMamaJama.hpp"
#include "IterationPack/src/AlgorithmTracker.hpp"
#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/BasisSystem.hpp"
#include "OptionsFromStream.hpp"

ReducedSpaceSQPPack::rSQPppSolver::ESolutionStatus
NLPInterfacePack::ExampleNLPFirstOrderDirectRun(
	const VectorSpace&   vec_space
	,value_type          xo
	,bool                has_bounds
	,bool                dep_bounded
	,std::ostream*       console_out
	,std::ostream*       error_out
	,bool                throw_solve_exception
	,std::ostream*       algo_out
	,std::ostream*       summary_out
	,std::ostream*       journal_out
	)
{
	using std::endl;
	using std::setw;
	namespace rcp = MemMngPack;
	using rcp::ref_count_ptr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;
	namespace rsqp = ReducedSpaceSQPPack;
	using rsqp::rSQPppSolver;
	using rsqp::rSQPAlgo_ConfigMamaJama;

	rSQPppSolver::ESolutionStatus
		solve_return = rSQPppSolver::SOLVE_RETURN_EXCEPTION;

	int err = 0;
	
	int w = 15;
	int prec = 8;

	if(console_out)
		*console_out
			<< std::setprecision(prec)
			<< std::scientific
			<< "***************************************************\n"
			<< "*** Running Tests on ExampleNLPFirstOrderDirect ***\n"
			<< "***************************************************\n"
			<< "\nUsing a vector space of type \'" << typeid(vec_space).name() << "\'"
			<< "\nwith a dimension of vec_space.dim() = " << vec_space.dim()
			<< std::endl;

	// Create the nlp
	ExampleNLPFirstOrderDirect
		nlp(VectorSpace::space_ptr_t(&vec_space,false),xo,has_bounds,dep_bounded);

	// Create the solver object and set it up
	rSQPppSolver solver;
	solver.set_nlp(rcp::rcp(&nlp,false));                  // Set the NLP!
	solver.set_error_handling(                             // set up outputting
		throw_solve_exception
		,rcp::rcp(error_out,false)
		);
	solver.set_console_out(rcp::rcp(console_out,false));
	solver.set_summary_out(rcp::rcp(summary_out,false));
	solver.set_journal_out(rcp::rcp(journal_out,false));
	solver.set_algo_out(   rcp::rcp(algo_out,false)   );

	// Run rSQP++ using the MamaJama configuration
	solve_return = solver.solve_nlp();

	return solve_return;
}
