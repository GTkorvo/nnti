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

#ifndef DECOMPOSITION_SYSTEM_STATE_STEP_BUILDER_STD_H
#define DECOMPOSITION_SYSTEM_STATE_STEP_BUILDER_STD_H

#include "MoochoPack_Types.hpp"
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
#include "AbstractLinAlgPack_BasisSystemPerm.hpp"
#endif
#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"

namespace OptionsFromStreamPack {
	class OptionsFromStream;
}

namespace MoochoPack {

///
/** Standard builder object for creating DecompositionSystem, EvalNewPoint Step and other objects
 * and setting up some of the state object.
 *
 * This class is designed to be used by <tt>NLPAlgoConfig</tt> subclasses based on SQP
 * and performs many different tasks that are common to all of these algorithms.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemStateStepBuilderStd
{
public:
	
	/** @name Enums for variaous options categories */
	//@{

	///
	enum ENullSpaceMatrixType {
		NULL_SPACE_MATRIX_AUTO, NULL_SPACE_MATRIX_EXPLICIT
		, NULL_SPACE_MATRIX_IMPLICIT };
	///
	enum ERangeSpaceMatrixType {
		RANGE_SPACE_MATRIX_AUTO, RANGE_SPACE_MATRIX_COORDINATE
		, RANGE_SPACE_MATRIX_ORTHOGONAL };

	//@}

	/** @name Struct for options values */
	//@{

	///
	struct SOptionValues {
		// Constructor (sets default values)
		SOptionValues();
		ENullSpaceMatrixType	null_space_matrix_type_;
		ERangeSpaceMatrixType	range_space_matrix_type_;
		int						max_dof_quasi_newton_dense_;    // If < 0, don't change default
	};

	//@}
	
	///
	typedef Teuchos::RefCountPtr<
		const OptionsFromStreamPack::OptionsFromStream>             options_ptr_t;

	///
	DecompositionSystemStateStepBuilderStd();

	///
	/** Set the options that will be used to configure the algorithmic objects.
	 *
	 *  @param  options
	 *               [in] If \c NULL then no options will be set.  If <tt>!=NULL</tt> then
	 *               this is the \c OptionsFromStream object that will be used to extract the
	 *               options to use for the algorithm.  The state of this object must
	 *               be maintained by the client until this object is no longer needed.
	 */
	void set_options( const options_ptr_t& options );
	///
	const options_ptr_t& get_options() const;
	///
	/** Process the %NLP and process the options passed in from <tt>set_options()</tt>.
	 * Postconditions:<ul>
	 * <li> <tt>this->current_option_values()</tt> returns the options that will be
	 *      used in all of the following method.
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	void process_nlp_and_options(
		std::ostream          *trase_out
		,NLP                  &nlp
		,NLPFirstOrder        **nlp_foi
		,NLPSecondOrder       **nlp_soi
		,NLPDirect            **nlp_fod
		,bool                 *tailored_approach
		);
	///
	/** Create the decomposition system object.
	 *
	 * ToDo: Finish documentation!
	 */
	void create_decomp_sys(
		std::ostream                                       *trase_out
		,NLP                                               &nlp
		,NLPFirstOrder                                     *nlp_foi
		,NLPSecondOrder                                    *nlp_soi
		,NLPDirect                                         *nlp_fod
		,bool                                              tailored_approach
		,Teuchos::RefCountPtr<DecompositionSystem>         *decomp_sys
		);
	///
	/** Add the common iteration quantities to the state object.
	 * 
	 * ToDo: Finish documentation!
	 */
	void add_iter_quantities(
		std::ostream                                           *trase_out
		,NLP                                                   &nlp
		,NLPFirstOrder                                         *nlp_foi
		,NLPSecondOrder                                        *nlp_soi
		,NLPDirect                                             *nlp_fod
		,bool                                                  tailored_approach
		,const Teuchos::RefCountPtr<DecompositionSystem>       &decomp_sys
		,const Teuchos::RefCountPtr<NLPAlgoState>              &state
		);

	///
	/** Create the EvalNewPoint step object and allocated objects.
	 *
	 * ToDo: Finish documentation!
	 */
	void create_eval_new_point(
		std::ostream                                                 *trase_out
		,NLP                                                         &nlp
		,NLPFirstOrder                                               *nlp_foi
		,NLPSecondOrder                                              *nlp_soi
		,NLPDirect                                                   *nlp_fod
		,bool                                                        tailored_approach
		,const Teuchos::RefCountPtr<DecompositionSystem>             &decomp_sys
		,Teuchos::RefCountPtr<IterationPack::AlgorithmStep>          *eval_new_point_step
		,Teuchos::RefCountPtr<CalcFiniteDiffProd>                    *calc_fd_prod
		,Teuchos::RefCountPtr<VariableBoundsTester>                  *bounds_tester
		,Teuchos::RefCountPtr<NewDecompositionSelection_Strategy>    *new_decomp_selection_strategy
		);

	///
	/** Return the current option values being used.
	 */
	SOptionValues& current_option_values();

private:
	
	// ///////////////////////////
	// Private data members

	/// Smart pointer to options
	options_ptr_t   options_;

	/// Options structs
	SOptionValues       uov_; // options set by user
	SOptionValues       cov_; // current option values actually used

#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
	Teuchos::RefCountPtr<BasisSystemPerm> basis_sys_perm_;
#endif

	// /////////////////////////
	// Private member functions

	/// Read in the options from a stream
	static void readin_options(
		const OptionsFromStreamPack::OptionsFromStream& options
		,SOptionValues *option_values, std::ostream* trase_out
		);

	/// Set the defaults for options not set by the user
	static void set_default_options(
		const SOptionValues& user_option_values
		,SOptionValues *current_option_values
		,std::ostream* trase_out
		);

}; // end class DecompositionSystemStateStepBuilderStd

// ///////////////////////////////
// Inline members

inline
DecompositionSystemStateStepBuilderStd::SOptionValues&
DecompositionSystemStateStepBuilderStd::current_option_values()
{
	return cov_;
}

}  // end namespace MoochoPack

#endif // DECOMPOSITION_SYSTEM_STATE_STEP_BUILDER_STD_H
