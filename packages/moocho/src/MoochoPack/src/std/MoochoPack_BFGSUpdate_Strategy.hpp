///////////////////////////////////////////////////////
// BFGSUpdate_Strategy.h

#ifndef BFGS_UPDATE_STRATEGY_H
#define BFGS_UPDATE_STRATEGY_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "ReducedSpaceSQPPack/include/std/QuasiNewtonStats.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Strategy interface which contains the guts for a dampened BFGS update.
 *
 * This object can not change the flow of control or do anything fancy.  It just
 * performs the dampened update or skips it if the update is not 
 * sufficiently positive definite.
 *
 * See the printed documentation generated by #this->print_step(...)#.
 */
class BFGSUpdate_Strategy {
public:

	///
	/** <<std member comp>> members for whether to rescale
	 * the initial identity Hessian or not.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, rescale_init_identity )

	///
	/** <<std member comp>> members for whether to perform
	  * dampended quasi-newton updating or not.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, use_dampening )

	///
	enum ESecantTesting { SECANT_TEST_DEFAULT, SECANT_TEST_ALWAYS, SECANT_NO_TEST };

	///
	/** <<std member comp>> members how and if the secant property of the BFGS
	  * update is tested.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ESecantTesting, secant_testing )

	///
	/** <<std member comp>> members for the warning tolerance for
	  * the check of the secant property.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, secant_warning_tol )

	///
	/** <<std member comp>> members for the error tolerance for
	  * the check of the secant property.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, secant_error_tol )

	///
	BFGSUpdate_Strategy(
		bool               rescale_init_identity  = false
		,bool              use_dampening          = false
		,ESecantTesting    secant_testing         = SECANT_TEST_DEFAULT
		,value_type        secant_warning_tol     = 1e-6
		,value_type        secant_error_tol       = 1e-1
		);

	///
	/** Perform the BFGS update.
	 *
	 * The function performs a straight forward (possibly dampended) BFGS update
	 * that strictly satisfies the secant property #B * s_bfgs = y_bfgs#.
	 *
	 * See the printed documentation generated by #this->print_step(...)#.
	 *
	 * Preconditions\begin{itemize}
	 * \item #s_bfgs->size() == y_bfgs->size() == B->rows() == B->cols()# (throws ???)
	 * \end{itemize}
	 *
	 * @param s_bfgs        [in/w] Secant change vector on input.  May be modified as
	 *                      modified as workspace.
	 * @param y_bfgs        [in/w] Secant change vector on input.  May be modified as
 	 *                      modified as workspace.
	 * @param first_update  [in] If true then this is the first update after #B# was
	 *                      initialized to identity.  In this case #B# will be rescaled
	 *                      as #B = Iscale * I# before the update is performed if
	 *                      #this->rescale_init_identity() == true#.
	 * @param out           [out] Output stream journal data is written to.
	 * @param olevel        [in] Output level for printing to #out#.
	 * @param check_results [in] Helps determine if the secant property is tested or not
	 *                      after the update (see the printed documentation).
	 * @param B             [in/out] The matrix to be updated.  #B# must support the
	 *                      #MatrixSymSecantUpdateable# interface or an exception will be thrown.
	 * @param quasi_newton_stats
	 *                      [out] The quasi-newton statistics object that is updated to
	 *                      inform what happened durring the update.
	 */
	void perform_update(
		VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
		,std::ostream& out, EJournalOutputLevel olevel, bool check_results
		,MatrixWithOp *B, QuasiNewtonStats* quasi_newton_stats 
		);
	
	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

}; // end class BFGSUpdate_Strategy

}  // end namespace ReducedSpaceSQPPack

#endif BFGS_UPDATE_STRATEGY_H
