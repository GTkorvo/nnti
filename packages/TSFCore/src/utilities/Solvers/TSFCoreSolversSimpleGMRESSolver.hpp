// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversSimpleGMRESSolver.hpp

#ifndef TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_HPP
#define TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_HPP

#include "TSFCoreSolversSimpleGMRESSolverDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreSolversNorm.hpp"
#include "TSFCoreSolversConvergenceTester.hpp"
#include "TSFCoreTestingTools.hpp"
#include "check_nan_inf.h"
#include "ThrowException.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
SimpleGMRESSolver<Scalar>::SimpleGMRESSolver(
	const out_ptr_t   &out
	,bool             dump_all
	,int              default_max_iter
	,Scalar           default_tol
	)
	:out_(out)
	,dump_all_(dump_all)
	,default_max_iter_(default_max_iter)
	,default_tol_(default_tol)
{}

// Overridden from SolverState

template<class Scalar>
Index SimpleGMRESSolver<Scalar>::totalNumSystems() const
{
	assert(0);
	return 0;
}

template<class Scalar>
Index SimpleGMRESSolver<Scalar>::currNumSystems() const
{
	assert(0);
	return 0;
}

template<class Scalar>
int SimpleGMRESSolver<Scalar>::currIteration() const
{
	assert(0);
	return 0;
}

template<class Scalar>
void SimpleGMRESSolver<Scalar>::currActiveSystems( Index activeSystems[] ) const
{
	assert(0);
}

template<class Scalar>
void SimpleGMRESSolver<Scalar>::currEstRelResidualNorms( Scalar norms[] ) const
{
	assert(0);
}

// Overridden from IterativeLinearSolver

template<class Scalar>
bool SimpleGMRESSolver<Scalar>::adjointRequired() const
{
	return false;
}

//
// Here we will transform the problem from
//
//    1/a*op(M)*X[j] - Y[j] == 0
//
// to
//
//    Y_hat[j] - op(M_hat)*X_hat[j] == 0
//
// with preconditioner M_tilde_inv
//
//    where:
//        if M_tilde_right_inv == NULL
//            Y_hat       = a*Y
//            M_hat       = M
//            X_hat       = X
//            M_tilde_inv = M_tilde_left_inv
//        elif M_tilde_left_inv == NULL && M_tilde_right_inv != NULL
//            Y_hat       = a*Y
//            M_hat       = M*M_tilde_right_inv
//            X_hat       = inv(M_tilde_right_inv)*X
//            M_tilde_inv = NULL
//        elif M_tilde_left_inv != NULL && M_tilde_right_inv != NULL
//            Y_hat       = a*M_tilde_left_inv*Y
//            M_hat       = M_tilde_left_inv*M*M_tilde_right_inv
//            X_hat       = inv(M_tilde_right_inv)*X
//            M_tilde_inv = NULL
//        endif
//
// and then apply the CG method as exactly described
// in "Templates for the Solution of Linear Systems" except
// this is a simple multi-vector version.
//
template<class Scalar>
SolveReturn SimpleGMRESSolver<Scalar>::solve(
	const LinearOp<Scalar> &M, ETransp M_trans, const MultiVector<Scalar> &Y, MultiVector<Scalar> *X
	,const Scalar a, const int max_iter_in, ConvergenceTester<Scalar> *convTester
	,const LinearOp<Scalar> *M_tilde_left_inv, const ETransp M_tilde_left_inv_trans
	,const LinearOp<Scalar> *M_tilde_right_inv, ETransp M_tilde_right_inv_trans
	) const
{
	THROW_EXCEPTION(
		M_tilde_right_inv!=NULL || M_tilde_right_inv != NULL, std::logic_error
		,"Error, we can not handle preconditioners yet!"
		);
	namespace mmp = MemMngPack;
	const VectorSpace<Scalar> &opM_domain     = ( M_trans == NOTRANS                ? *M.domain() : *M.range()  );
	const VectorSpace<Scalar> &opM_range      = ( M_trans == NOTRANS                ? *M.range()  : *M.domain() );
	const ETransp      opM_notrans            = ( M_trans == NOTRANS                ? NOTRANS     : TRANS       );
	const ETransp      opM_trans              = ( M_trans == NOTRANS                ? TRANS       : NOTRANS     );
	const ETransp      opM_tilde_inv_notrans  = ( M_tilde_left_inv_trans == NOTRANS ? NOTRANS     : TRANS       );
	const ETransp      opM_tilde_inv_trans    = ( M_tilde_left_inv_trans == NOTRANS ? TRANS       : NOTRANS     );
	//
	const int totalNumSystems = Y.domain()->dim();
	//
	int j;
	//
	// Validate input
	//
#ifdef _DEBUG
	const char func_name[] = "SimpleGMRESSolver<Scalar>::solve(...)";
	THROW_EXCEPTION(X==NULL,std::invalid_argument,": Error!");
	if(M_tilde_left_inv) {
		const VectorSpace<Scalar> &opM_tilde_inv_domain = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->domain() : *M_tilde_left_inv->range()  );
		const VectorSpace<Scalar> &opM_tilde_inv_range  = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->range()  : *M_tilde_left_inv->domain() );
		const bool
			domain_compatible = opM_domain.isCompatible(opM_tilde_inv_range),
			range_compatible =  opM_range.isCompatible(opM_tilde_inv_domain);
		THROW_EXCEPTION(
			!(domain_compatible && range_compatible), TSFCore::Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_inv) do are not compatible!");
	}
	bool is_compatible = opM_domain.isCompatible(*X->range());
	THROW_EXCEPTION(
		!is_compatible, TSFCore::Exceptions::IncompatibleVectorSpaces
		,func_name<<": Error, the op(M).domain() not compatible with X->range()!");
	is_compatible = opM_range.isCompatible(*Y.range());
	THROW_EXCEPTION(
		!is_compatible, TSFCore::Exceptions::IncompatibleVectorSpaces
		,func_name<<": Error, the op(M).range() not compatible with Y.range()!");
	is_compatible = X->domain()->isCompatible(*Y.domain());
	THROW_EXCEPTION(
		!is_compatible, TSFCore::Exceptions::IncompatibleVectorSpaces
		,func_name<<": Error, the X->domain() not compatible with Y.domain()!");
#endif
	if(get_out().get()) {
		*get_out() << "\n*** Entering SimpleGMRESSolver<Scalar>::solve(...)\n" << std::setprecision(16);
		if(dump_all()) {
			*get_out() << "\nM =\n" << M;
			*get_out() << "\nM_trans = " << toString(M_trans) << std::endl;
			*get_out() << "\nY =\n" << Y;
			*get_out() << "\nX =\n" << *X;
			*get_out() << "\na = " << a << std::endl;
		}
	}
	//
	// Resolve default parameters
	//
	const int max_iter = ( max_iter_in == DEFAULT_MAX_ITER ? default_max_iter() : max_iter_in );
//	if (convTester) norm_ = convTester->norm(); else norm_ = mmp::rcp(new Solvers::Norm<Scalar>());
	//
	// Solve each linear system one at a time
	//
	bool all_solved = true;
	int max_iter_taken = 0;
	mmp::ref_count_ptr<Vector<Scalar> >
		y = Y.range()->createMember(),
		x = X->range()->createMember();
	for( int k = 1; k <= totalNumSystems; ++k ) {
		assign( &*y, *Y.col(k) );  Vt_S( &*y, a );  // y = a*Y.col(k)  (copy)
		x = X->col(k);                              // x = X.col(k)    (view)
		const SolveReturn single_solve_return
			= solver_.solve(
				M                       // Op
				,*y                     // b
				,&*x                    // curr_soln
				,M_trans                // Op_trans
				,max_iter               // max_iter_in
				,default_tol_           // tol_in
				);
		if( single_solve_return.solve_status == MAX_ITER_EXCEEDED ) all_solved = false;
		if( single_solve_return.num_iter > max_iter_taken ) max_iter_taken = single_solve_return.num_iter;
	}
	//
	// Return the solution
	//
	return SolveReturn( all_solved ? SOLVED_TO_TOL : MAX_ITER_EXCEEDED , max_iter_taken );
}

template<class Scalar>
MemMngPack::ref_count_ptr<const IterativeLinearSolver<Scalar> >
SimpleGMRESSolver<Scalar>::clone() const
{
	return MemMngPack::rcp( new SimpleGMRESSolver<Scalar>(*this) );
}


} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_HPP
