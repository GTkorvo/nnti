// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLPModL1.cpp

#include "../include/MeritFuncNLPModL1.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"

namespace {
typedef LinAlgPack::value_type	value_type;
using LinAlgPack::VectorSlice;
using LinAlgPack::Vector;

// Compute the term sum( mu(j) * abs(c(j)), j = 1,...,m )
value_type local_constr_term( const Vector& mu, const VectorSlice& c
	, const char func_name[] );

}	// end namespace

namespace ConstrainedOptimizationPack {

MeritFuncNLPModL1::MeritFuncNLPModL1()
	: deriv_(0.0), mu_(0.0)
{}

// Overridden from MeritFuncNLP

value_type MeritFuncNLPModL1::value(value_type f, const VectorSlice& c) const
{
	using LinAlgPack::norm_1;
	return f + local_constr_term( mu_, c, "calc_deriv" );
}

value_type MeritFuncNLPModL1::deriv() const
{
	return deriv_;
}

void MeritFuncNLPModL1::print_merit_func(std::ostream& out
	, const std::string& L ) const
{
	out
		<< L << "*** Define a modified L1 merit funciton that uses different\n"
		<< L << "*** penalty parameters for each constriant.\n"
		<< L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
		<< L << "phi(f,c) = f + sum( mu(j) * abs(c(j)), j = 1,...,m )\n"
		<< L << "Dphi(x_k,d_k) = Gf_k' * d_k - sum( mu(j) * abs(c(j)), j = 1,...,m )\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPModL1::calc_deriv( const VectorSlice& Gf_k, const VectorSlice& c_k
	, const VectorSlice& d_k )
{
	using LinAlgPack::dot; using LinAlgPack::norm_1;
	return deriv_ = dot( Gf_k, d_k ) - local_constr_term( mu_, c_k, "calc_deriv" );
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPModL1::resize(size_type n)
{
	mu_.resize(n);
	mu_ = 0.0;
}

VectorSlice MeritFuncNLPModL1::mu()
{
	return mu_();
}

const VectorSlice MeritFuncNLPModL1::mu() const
{
	return mu_();
}

}	// end namespace ConstrainedOptimizationPack

namespace {

value_type local_constr_term( const Vector& mu, const VectorSlice& c
	, const char func_name[] )
{
	if( mu.size() != c.size() ) {
		std::ostringstream omsg;
		omsg
			<< "MeritFuncNLPModL1::" << func_name << "(...) : "
			<< "Error, the sizes mu.size() == " << mu.size()
			<< " != c.size() == " << c.size();
		throw ConstrainedOptimizationPack::MeritFuncNLP::InvalidInitialization(omsg.str());
	}
	value_type r = 0.0;
	Vector::const_iterator
		mu_itr = mu.begin();
	VectorSlice::const_iterator
		c_itr = c.begin();
	while( mu_itr != mu.end() )
		r += *mu_itr++ * ::fabs( *c_itr++ );
	return r;
}

}	// end namespace