// ////////////////////////////////////////////////////////////////////
// rSQPState.cpp

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <sstream>
#include <typeinfo>

#include "../include/rSQPState.h"

// local implementation functions.
namespace {

using ReducedSpaceSQPPack::IterQuantityAccess;
using ReducedSpaceSQPPack::AlgorithmState;
using ReducedSpaceSQPPack::rSQPState;

// Lookup and cast an iteration quantity.
template<class Type>
IterQuantityAccess<Type>& cast_to_type(AlgorithmState& state, AlgorithmState::iq_id_type iq_id
	, const std::string& iq_name, bool use_iq_id)
{
	GeneralIterationPack::IterQuantity
		&iter_quant = use_iq_id ? state.iter_quant(iq_id) : state.iter_quant(iq_name);

	IterQuantityAccess<Type> *iqa = dynamic_cast<IterQuantityAccess<Type>*>( &iter_quant );
	if(!iqa) {
		std::ostringstream omsg;
		omsg	<< "cast_to_type(...) : The concrete type for the rSQPState iteration quantity "
				<<	"iq_name = \"" << iq_name << "\" is "
				<< typeid(iter_quant).name()
				<< "and is not a subclass of type "
				<< typeid(IterQuantityAccess<Type>).name();
		throw rSQPState::InvalidType(omsg.str());
	}
	return *iqa;
}

template<class Type>
const IterQuantityAccess<Type>& cast_to_type(const AlgorithmState& state, AlgorithmState::iq_id_type iq_id
	, const std::string& iq_name, bool use_iq_id)
{
	return cast_to_type<Type>(const_cast<AlgorithmState&>(state),iq_id,iq_name,use_iq_id);
}

}	// end namespace


// rSQPState iteration quantities names

// NLP Problem Info 
const std::string ReducedSpaceSQPPack::x_name				= "x";
const std::string ReducedSpaceSQPPack::f_name				= "f";
const std::string ReducedSpaceSQPPack::Gf_name				= "Gf";
const std::string ReducedSpaceSQPPack::HL_name				= "HL";
const std::string ReducedSpaceSQPPack::c_name				= "c";
const std::string ReducedSpaceSQPPack::Gc_name				= "Gc";
// Constraint Gradient Null Space / Range Space Decomposition Info
const std::string ReducedSpaceSQPPack::Y_name				= "Y";
const std::string ReducedSpaceSQPPack::Z_name				= "Z";
const std::string ReducedSpaceSQPPack::U_name				= "U";
const std::string ReducedSpaceSQPPack::V_name				= "V";
// Search Direction Info
const std::string ReducedSpaceSQPPack::py_name				= "py";
const std::string ReducedSpaceSQPPack::Ypy_name				= "Ypy";
const std::string ReducedSpaceSQPPack::pz_name				= "pz";
const std::string ReducedSpaceSQPPack::Zpz_name				= "Zpz";
const std::string ReducedSpaceSQPPack::d_name				= "d";
// Reduced QP Subproblem Info
const std::string ReducedSpaceSQPPack::rGf_name				= "rGf";
const std::string ReducedSpaceSQPPack::rHL_name				= "rHL";
const std::string ReducedSpaceSQPPack::w_name				= "w";
const std::string ReducedSpaceSQPPack::zeta_name				= "zeta";
const std::string ReducedSpaceSQPPack::qp_grad_name			= "qp_grad";
const std::string ReducedSpaceSQPPack::eta_name				= "eta";
// Global Convergence Info
const std::string ReducedSpaceSQPPack::alpha_name			= "alpha";
const std::string ReducedSpaceSQPPack::mu_name				= "mu";
const std::string ReducedSpaceSQPPack::Delta_name			= "Delta";
const std::string ReducedSpaceSQPPack::phi_name				= "phi";
// KKT Info
const std::string ReducedSpaceSQPPack::opt_kkt_err_name		= "opt_kkt_err";
const std::string ReducedSpaceSQPPack::feas_kkt_err_name		= "feas_kkt_err";
const std::string ReducedSpaceSQPPack::GL_name				= "GL";
const std::string ReducedSpaceSQPPack::rGL_name				= "rGL";
const std::string ReducedSpaceSQPPack::lambda_name			= "lambda";
const std::string ReducedSpaceSQPPack::nu_name				= "nu";

void ReducedSpaceSQPPack::initialize_rsqp_names()
{
	// NLP Problem Info 
	const_cast<std::string&>(x_name)				= "x";
	const_cast<std::string&>(f_name)				= "f";
	const_cast<std::string&>(Gf_name)				= "Gf";
	const_cast<std::string&>(HL_name)				= "HL";
	const_cast<std::string&>(c_name)				= "c";
	const_cast<std::string&>(Gc_name)				= "Gc";
	// Constraint Gradient Null Space / Range Space Decomposition Info
	const_cast<std::string&>(Y_name)				= "Y";
	const_cast<std::string&>(Z_name)				= "Z";
	const_cast<std::string&>(U_name)				= "U";
	const_cast<std::string&>(V_name)				= "V";
	// Search Direction Info
	const_cast<std::string&>(py_name)				= "py";
	const_cast<std::string&>(Ypy_name)				= "Ypy";
	const_cast<std::string&>(pz_name)				= "pz";
	const_cast<std::string&>(Zpz_name)				= "Zpz";
	const_cast<std::string&>(d_name)				= "d";
	// Reduced QP Subproblem Info
	const_cast<std::string&>(rGf_name)				= "rGf";
	const_cast<std::string&>(rHL_name)				= "rHL";
	const_cast<std::string&>(w_name)				= "w";
	const_cast<std::string&>(zeta_name)				= "zeta";
	const_cast<std::string&>(qp_grad_name)			= "qp_grad";
	const_cast<std::string&>(eta_name)				= "eta";
	// Global Convergence Info
	const_cast<std::string&>(alpha_name)			= "alpha";
	const_cast<std::string&>(mu_name)				= "mu";
	const_cast<std::string&>(Delta_name)			= "Delta";
	const_cast<std::string&>(phi_name)				= "phi";
	// KKT Info
	const_cast<std::string&>(opt_kkt_err_name)		= "opt_kkt_err";
	const_cast<std::string&>(feas_kkt_err_name)		= "feas_kkt_err";
	const_cast<std::string&>(GL_name)				= "GL";
	const_cast<std::string&>(rGL_name)				= "rGL";
	const_cast<std::string&>(lambda_name)			= "lambda";
	const_cast<std::string&>(nu_name)				= "nu";
}

namespace ReducedSpaceSQPPack {

// /////////////////////////////////////////////////
// Static data members

// Note to maintainers.  It is important that these names match up with
// the enumerations that index into them.

const std::string * const rSQPState::iq_name_all_[rSQPState::num_quantites] = {

	// value_type names

	&f_name,				&zeta_name,					&eta_name,					&alpha_name,
	&mu_name,				&phi_name,					&opt_kkt_err_name,			&feas_kkt_err_name,

	// Vector names

	&x_name,				&Gf_name,					&c_name,					&py_name,
	&Ypy_name,				&pz_name,					&Zpz_name,					&d_name,
	&rGf_name,				&w_name,					&qp_grad_name,				&Delta_name,
	&GL_name,				&rGL_name,					&lambda_name,

	// SpVector names

	&nu_name,

	// MatrixWithOp names

	&HL_name,				&Gc_name,					&Y_name,
	&Z_name,				&U_name,					&V_name,					&rHL_name
};

const std::string * const  * const rSQPState::iq_name_[rSQPState::num_quantity_types] = {
	iq_name_all_,
	iq_name_all_ + num_value_type_quantities,
	iq_name_all_ + num_value_type_quantities + num_Vector_quantities,
	iq_name_all_ + num_value_type_quantities + num_Vector_quantities + num_SpVector_quantities
};

// ///////////////////////////////////////////////////////////////////
// Inline private functions

inline
rSQPState::IQA_value_type& rSQPState::iqa_value_type(E_IterQuantities_value_type iq) {
	return cast_to_type<value_type>( *this, iq_id_[VALUE_TYPE][iq], *iq_name_[VALUE_TYPE][iq]
		, iq_ids_initialized_ );
}

inline
const rSQPState::IQA_value_type& rSQPState::iqa_value_type(E_IterQuantities_value_type iq) const {
	return cast_to_type<value_type>( *this, iq_id_[VALUE_TYPE][iq], *iq_name_[VALUE_TYPE][iq]
		, iq_ids_initialized_ );
}


inline
rSQPState::IQA_Vector& rSQPState::iqa_Vector(E_IterQuantities_Vector iq) {
	return cast_to_type<VectorWithNorms>( *this, iq_id_[VECTOR][iq], *iq_name_[VECTOR][iq]
		, iq_ids_initialized_ );
}

inline
const rSQPState::IQA_Vector& rSQPState::iqa_Vector(E_IterQuantities_Vector iq) const {
	return cast_to_type<VectorWithNorms>( *this, iq_id_[VECTOR][iq], *iq_name_[VECTOR][iq]
		, iq_ids_initialized_ );
}

inline
rSQPState::IQA_SpVector& rSQPState::iqa_SpVector(E_IterQuantities_SpVector iq) {
	return cast_to_type<SpVector>( *this, iq_id_[SP_VECTOR][iq], *iq_name_[SP_VECTOR][iq]
		, iq_ids_initialized_ );
}

inline
const rSQPState::IQA_SpVector& rSQPState::iqa_SpVector(E_IterQuantities_SpVector iq) const {
	return cast_to_type<SpVector>( *this, iq_id_[SP_VECTOR][iq], *iq_name_[SP_VECTOR][iq]
		, iq_ids_initialized_ );
}

inline
rSQPState::IQA_MatrixWithOp& rSQPState::iqa_MatrixWithOp(E_IterQuantities_MatrixWithOp iq) {
	return cast_to_type<MatrixWithOp>( *this, iq_id_[MATRIX_WITH_OP][iq], *iq_name_[MATRIX_WITH_OP][iq]
		, iq_ids_initialized_ );
}

inline
const rSQPState::IQA_MatrixWithOp& rSQPState::iqa_MatrixWithOp(E_IterQuantities_MatrixWithOp iq) const {
	return cast_to_type<MatrixWithOp>( *this, iq_id_[MATRIX_WITH_OP][iq], *iq_name_[MATRIX_WITH_OP][iq]
		, iq_ids_initialized_ );
}

// /////////////////////////////////////////////////
// Member functions

rSQPState::rSQPState()
	: AlgorithmState(num_quantites)
		, num_basis_(0), iq_ids_initialized_(0)
{
	// Initialize array to iq_id arrays.
	typedef const AlgorithmState::iq_id_type ** iq_id_t;
	const_cast<iq_id_t>(iq_id_)[VALUE_TYPE]		= iq_id_all_;
	const_cast<iq_id_t>(iq_id_)[VECTOR]			= iq_id_all_ + num_value_type_quantities;
	const_cast<iq_id_t>(iq_id_)[SP_VECTOR]		= iq_id_all_ + num_value_type_quantities + num_Vector_quantities;
	const_cast<iq_id_t>(iq_id_)[MATRIX_WITH_OP]	= iq_id_all_ + num_value_type_quantities + num_Vector_quantities
													+ num_SpVector_quantities;
}

// Overridden from AlgorithmState

void rSQPState::erase_iter_quant(const std::string& iq_name) {
	iq_ids_initialized_ = false;
	AlgorithmState::erase_iter_quant(iq_name);
}


void rSQPState::dump_iter_quant(std::ostream& out) const {
	if( get_decomp_sys().get() ) {	
		out
			<< "\n*** DecompositionSystem ***\n"
			<< typeid(*get_decomp_sys().get()).name() << "\n";
		}
	AlgorithmState::dump_iter_quant(out);
}

// �std comp� members for decomp_sys.

void rSQPState::set_decomp_sys(const decomp_sys_ptr_t& decomp_sys) {
	decomp_sys_ = decomp_sys;
}

rSQPState::decomp_sys_ptr_t& rSQPState::get_decomp_sys() {
	return decomp_sys_;
}

const rSQPState::decomp_sys_ptr_t& rSQPState::get_decomp_sys() const {
	return decomp_sys_;
}

DecompositionSystem& rSQPState::decomp_sys() {
	return *decomp_sys_;
}

const DecompositionSystem& rSQPState::decomp_sys() const {
	return *decomp_sys_;
}


bool rSQPState::is_fully_setup() const {
	const AlgorithmState::iq_id_type *	iq_id_all_itr	= iq_id_all_;
	while(iq_id_all_itr != iq_id_all_ + num_quantites)
		if(*iq_id_all_itr++ == AlgorithmState::DOES_NOT_EXIST) return false;
	return true;
}

void rSQPState::initialize_fast_access() {
	const std::string * const *		iq_name_all_itr	= iq_name_all_;
	AlgorithmState::iq_id_type *	iq_id_all_itr	= iq_id_all_;
	while(iq_name_all_itr != iq_name_all_ + num_quantites)
		*iq_id_all_itr++ = get_iter_quant_id(**iq_name_all_itr++);
}
	
// Iteration Info

void rSQPState::set_num_basis(int num_basis) {
	num_basis_ = num_basis_;
}

int rSQPState::incr_num_basis() {
	return ++num_basis_;
}

int rSQPState::num_basis() const {
	return num_basis_;
}
	
// NLP Problem Info 

rSQPState::IQA_Vector& rSQPState::x() {
	return iqa_Vector(Q_x);
}

const rSQPState::IQA_Vector& rSQPState::x() const {
	return iqa_Vector(Q_x);
}


rSQPState::IQA_value_type& rSQPState::f() {
	return iqa_value_type(Q_f);
}

const rSQPState::IQA_value_type& rSQPState::f() const {
	return iqa_value_type(Q_f);
}


rSQPState::IQA_Vector& rSQPState::Gf() {
	return iqa_Vector(Q_Gf);
}

const rSQPState::IQA_Vector& rSQPState::Gf() const {
	return iqa_Vector(Q_Gf);
}


rSQPState::IQA_MatrixWithOp& rSQPState::HL() {
	return iqa_MatrixWithOp(Q_HL);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::HL() const {
	return iqa_MatrixWithOp(Q_HL);
}


rSQPState::IQA_Vector& rSQPState::c() {
	return iqa_Vector(Q_c);
}

const rSQPState::IQA_Vector& rSQPState::c() const {
	return iqa_Vector(Q_c);
}


rSQPState::IQA_MatrixWithOp& rSQPState::Gc() {
	return iqa_MatrixWithOp(Q_Gc);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::Gc() const {
	return iqa_MatrixWithOp(Q_Gc);
}
	
// Constraint Gradient Null Space / Range Space Decomposition Info.


rSQPState::IQA_MatrixWithOp& rSQPState::Y() {
	return iqa_MatrixWithOp(Q_Y);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::Y() const {
	return iqa_MatrixWithOp(Q_Y);
}


rSQPState::IQA_MatrixWithOp& rSQPState::Z() {
	return iqa_MatrixWithOp(Q_Z);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::Z() const {
	return iqa_MatrixWithOp(Q_Z);
}


rSQPState::IQA_MatrixWithOp& rSQPState::U() {
	return iqa_MatrixWithOp(Q_U);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::U() const {
	return iqa_MatrixWithOp(Q_U);
}


rSQPState::IQA_MatrixWithOp& rSQPState::V() {
	return iqa_MatrixWithOp(Q_V);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::V() const {
	return iqa_MatrixWithOp(Q_V);
}

// Search Direction Info

rSQPState::IQA_Vector& rSQPState::py() {
	return iqa_Vector(Q_py);
}

const rSQPState::IQA_Vector& rSQPState::py() const {
	return iqa_Vector(Q_py);
}


rSQPState::IQA_Vector& rSQPState::Ypy() {
	return iqa_Vector(Q_Ypy);
}

const rSQPState::IQA_Vector& rSQPState::Ypy() const {
	return iqa_Vector(Q_Ypy);
}

rSQPState::IQA_Vector& rSQPState::pz() {
	return iqa_Vector(Q_pz);
}

const rSQPState::IQA_Vector& rSQPState::pz() const {
	return iqa_Vector(Q_pz);
}


rSQPState::IQA_Vector& rSQPState::Zpz() {
	return iqa_Vector(Q_Zpz);
}

const rSQPState::IQA_Vector& rSQPState::Zpz() const {
	return iqa_Vector(Q_Zpz);
}

rSQPState::IQA_Vector& rSQPState::d() {
	return iqa_Vector(Q_d);
}

const rSQPState::IQA_Vector& rSQPState::d() const {
	return iqa_Vector(Q_d);
}


// QP Subproblem Info


rSQPState::IQA_Vector& rSQPState::rGf() {
	return iqa_Vector(Q_rGf);
}

const rSQPState::IQA_Vector& rSQPState::rGf() const {
	return iqa_Vector(Q_rGf);
}


rSQPState::IQA_MatrixWithOp& rSQPState::rHL() {
	return iqa_MatrixWithOp(Q_rHL);
}

const rSQPState::IQA_MatrixWithOp& rSQPState::rHL() const {
	return iqa_MatrixWithOp(Q_rHL);
}


rSQPState::IQA_Vector& rSQPState::w() {
	return iqa_Vector(Q_w);
}

const rSQPState::IQA_Vector& rSQPState::w() const {
	return iqa_Vector(Q_w);
}


rSQPState::IQA_value_type& rSQPState::zeta() {
	return iqa_value_type(Q_zeta);
}

const rSQPState::IQA_value_type& rSQPState::zeta() const {
	return iqa_value_type(Q_zeta);
}


rSQPState::IQA_Vector& rSQPState::qp_grad() {
	return iqa_Vector(Q_qp_grad);
}

const rSQPState::IQA_Vector& rSQPState::qp_grad() const {
	return iqa_Vector(Q_qp_grad);
}


rSQPState::IQA_value_type& rSQPState::eta() {
	return iqa_value_type(Q_eta);
}

const rSQPState::IQA_value_type& rSQPState::eta() const {
	return iqa_value_type(Q_eta);
}

// Global Convergence Info


rSQPState::IQA_value_type& rSQPState::alpha() {
	return iqa_value_type(Q_alpha);
}

const rSQPState::IQA_value_type& rSQPState::alpha() const {
	return iqa_value_type(Q_alpha);
}


rSQPState::IQA_value_type& rSQPState::mu() {
	return iqa_value_type(Q_mu);
}

const rSQPState::IQA_value_type& rSQPState::mu() const {
	return iqa_value_type(Q_mu);
}


rSQPState::IQA_Vector& rSQPState::Delta() {
	return iqa_Vector(Q_Delta);
}

const rSQPState::IQA_Vector& rSQPState::Delta() const {
	return iqa_Vector(Q_Delta);
}


rSQPState::IQA_value_type& rSQPState::phi() {
	return iqa_value_type(Q_phi);
}

const rSQPState::IQA_value_type& rSQPState::phi() const {
	return iqa_value_type(Q_phi);
}


// KKT Info

rSQPState::IQA_value_type& rSQPState::opt_kkt_err() {
	return iqa_value_type(Q_opt_kkt_err);
}

const rSQPState::IQA_value_type& rSQPState::opt_kkt_err() const {
	return iqa_value_type(Q_opt_kkt_err);
}

rSQPState::IQA_value_type& rSQPState::feas_kkt_err() {
	return iqa_value_type(Q_feas_kkt_err);
}

const rSQPState::IQA_value_type& rSQPState::feas_kkt_err() const {
	return iqa_value_type(Q_feas_kkt_err);
}

rSQPState::IQA_Vector& rSQPState::GL() {
	return iqa_Vector(Q_GL);
}

const rSQPState::IQA_Vector& rSQPState::GL() const {
	return iqa_Vector(Q_GL);
}

rSQPState::IQA_Vector& rSQPState::rGL() {
	return iqa_Vector(Q_rGL);
}

const rSQPState::IQA_Vector& rSQPState::rGL() const {
	return iqa_Vector(Q_rGL);
}

rSQPState::IQA_Vector& rSQPState::lambda() {
	return iqa_Vector(Q_lambda);
}

const rSQPState::IQA_Vector& rSQPState::lambda() const {
	return iqa_Vector(Q_lambda);
}


rSQPState::IQA_SpVector& rSQPState::nu() {
	return iqa_SpVector(Q_nu);
}

const rSQPState::IQA_SpVector& rSQPState::nu() const {
	return iqa_SpVector(Q_nu);
}

// Basis Pivot Info

IVector& rSQPState::var_perm_new() {
	return var_perm_new_;
}

const IVector& rSQPState::var_perm_new() const {
	return var_perm_new_;
}


IVector& rSQPState::con_perm_new() {
	return con_perm_new_;
}

const IVector& rSQPState::con_perm_new() const {
	return con_perm_new_;
}


IVector& rSQPState::var_perm_old() {
	return var_perm_old_;
}

const IVector& rSQPState::var_perm_old() const {
	return var_perm_old_;
}


IVector& rSQPState::con_perm_old() {
	return con_perm_old_;
}

const IVector& rSQPState::con_perm_old() const {
	return con_perm_old_;
}


Range1D& rSQPState::var_dep() {
	return var_dep_;
}

const Range1D& rSQPState::var_dep() const {
	return var_dep_;
}

Range1D& rSQPState::var_indep() {
	return var_indep_;
}

const Range1D& rSQPState::var_indep() const {
	return var_indep_;
}
	

Range1D& rSQPState::con_decomp() {
	return con_decomp_;
}

const Range1D& rSQPState::con_decomp() const {
	return con_decomp_;
}

Range1D& rSQPState::con_undecomp() {
	return con_undecomp_;
}

const Range1D& rSQPState::con_undecomp() const {
	return con_undecomp_;
}

}	// end namespace ReducedSpaceSQPPack
