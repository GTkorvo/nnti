// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrder.cpp
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

#include "NLPInterfacePack/src/abstract/interfaces/NLPFirstOrder.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixOp.hpp"
#include "Range1D.hpp"
#include "ThrowException.hpp"

namespace {
	const char name_Gc[] = "Gc";
	const char name_Gh[] = "Gh";
} // end namespace

namespace NLPInterfacePack {

// constructors

NLPFirstOrder::NLPFirstOrder()
	: Gc_(NULL), Gh_(NULL)
{}

void NLPFirstOrder::initialize(bool test_setup) {
	num_Gc_evals_ = 0;
	num_Gh_evals_ = 0;
	NLPObjGrad::initialize(test_setup);
}

// BasisSystem

const NLPFirstOrder::basis_sys_ptr_t
NLPFirstOrder::basis_sys() const
{
	return MemMngPack::null;
}

// <<std aggr>> members for Gc

void NLPFirstOrder::set_Gc(MatrixOp* Gc)
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	Gc_ = Gc;
}

MatrixOp* NLPFirstOrder::get_Gc()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(Gc_, false, name_Gc);
}

MatrixOp& NLPFirstOrder::Gc()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

const MatrixOp& NLPFirstOrder::Gc() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gc_, false, name_Gc);
}

// <<std aggr>> members for Gh

void NLPFirstOrder::set_Gh(MatrixOp* Gh)
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	Gh_ = Gh;
}

MatrixOp* NLPFirstOrder::get_Gh()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(Gh_, false, name_Gh);
}

MatrixOp& NLPFirstOrder::Gh()
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gh_, false, name_Gh);
}

const MatrixOp& NLPFirstOrder::Gh() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(Gh_, false, name_Gh);
}

// calculations

void NLPFirstOrder::calc_Gc(const Vector& x, bool newx) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(Gc_, "NLP::calc_Gc()", name_Gc);
	imp_calc_Gc(x,newx,first_order_info());
	num_Gc_evals_++;
}

void NLPFirstOrder::calc_Gh(const Vector& x, bool newx) const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(Gh_, "NLP::calc_Gh()", name_Gh);
	imp_calc_Gh(x,newx,first_order_info());
	num_Gh_evals_++;
}

size_type NLPFirstOrder::num_Gc_evals() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return num_Gc_evals_;
}

size_type NLPFirstOrder::num_Gh_evals() const
{
#ifdef _DEBUG
	THROW_EXCEPTION( this->mI() == 0, std::logic_error, "" );
#endif
	return num_Gh_evals_;
}

}	// end namespace NLPInterfacePack 
