// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// /////////////////////////////////////////////////////////////////////////////
// TSFCore_apply_op_helper.hpp

#ifndef TSFCORE_APPLY_OP_HELPER_HPP
#define TSFCORE_APPLY_OP_HELPER_HPP

#include "TSFCore_apply_op_helper_decl.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreAssertOp.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"

template<class Scalar>
void TSFCore::apply_op_validate_input(
	const char                      func_name[]
	,const VectorSpace<Scalar>      &space
	,const RTOpPack::RTOpT<Scalar>  &op
	,const int                      num_vecs
	,const Vector<Scalar>*          vecs[]
	,const int                      num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOpPack::ReductTarget         *reduct_obj
	,const Index                    first_ele_in
	,const Index                    sub_dim_in
	,const Index                    global_offset_in
	)
{
	int k;
	const Index
		dim = space.dim();
	TEST_FOR_EXCEPTION(
		global_offset_in < 0, std::logic_error
		,func_name << " : Error!  global_offset_in = "
		<<global_offset_in<<" is not valid" );
	TEST_FOR_EXCEPTION(
		first_ele_in > dim, std::logic_error
		,func_name << " : Error!  first_ele_in = "
		<<first_ele_in<<" is not compatible with space.dim() = " << dim );
	TEST_FOR_EXCEPTION(
		sub_dim_in < 0 || (sub_dim_in > 0 && sub_dim_in > dim-(first_ele_in-1)), std::logic_error
		,func_name << " : Error!  first_ele_in = "
		<<first_ele_in<<" and sub_dim_in = "<<sub_dim_in
		<<" is not compatible with space.dim() = " << dim );
	for(k = 0; k < num_vecs; ++k)
		TSFCORE_ASSERT_VEC_SPACES(func_name,space,*vecs[k]->space());
	for(k = 0; k < num_targ_vecs; ++k)
		TSFCORE_ASSERT_VEC_SPACES(func_name,space,*targ_vecs[k]->space());
}

template<class Scalar>
void TSFCore::apply_op_validate_input(
	const char                      func_name[]
	,const VectorSpace<Scalar>      &domain
	,const VectorSpace<Scalar>      &range
	,const RTOpPack::RTOpT<Scalar>  &primary_op
	,const int                      num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const int                      num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOpPack::ReductTarget*        reduct_objs[]
	,const Index                    primary_first_ele_in
	,const Index                    primary_sub_dim_in
	,const Index                    primary_global_offset_in
	,const Index                    secondary_first_ele_in
	,const Index                    secondary_sub_dim_in
	)
{
	int k;
	// Validate primary range arguments
	const Index
		range_dim = range.dim();
	TEST_FOR_EXCEPTION(
		primary_global_offset_in < 0, std::logic_error
		,func_name << " : Error!  primary_global_offset_in = "
		<<primary_global_offset_in<<" is not valid" );
	TEST_FOR_EXCEPTION(
		primary_first_ele_in <= 0 || range_dim < primary_first_ele_in, std::logic_error
		,func_name << " : Error!  primary_first_ele_in = "
		<<primary_first_ele_in<<" is not compatible with range.dim() = " << range_dim );
	TEST_FOR_EXCEPTION(
		primary_sub_dim_in < 0 || (primary_sub_dim_in > 0 && primary_sub_dim_in > range_dim-(primary_first_ele_in-1))
		, std::logic_error
		,func_name << " : Error!  primary_first_ele_in = "
		<<primary_first_ele_in<<" and primary_sub_dim_in = "<<primary_sub_dim_in
		<<" are not compatible with range.dim() = " << range_dim );
	// Validate secondary domain arguments
	const Index
		domain_dim = domain.dim();
	TEST_FOR_EXCEPTION(
		secondary_first_ele_in <= 0 || domain_dim < secondary_first_ele_in, std::logic_error
		,func_name << " : Error!  secondary_first_ele_in = "
		<<secondary_first_ele_in<<" is not compatible with domain.dim() = " << domain_dim );
	TEST_FOR_EXCEPTION(
		secondary_sub_dim_in < 0 || (secondary_sub_dim_in > 0 && secondary_sub_dim_in > domain_dim-(secondary_first_ele_in-1))
		, std::logic_error
		,func_name << " : Error!  secondary_first_ele_in = "
		<<secondary_first_ele_in<<" and secondary_sub_dim_in = "<<secondary_sub_dim_in
		<<" are not compatible with domain.dim() = " << domain_dim );
	// Validate spaces
	for(k = 0; k < num_multi_vecs; ++k) {
		TSFCORE_ASSERT_VEC_SPACES(func_name,domain,*multi_vecs[k]->domain());
		TSFCORE_ASSERT_VEC_SPACES(func_name,range,*multi_vecs[k]->range());
	}
	for(k = 0; k < num_targ_multi_vecs; ++k) {
		TSFCORE_ASSERT_VEC_SPACES(func_name,domain,*targ_multi_vecs[k]->domain());
		TSFCORE_ASSERT_VEC_SPACES(func_name,range,*targ_multi_vecs[k]->range());
	}
}

template<class Scalar>
void TSFCore::apply_op_serial(
	const VectorSpace<Scalar>      &space
	,const RTOpPack::RTOpT<Scalar> &op
	,const int                     num_vecs
	,const Vector<Scalar>*         vecs[]
	,const int                     num_targ_vecs
	,Vector<Scalar>*               targ_vecs[]
	,RTOpPack::ReductTarget        *reduct_obj
	,const Index                   first_ele_in
	,const Index                   sub_dim_in
	,const Index                   global_offset_in
	)
{
 	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
	
	// Dimension of global sub-vector
	const Index
		full_dim       = space.dim(),
		global_sub_dim = sub_dim_in ? sub_dim_in : full_dim - (first_ele_in-1);
	const Range1D
		global_sub_rng = Range1D(first_ele_in,(first_ele_in-1)+global_sub_dim);

	//
	// Get explicit views of the vector elements
	//

	Workspace<RTOpPack::SubVectorT<Scalar> >         local_vecs(wss,num_vecs);
	Workspace<RTOpPack::MutableSubVectorT<Scalar> >  local_targ_vecs(wss,num_targ_vecs);
	int k;
	for(k = 0; k < num_vecs; ++k) {
		RTOpPack::SubVectorT<Scalar> &v = local_vecs[k];
		vecs[k]->getSubVector( global_sub_rng, &v );
		v.setGlobalOffset( global_offset_in );
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		RTOpPack::MutableSubVectorT<Scalar> &v = local_targ_vecs[k];
		targ_vecs[k]->getSubVector( global_sub_rng, &v );
		v.setGlobalOffset( global_offset_in );
	}

	//
	// Apply the reduction/transformation operator on all elements all at once!
	//

	op.apply_op(
		num_vecs,       num_vecs      ? &local_vecs[0]      : NULL
		,num_targ_vecs, num_targ_vecs ? &local_targ_vecs[0] : NULL
		,reduct_obj
		);

	//
	// Free (and commit) the explicit views of the vector elements which
	// should also inform the vectors that they have changed.
	//

	for(k = 0; k < num_vecs; ++k) {
		RTOpPack::SubVectorT<Scalar> &v = local_vecs[k];
		v.setGlobalOffset( global_sub_rng.lbound() - 1);
		vecs[k]->freeSubVector(&v);
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		RTOpPack::MutableSubVectorT<Scalar> &v = local_targ_vecs[k];
		v.setGlobalOffset( global_sub_rng.lbound() - 1);
		targ_vecs[k]->commitSubVector(&v);
	}

}

template<class Scalar>
void TSFCore::apply_op_serial(
	const VectorSpace<Scalar>       &domain
	,const VectorSpace<Scalar>      &range
	,const RTOpPack::RTOpT<Scalar>  &pri_op
	,const int                      num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const int                      num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOpPack::ReductTarget*        reduct_objs[]
	,const Index                    pri_first_ele_in
	,const Index                    pri_sub_dim_in
	,const Index                    pri_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	)
{

 	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
	
	// Primary range global sub-vector
	const Index
		range_dim          = range.dim(),
		pri_global_sub_dim = pri_sub_dim_in ? pri_sub_dim_in : range_dim - (pri_first_ele_in-1);
	const Range1D
		pri_global_sub_rng = Range1D(pri_first_ele_in,(pri_first_ele_in-1)+pri_global_sub_dim);
	// Secondary domain
	const Index
		domain_dim         = domain.dim(),
		sec_global_sub_dim = sec_sub_dim_in ? sec_sub_dim_in : domain_dim - (sec_first_ele_in-1);
	const Range1D
		sec_global_sub_rng = Range1D(sec_first_ele_in,(sec_first_ele_in-1)+sec_global_sub_dim);

	//
	// Get explicit views of the multi-vector elements
	//

	Workspace<RTOpPack::SubMultiVectorT<Scalar> >         local_multi_vecs(wss,num_multi_vecs);
	Workspace<RTOpPack::MutableSubMultiVectorT<Scalar> >  local_targ_multi_vecs(wss,num_targ_multi_vecs);
	int k;
	for(k = 0; k < num_multi_vecs; ++k) {
		RTOpPack::SubMultiVectorT<Scalar> &mv = local_multi_vecs[k];
		multi_vecs[k]->getSubMultiVector( pri_global_sub_rng, sec_global_sub_rng, &mv );
		mv.setGlobalOffset( pri_global_offset_in );
	}
	for(k = 0; k < num_targ_multi_vecs; ++k) {
		RTOpPack::MutableSubMultiVectorT<Scalar> &mv = local_targ_multi_vecs[k];
		targ_multi_vecs[k]->getSubMultiVector( pri_global_sub_rng, sec_global_sub_rng, &mv );
		mv.setGlobalOffset( pri_global_offset_in );
	}

	//
	// Apply the reduction/transformation operator one column at a time
	//

	Workspace<RTOpPack::SubVectorT<Scalar> >         local_vecs(wss,num_multi_vecs);
	Workspace<RTOpPack::MutableSubVectorT<Scalar> >  local_targ_vecs(wss,num_targ_multi_vecs);

	for(int j = 0; j < sec_global_sub_dim; ++j ) {
		for(k = 0; k < num_multi_vecs; ++k)       local_vecs[k]      = local_multi_vecs[k].col(j+1);
		for(k = 0; k < num_targ_multi_vecs; ++k)  local_targ_vecs[k] = local_targ_multi_vecs[k].col(j+1);
		pri_op.apply_op(
			num_multi_vecs,       num_multi_vecs      ? &local_vecs[0]      : NULL
			,num_targ_multi_vecs, num_targ_multi_vecs ? &local_targ_vecs[0] : NULL
			,reduct_objs ? reduct_objs[j] : NULL
			);
	}

	//
	// Free (and commit) the explicit views of the multi-vector elements
	// which should also inform the vectors that they have changed.
	//

	for(k = 0; k < num_multi_vecs; ++k) {
		RTOpPack::SubMultiVectorT<Scalar> &mv = local_multi_vecs[k];
		mv.setGlobalOffset( pri_global_sub_rng.lbound() - 1);
		multi_vecs[k]->freeSubMultiVector(&mv);
	}
	for(k = 0; k < num_targ_multi_vecs; ++k) {
		RTOpPack::MutableSubMultiVectorT<Scalar> &mv = local_targ_multi_vecs[k];
		mv.setGlobalOffset( pri_global_sub_rng.lbound() - 1);
		targ_multi_vecs[k]->commitSubMultiVector(&mv);
	}

}

#endif // TSFCORE_APPLY_OP_HELPER_HPP
