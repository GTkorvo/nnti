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

// ///////////////////////////////////////////////////////////
// TSFCoreMultiVector.hpp

#ifndef TSFCORE_MULTI_VECTOR_HPP
#define TSFCORE_MULTI_VECTOR_HPP

#include "TSFCoreMultiVectorDecl.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreLinearOp.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreAssertOp.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "RTOp_TOp_assign_vectors.h"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

// Provide access to the columns as Vector objects

template<class Scalar>
Teuchos::RefCountPtr<const Vector<Scalar> >
MultiVector<Scalar>::col(Index j) const
{
	return const_cast<MultiVector*>(this)->col(j);
}

// Cloning

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MultiVector<Scalar>::clone_mv() const
{
	const VectorSpace<Scalar>
		&domain = *this->domain(),
		&range  = *this->range();
	Teuchos::RefCountPtr<MultiVector<Scalar> >
		copy = range.createMembers(domain.dim());
	assign( &*copy, *this );
	return copy;
}

// Sub-view methods

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
MultiVector<Scalar>::subView( const Range1D& colRng_in ) const
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore    *wss      = Teuchos::get_default_workspace_store().get();
	const VectorSpace<Scalar>  &domain   = *this->domain();
	const VectorSpace<Scalar>  &range    = *this->range();
	const Index                dimDomain = domain.dim();
	const Range1D              colRng    = RangePack::full_range(colRng_in,1,dimDomain);
	if( colRng.lbound() == 1 && static_cast<Index>(colRng.ubound()) == dimDomain )
		return Teuchos::rcp(this,false); // Takes all of the colunns!
	if( colRng.size() ) {
		// We have to create a view of a subset of the columns
		Workspace< Teuchos::RefCountPtr< Vector<Scalar> > >  col_vecs(wss,colRng.size());
		for( Index j = colRng.lbound(); j <= colRng.ubound(); ++j )
			col_vecs[j-colRng.lbound()] = Teuchos::rcp_const_cast<Vector<Scalar> >(this->col(j));
		return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(colRng.size()),&col_vecs[0]));
	}
	return Teuchos::null; // There was an empty set in colRng_in!
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MultiVector<Scalar>::subView( const Range1D& colRng_in )
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore    *wss      = Teuchos::get_default_workspace_store().get();
	const VectorSpace<Scalar>  &domain   = *this->domain();
	const VectorSpace<Scalar>  &range    = *this->range();
	const Index                dimDomain = domain.dim();
	const Range1D              colRng    = RangePack::full_range(colRng_in,1,dimDomain);
	if( colRng.lbound() == 1 && static_cast<Index>(colRng.ubound()) == dimDomain )
		return Teuchos::rcp(this,false); // Takes all of the colunns!
	if( colRng.size() ) {
		// We have to create a view of a subset of the columns
		Workspace< Teuchos::RefCountPtr< Vector<Scalar> > >  col_vecs(wss,colRng.size());
		for( Index j = colRng.lbound(); j <= colRng.ubound(); ++j )
			col_vecs[j-colRng.lbound()] = this->col(j);
		return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(colRng.size()),&col_vecs[0]));
	}
	return Teuchos::null; // There was an empty set in colRng_in!
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
MultiVector<Scalar>::subView( const int numCols, const int cols[] ) const
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore    *wss      = Teuchos::get_default_workspace_store().get();
	const VectorSpace<Scalar>  &domain   = *this->domain();
	const VectorSpace<Scalar>  &range    = *this->range();
#ifdef _DEBUG
	const Index                dimDomain = domain.dim();
	const char msg_err[] = "MultiVector<Scalar>::subView(numCols,cols[]): Error!";
 	TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
	// We have to create a view of a subset of the columns
	Workspace< Teuchos::RefCountPtr< Vector<Scalar> > > col_vecs(wss,numCols);
	for( int k = 0; k < numCols; ++k ) {
		const int col_k = cols[k];
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			col_k < 1 || dimDomain < col_k, std::invalid_argument
			,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [1,"<<dimDomain<<"]!"
			);
#endif
		col_vecs[k] = Teuchos::rcp_const_cast<Vector<Scalar> >(this->col(col_k));
	}
	return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(numCols),&col_vecs[0]));
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MultiVector<Scalar>::subView( const int numCols, const int cols[] )
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore    *wss      = Teuchos::get_default_workspace_store().get();
	const VectorSpace<Scalar>  &domain   = *this->domain();
	const VectorSpace<Scalar>  &range    = *this->range();
#ifdef _DEBUG
	const Index                dimDomain = domain.dim();
	const char msg_err[] = "MultiVector<Scalar>::subView(numCols,cols[]): Error!";
 	TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
	// We have to create a view of a subset of the columns
	Workspace< Teuchos::RefCountPtr< Vector<Scalar> > > col_vecs(wss,numCols);
	for( int k = 0; k < numCols; ++k ) {
		const int col_k = cols[k];
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			col_k < 1 || dimDomain < col_k, std::invalid_argument
			,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [1,"<<dimDomain<<"]!"
			);
#endif
		col_vecs[k] = this->col(col_k);
	}
	return Teuchos::rcp(new MultiVectorCols<Scalar>(this->range(),range.smallVecSpcFcty()->createVecSpc(numCols),&col_vecs[0]));
}

// Collective applyOp() methods

template<class Scalar>
void MultiVector<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &prim_op
	,const int                      num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const int                      num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOpPack::ReductTarget*        reduct_objs[]
	,const Index                    prim_first_ele_in
	,const Index                    prim_sub_dim_in
	,const Index                    prim_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	) const
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	// ToDo: Validate the input!

	const VectorSpace<Scalar>  &domain = *this->domain();
  const VectorSpace<Scalar>  &range  = *this->range();

	// Get the primary and secondary dimmensions.

	const Index	prim_dim = range.dim();
  const Index sec_dim = domain.dim();
	const Index prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim - prim_first_ele_in + 1 );
	const Index sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  -  sec_first_ele_in + 1  );
#ifdef _DEBUG
	const char err_msg[] = "MultiVector<Scalar>::applyOp(...): Error!";
	TEST_FOR_EXCEPTION( !(0 < prim_sub_dim && prim_sub_dim <= prim_dim), std::invalid_argument, err_msg );
	TEST_FOR_EXCEPTION( !(0 < sec_sub_dim  && sec_sub_dim  <= sec_dim),  std::invalid_argument, err_msg );
#endif

	//
	// Apply the reduction/transformation operator and transform the
	// target vectors and reduce each of the reduction objects.
	//

	Workspace< Teuchos::RefCountPtr<const Vector<Scalar> > >   vecs_s(wss,num_multi_vecs);
	Workspace<const Vector<Scalar>*>                           vecs(wss,num_multi_vecs,false);
	Workspace< Teuchos::RefCountPtr<Vector<Scalar> > >         targ_vecs_s(wss,num_targ_multi_vecs);
	Workspace<Vector<Scalar>*>                                 targ_vecs(wss,num_targ_multi_vecs,false);

	for(Index j = sec_first_ele_in; j <= sec_first_ele_in - 1 + sec_sub_dim; ++j) {
		// Fill the arrays of vector arguments
		{for(Index k = 0; k < static_cast<Index>(num_multi_vecs); ++k) {
			vecs_s[k] = multi_vecs[k]->col(j);
			vecs[k] = vecs_s[k].get();
		}}
		{for(Index k = 0; k < static_cast<Index>(num_targ_multi_vecs); ++k) {
			targ_vecs_s[k] = targ_multi_vecs[k]->col(j);
			targ_vecs[k] = targ_vecs_s[k].get();
		}}
		// Apply the reduction/transformation operator
		TSFCore::applyOp(
			prim_op
			,num_multi_vecs,      (num_multi_vecs      ? &vecs[0]      : NULL)
			,num_targ_multi_vecs, (num_targ_multi_vecs ? &targ_vecs[0] : NULL)
			,reduct_objs ? reduct_objs[j-1] : NULL
			,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
			);
	}
	// At this point all of the designated targ vectors in the target multi-vectors have
	// been transformed and all the reduction objects in reduct_obj[] have accumulated
	// the reductions.
}

template<class Scalar>
void MultiVector<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &prim_op
	,const RTOpPack::RTOpT<Scalar>  &sec_op
	,const int                      num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const int                      num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOpPack::ReductTarget         *reduct_obj
	,const Index                    prim_first_ele_in
	,const Index                    prim_sub_dim_in
	,const Index                    prim_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	) const
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	// ToDo: Validate the input!

	const VectorSpace<Scalar> &domain = *this->domain();
  const VectorSpace<Scalar> &range = *this->range();

	// Get the primary and secondary dimmensions.
  const Index prim_dim = range.dim();
  const Index sec_dim = domain.dim();
  const Index prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim - prim_first_ele_in + 1 );
	const Index sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  -  sec_first_ele_in + 1  );
#ifdef _DEBUG
	const char err_msg[] = "MultiVector<Scalar>::applyOp(...): Error!";
	TEST_FOR_EXCEPTION( !(0 < prim_sub_dim && prim_sub_dim <= prim_dim), std::invalid_argument, err_msg );
	TEST_FOR_EXCEPTION( !(0 < sec_sub_dim  && sec_sub_dim  <= sec_dim),  std::invalid_argument, err_msg );
#endif

	// Create a temporary buffer for the reduction objects of the primary reduction
	// so that we can call the companion version of this method.
	Workspace<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >
    rcp_reduct_objs(wss,reduct_obj!=NULL?sec_sub_dim:0);
	Workspace<RTOpPack::ReductTarget*>
    reduct_objs(wss,reduct_obj!=NULL?sec_sub_dim:0,false);
	if(reduct_obj) {
		for(Index k = 0; k < sec_sub_dim; ++k) {
      rcp_reduct_objs[k] = prim_op.reduct_obj_create();
      reduct_objs[k] = &*rcp_reduct_objs[k];
		}
	}
	
	// Call the campanion version that accepts an array of reduction objects
	this->applyOp(
		prim_op
		,num_multi_vecs,       multi_vecs
		,num_targ_multi_vecs,  targ_multi_vecs
		,reduct_obj ? &reduct_objs[0] : NULL
		,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
		,sec_first_ele_in,  sec_sub_dim_in
		);

	// Reduce all the reduction objects using the secondary reduction operator
	// into one reduction object and free the intermedate reduction objects.
	if(reduct_obj) {
		for(Index k = 0; k < sec_sub_dim; ++k) {
			sec_op.reduce_reduct_objs( *reduct_objs[k], reduct_obj );
		}
	}
}

// Explicit sub-multi-vector access

template<class Scalar>
void MultiVector<Scalar>::getSubMultiVector(
	const Range1D                       &rowRng_in
	,const Range1D                      &colRng_in
	,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
	) const
{
	const Index
		rangeDim  = this->range()->dim(),
		domainDim = this->domain()->dim();
	const Range1D
		rowRng = rowRng_in.full_range() ? Range1D(1,rangeDim)  : rowRng_in,
		colRng = colRng_in.full_range() ? Range1D(1,domainDim) : colRng_in;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		rowRng.ubound() > rangeDim, std::out_of_range
		,"MultiVector<Scalar>::getSubMultiVector(...): Error, rowRng = ["
		<<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not in the range = [1,"
		<<rangeDim<<"]!"
		);
	TEST_FOR_EXCEPTION(
		colRng.ubound() > domainDim, std::out_of_range
		,"MultiVector<Scalar>::getSubMultiVector(...): Error, colRng = ["
		<<colRng.lbound()<<","<<colRng.ubound()<<"] is not in the range = [1,"
		<<domainDim<<"]!"
		);
#endif
	// Allocate storage for the multi-vector (stored column-major)
	Scalar *values = new Scalar[ rowRng.size() * colRng.size() ];
	// Extract multi-vector values colum by column
	RTOpPack::SubVectorT<Scalar> sv; // uninitiaized by default
	for( int k = colRng.lbound(); k <= colRng.ubound(); ++k ) {
		Teuchos::RefCountPtr<const Vector<Scalar> > col_k = this->col(k);
		col_k->getSubVector( rowRng, &sv );
		for( int i = 0; i < rowRng.size(); ++i )
			values[ i + (k-1)*rowRng.size() ] = sv[i];
		col_k->freeSubVector( &sv );
	}
	// Initialize the multi-vector view object
	sub_mv->initialize(
		rowRng.lbound()-1            // globalOffset
		,rowRng.size()               // subDim
		,colRng.lbound()-1           // colOffset
		,colRng.size()               // numSubCols
		,values                      // values
		,rowRng.size()               // leadingDim
		);
}

template<class Scalar>
void MultiVector<Scalar>::freeSubMultiVector(
	RTOpPack::SubMultiVectorT<Scalar>* sub_mv
	) const
{
	// Here we just need to free the view and that is it!
	delete [] const_cast<Scalar*>(sub_mv->values());
	sub_mv->set_uninitialized();
}

template<class Scalar>
void MultiVector<Scalar>::getSubMultiVector(
	const Range1D                                &rowRng
	,const Range1D                               &colRng
	,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
	)
{
	// Use the non-const implementation since it does exactly the
	// correct thing in this case also!
	MultiVector<Scalar>::getSubMultiVector(
		rowRng, colRng
		,static_cast<RTOpPack::SubMultiVectorT<Scalar>*>(sub_mv)
		// This cast will work as long as MutableSubMultiVectorT
		// maintains no extra state over SubMultiVectorT (which it
		// currently does not) but this is something that I should
		// technically check for some how.
		);
}

template<class Scalar>
void MultiVector<Scalar>::commitSubMultiVector(
	RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		sub_mv==NULL, std::logic_error, "MultiVector<Scalar>::commitSubMultiVector(...): Error!"
		);
#endif
	// Set back the multi-vector values colum by column
	const Range1D rowRng(sub_mv->globalOffset()+1,sub_mv->globalOffset()+sub_mv->subDim());
	RTOpPack::MutableSubVectorT<Scalar> msv; // uninitiaized by default
	for( int k = sub_mv->colOffset()+1; k <= sub_mv->numSubCols(); ++k ) {
		Teuchos::RefCountPtr<Vector<Scalar> > col_k = this->col(k);
		col_k->getSubVector( rowRng, &msv );
		for( int i = 0; i < rowRng.size(); ++i )
			msv[i] = sub_mv->values()[ i + (k-1)*rowRng.size() ];
		col_k->commitSubVector( &msv );
	}
	// Free the memory
	delete [] const_cast<Scalar*>(sub_mv->values());
	// Zero out the view
	sub_mv->set_uninitialized();
}

// Overridden methods from LinearOp

template<class Scalar>
void MultiVector<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
#ifdef _DEBUG
	TSFCORE_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("MultiVector<Scalar>::apply()",*this,M_trans,x,y);
#endif
	const Index nc = this->domain()->dim();
	// y *= beta
	Vt_S(y,beta);
	// y += alpha*op(M)*x
	if(M_trans == NOTRANS) {
		//
		// y += alpha*M*x = alpha*M.col(1)*x(1) + ... + alpha*M.col(nc)*x(nc)
		//
		// Extract an explicit view of x
		RTOpPack::SubVectorT<Scalar> x_sub_vec;               
		x.getSubVector(Range1D(),&x_sub_vec);
		// Loop through and add the multiple of each column
		for(Index j = 1; j <= nc; ++j )
			Vp_StV( y, Scalar(alpha*x_sub_vec(j)), *this->col(j) );
		// Release the view of x
		x.freeSubVector(&x_sub_vec);
	}
	else {
		//
		//                   [ alpha*dot(M.col(1),x)  ]
		// y += alpha*M'*x = [ alpha*dot(M.col(2),x)  ]
		//                   [ ...                    ]
		//                   [ alpha*dot(M.col(nc),x) ]
		//
		// Extract an explicit view of y
		RTOpPack::MutableSubVectorT<Scalar> y_sub_vec;               
		y->getSubVector(Range1D(),&y_sub_vec);
		// Loop through and add to each element in y
		for(Index j = 1; j <= nc; ++j )
			y_sub_vec(j) += alpha*dot(*this->col(j),x);
		// Commit explicit view of y
		y->commitSubVector(&y_sub_vec);
	}
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
MultiVector<Scalar>::clone() const
{
	return this->clone_mv();
}

} // end namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_HPP
