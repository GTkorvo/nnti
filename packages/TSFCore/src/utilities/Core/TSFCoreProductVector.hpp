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

// //////////////////////////////////////////////////////////////
// TSFCoreProductVector.hpp

#ifndef TSFCORE_PRODUCT_VECTOR_HPP
#define TSFCORE_PRODUCT_VECTOR_HPP

#include "TSFCoreProductVectorDecl.hpp"
#include "TSFCoreProductVectorSpace.hpp"
#include "WorkspacePack.hpp"

namespace TSFCore {

// Constructors/initializers/accessors

template <class Scalar>
ProductVector<Scalar>::ProductVector()
{
	uninitialize();
}

template <class Scalar>
ProductVector<Scalar>::ProductVector(
	const Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >  &productSpace
	,const Teuchos::RefCountPtr<Vector<Scalar> >                   vecs[]
	)
{
	initialize(productSpace,vecs);
}

template <class Scalar>
void ProductVector<Scalar>::initialize(
	const Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >  &productSpace
	,const Teuchos::RefCountPtr<Vector<Scalar> >                   vecs[]
	)
{
	// ToDo: Validate input!
	numBlocks_ = productSpace->numBlocks();
	productSpace_ = productSpace;
	vecs_.resize(numBlocks_);
	if(vecs) {
		std::copy( vecs, vecs + numBlocks_, &vecs_[0] );
	}
	else {
		for( int k = 0; k < numBlocks_; ++k )
			vecs_[k] = productSpace->getBlock(k)->createMember();
	}
}

template <class Scalar>
void ProductVector<Scalar>::uninitialize(
	Teuchos::RefCountPtr<const ProductVectorSpace<Scalar> >  *productSpace
	,Teuchos::RefCountPtr<Vector<Scalar> >                   vecs[]
	)
{
	if(productSpace) *productSpace = productSpace_;
	if(vecs) std::copy( &vecs_[0], &vecs_[0]+numBlocks_, vecs );
	productSpace_ = Teuchos::null;
	vecs_.resize(0);
	numBlocks_ = 0;
}

// Overridden from ProductVectorBase

template <class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
ProductVector<Scalar>::productSpace() const
{
	return productSpace_;
}

template <class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
ProductVector<Scalar>::getBlock(const int k)
{
	TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
	return vecs_[k];
}

template <class Scalar>
Teuchos::RefCountPtr<const Vector<Scalar> >
ProductVector<Scalar>::getBlock(const int k) const
{
	TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
	return vecs_[k];
}

// Overridden from Vector

template <class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
ProductVector<Scalar>::space() const
{
	return productSpace_;
}

template <class Scalar>
void ProductVector<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>    &op
	,const int                       num_vecs
	,const Vector<Scalar>*           vecs[]
	,const int                       num_targ_vecs
	,Vector<Scalar>*                 targ_vecs[]
	,RTOpPack::ReductTarget          *reduct_obj
	,const Index                     first_ele_in
	,const Index                     sub_dim_in
	,const Index                     global_offset_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	//
	const Index	n = productSpace_->dim();
	// Validate the compatibility of the vectors!
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!(1 <= first_ele_in && first_ele_in <= n), std::out_of_range
		,"ProductVector::applyOp(...): Error, "
		"first_ele_in = " << first_ele_in << " is not in range [1,"<<n<<"]" );
	TEST_FOR_EXCEPTION(
		sub_dim_in < 0, std::invalid_argument
		,"ProductVector::applyOp(...): Error, "
		"sub_dim_in = " << sub_dim_in << " is not acceptable" );
	TEST_FOR_EXCEPTION(
		global_offset_in < 0, std::invalid_argument
		,"ProductVector::applyOp(...): Error, "
		"global_offset_in = " << global_offset_in << " is not acceptable" );
	TEST_FOR_EXCEPTION(
		sub_dim_in > 0 && sub_dim_in - (first_ele_in - 1) > n, std::length_error
		,"ProductVector::applyOp(...): Error, "
		"global_offset_in = " << global_offset_in << ", sub_dim_in = " << sub_dim_in
		<< "first_ele_in = " << first_ele_in << " and n = " << n
		<< " are not compatible" );
	bool test_failed;
	for(int k = 0; k < num_vecs; ++k) {
		test_failed = !this->space()->isCompatible(*vecs[k]->space());
		TEST_FOR_EXCEPTION(
			test_failed, Exceptions::IncompatibleVectorSpaces
			,"ProductVector::applyOp(...): Error vecs["<<k<<"]->space() "
			<<"of type \'"<<typeid(*vecs[k]->space()).name()<<"\' is not compatible with this "
			<<"\'VectorSpaceBlocked\' vector space!"
			);
	}
	for(int k = 0; k < num_targ_vecs; ++k) {
		test_failed = !this->space()->isCompatible(*targ_vecs[k]->space());
		TEST_FOR_EXCEPTION(
			test_failed, Exceptions::IncompatibleVectorSpaces
			,"ProductVector::applyOp(...): Error targ_vecs["<<k<<"]->space() "
			<<"of type \'"<<typeid(*vecs[k]->space()).name()<<"\' is not compatible with this "
			<<"\'VectorSpaceBlocked\' vector space!"
			);
	}
#endif
	// Get the index of an incore-only input vector and input/output vector?
	const bool this_isInCore = productSpace_->isInCore();
	int incore_vec_k = -1, incore_targ_vec_k = -1;
	// Dynamic cast the pointers for the vector arguments
	wsp::Workspace<const ProductVector<Scalar>*>
		vecs_args(wss,num_vecs);
	for(int k = 0; k < num_vecs; ++k) {
		vecs_args[k] = dynamic_cast<const ProductVector<Scalar>*>(vecs[k]);
		if( vecs_args[k] == NULL ) {
			const bool isInCore_k = vecs[k]->space()->isInCore();
			if( this_isInCore && isInCore_k ) {
				incore_vec_k = k;
				break;
			}
			TEST_FOR_EXCEPTION(
				!this_isInCore || (this_isInCore && !isInCore_k), Exceptions::IncompatibleVectorSpaces
				,"ProductVector::applyOp(...): Error vecs["<<k<<"] "
				<<"of type \'"<<typeid(*vecs[k]).name()<<"\' does not support the "
				<<"\'ProductVector<Scalar>\' interface and is not an incore vector or this is not an incore vector!"
				);
		}
	}
	wsp::Workspace<ProductVector<Scalar>*>
		targ_vecs_args(wss,num_targ_vecs);
	for(int k = 0; k < num_targ_vecs; ++k) {
		targ_vecs_args[k] = dynamic_cast<ProductVector<Scalar>*>(targ_vecs[k]);
		if( targ_vecs_args[k] == NULL ) {
			const bool isInCore_k = targ_vecs[k]->space()->isInCore();
			if( this_isInCore && isInCore_k ) {
				incore_targ_vec_k = k;
				break;
			}
			TEST_FOR_EXCEPTION(
				!this_isInCore || (this_isInCore && !isInCore_k), Exceptions::IncompatibleVectorSpaces
				,"ProductVector::applyOp(...): Error targ_vecs["<<k<<"] "
				<<"of type \'"<<typeid(*targ_vecs[k]).name()<<"\' does not support the "
				<<"\'ProductVector<Scalar>\' interface and is not an incore vector or this is not an incore vector!"
				);
		}
	}
	// Let a incore-only vector handle this through explicit vector access
	// if all else fails?
	if( incore_vec_k >= 0 ) {
		vecs[incore_vec_k]->applyOp(
			op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
			,first_ele_in,sub_dim_in,global_offset_in
			);
		return;
	}
	else if ( incore_targ_vec_k >= 0 ) {
		targ_vecs[incore_targ_vec_k]->applyOp(
			op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
			,first_ele_in,sub_dim_in,global_offset_in
			);
		return;
	}
	// Perform the reduction on each vector segment at a time.
	const Index this_dim = n;
	const Index sub_dim  = ( sub_dim_in == 0
													 ? this_dim - (first_ele_in - 1)
													 : sub_dim_in );
	Index num_elements_remaining = sub_dim;
	const int  numBlocks = productSpace_->numBlocks();
	wsp::Workspace<const Vector<Scalar>*>
		sub_vecs(wss,num_vecs);
	wsp::Workspace<Vector<Scalar>*>
		sub_targ_vecs(wss,num_targ_vecs);
	Index g_off = -first_ele_in + 1;
	for(int k = 0; k < numBlocks; ++k) {
		const Index local_dim = productSpace_->getBlock(k)->dim();
		if( g_off < 0 && -g_off > local_dim ) {
			g_off += local_dim;
			continue;
		}
		const Index
			local_sub_dim = ( g_off >= 0
												? std::min( local_dim, num_elements_remaining )
												: std::min( local_dim + g_off, num_elements_remaining ) );
		if( local_sub_dim <= 0 )
			break;
		for( int i = 0; i < num_vecs; ++i )       // Fill constituent vectors for block k
			sub_vecs[i] = &*vecs_args[i]->vecs_[k];
		for( int j = 0; j < num_targ_vecs; ++j )  // Fill constituent target vectors for block k
			sub_targ_vecs[j] = &*targ_vecs_args[j]->vecs_[k];
		TSFCore::applyOp( 
			op
			,num_vecs,num_vecs?&sub_vecs[0]:NULL
			,num_targ_vecs,num_targ_vecs?&sub_targ_vecs[0]:NULL
			,reduct_obj
			,g_off < 0 ? -g_off + 1 : 1                                // first_ele
			,local_sub_dim                                             // sub_dim
			,g_off < 0 ? global_offset_in : global_offset_in + g_off   // global_offset
			);
		g_off += local_dim;
		num_elements_remaining -= local_sub_dim;
	}
	TEST_FOR_EXCEPT(!(num_elements_remaining==0));
}

template <class Scalar>
void ProductVector<Scalar>::getSubVector(
	const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec
	) const
{
	const Range1D
		rng = rng_in.full_range() ? Range1D( 1, productSpace_->dim()) : rng_in;
	int    kth_vector_space  = -1;
	Index  kth_global_offset = 0;
	productSpace_->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
	if( rng.lbound() + rng.size() <= kth_global_offset + 1 + vecs_[kth_vector_space]->space()->dim() ) {
		// This involves only one sub-vector so just return it.
		const_cast<const Vector<Scalar>*>(&*vecs_[kth_vector_space])->getSubVector(
			rng - kth_global_offset, sub_vec );
		sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
	}
	else {
		// Just let the default implementation handle this.  ToDo: In the futrue
		// we could manually construct an explicit sub-vector that spanned
		// two or more consitituent vectors but this would be a lot of work.
		// However, this would require the use of temporary memory but
		// so what.
		Vector<Scalar>::getSubVector(rng_in,sub_vec);
	}
}

template <class Scalar>
void ProductVector<Scalar>::freeSubVector(
	RTOpPack::SubVectorT<Scalar>* sub_vec
	) const
{
	if( sub_vec->values() == NULL ) return;
	int    kth_vector_space  = -1;
	Index  kth_global_offset = 0;
	productSpace_->getVecSpcPoss(sub_vec->globalOffset()+1,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
	if( sub_vec->globalOffset() + sub_vec->subDim() <= kth_global_offset +  vecs_[kth_vector_space]->space()->dim() ) {
		// This sub_vec was extracted from a single constituent vector
		sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
		vecs_[kth_vector_space]->freeSubVector(sub_vec);
	}
	else {
		// This sub_vec was created by the default implementation!
		Vector<Scalar>::freeSubVector(sub_vec);
	}
}

template <class Scalar>
void ProductVector<Scalar>::getSubVector(
	const Range1D& rng_in, RTOpPack::MutableSubVectorT<Scalar>* sub_vec
	)
{
	const Range1D
		rng = rng_in.full_range() ? Range1D( 1, productSpace_->dim()) : rng_in;
	int    kth_vector_space  = -1;
	Index  kth_global_offset = 0;
	productSpace_->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
	if( rng.lbound() + rng.size() <= kth_global_offset + 1 + vecs_[kth_vector_space]->space()->dim() ) {
		// This involves only one sub-vector so just return it.
		vecs_[kth_vector_space]->getSubVector(
			rng - kth_global_offset, sub_vec );
		sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
	}
	else {
		// Just let the default implementation handle this.  ToDo: In the futrue
		// we could manually construct an explicit sub-vector that spanned
		// two or more consitituent vectors but this would be a lot of work.
		// However, this would require the use of temporary memory but
		// so what.
		Vector<Scalar>::getSubVector(rng_in,sub_vec);
	}
}

template <class Scalar>
void ProductVector<Scalar>::commitSubVector(
	RTOpPack::MutableSubVectorT<Scalar>* sub_vec
	)
{
	if( sub_vec->values() == NULL ) return;
	int    kth_vector_space  = -1;
	Index  kth_global_offset = 0;
	productSpace_->getVecSpcPoss(sub_vec->globalOffset()+1,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
	if( sub_vec->globalOffset() + sub_vec->subDim() <= kth_global_offset +  vecs_[kth_vector_space]->space()->dim() ) {
		// This sub_vec was extracted from a single constituent vector
		sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
		vecs_[kth_vector_space]->commitSubVector(sub_vec);
	}
	else {
		// This sub_vec was created by the default implementation!
		Vector<Scalar>::commitSubVector(sub_vec);
	}
}

template <class Scalar>
void ProductVector<Scalar>::setSubVector(
	const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
	)
{
	int    kth_vector_space  = -1;
	Index  kth_global_offset = 0;
	productSpace_->getVecSpcPoss(sub_vec.globalOffset()+1,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
	if( sub_vec.globalOffset() + sub_vec.subDim() <= kth_global_offset + vecs_[kth_vector_space]->space()->dim() ) {
		// This sub-vector fits into a single constituent vector
		RTOpPack::SparseSubVectorT<Scalar> sub_vec_g = sub_vec;
		sub_vec_g.setGlobalOffset( sub_vec_g.globalOffset() - kth_global_offset );
		vecs_[kth_vector_space]->setSubVector(sub_vec_g);
	}
	else {
		// Let the default implementation take care of this.  ToDo: In the futrue
		// it would be possible to manualy set the relavent constituent
		// vectors with no temp memory allocations.
		Vector<Scalar>::setSubVector(sub_vec);
	}
}

} // namespace TSFCore

#endif // TSFCORE_PRODUCT_VECTOR_HPP
