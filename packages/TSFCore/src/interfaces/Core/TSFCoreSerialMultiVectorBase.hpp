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

#ifndef TSFCORE_SERIAL_MULTI_VECTOR_BASE_HPP
#define TSFCORE_SERIAL_MULTI_VECTOR_BASE_HPP

#include "TSFCoreSerialMultiVectorBaseDecl.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreEuclideanLinearOpBase.hpp"
#include "TSFCoreSerialVectorSpaceBase.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "TSFCore_apply_op_helper.hpp"
#include "RTOp_parallel_helpers.h"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Time.hpp"

// Define to see some timing output!
//#define TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES

namespace TSFCore {

template<class Scalar>
SerialMultiVectorBase<Scalar>::SerialMultiVectorBase()
	:in_applyOp_(false)
	,numRows_(0)
	,numCols_(0)
{}

// Overridden from EuclideanLinearOpBase

template<class Scalar>
void SerialMultiVectorBase<Scalar>::euclideanApply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	this->euclideanApply(M_trans,static_cast<const MultiVector<Scalar>&>(x),static_cast<MultiVector<Scalar>*>(y),alpha,beta);
}

template<class Scalar>
void SerialMultiVectorBase<Scalar>::euclideanApply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{

#ifdef TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES
	Teuchos::Time timerTotal("dummy",true);
	Teuchos::Time timer("dummy");
#endif

	//
	// This function performs one of two operations.
	//
	// The first operation (M_trans == NOTRANS) is:

	//     Y = beta * Y + alpha * M * X
	//
	// The second operation (M_trans == TRANS) is:
	//
	//     Y = beta * Y + alpha * M' * X
	//

#ifdef _DEBUG
	TSFCORE_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES("SerialMultiVectorBase<Scalar>::euclideanApply()",*this,M_trans,X,Y);
#endif

	//
	// Get explicit views of Y, M and X
	//

#ifdef TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.start();
#endif
	ExplicitMutableMultiVectorView<Scalar>  Y_local(*Y);
	ExplicitMultiVectorView<Scalar>         M_local(*this);
	ExplicitMultiVectorView<Scalar>         X_local(X);
#ifdef TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nSerialMultiVectorBase<Scalar>::apply(...): Time for getting view = " << timer.totalElapsedTime() << " seconds\n";
#endif
		
	//
	// Perform the multiplication:
	//
	//     Y(local) = localBeta * Y(local) + alpha * op(M(local)) * X(local)
	//
	// or in BLAS lingo:
	//
	//     C        = beta      * C        + alpha * op(A)        * op(B)
	//

#ifdef TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.start();
#endif
	blas_.GEMM(
		M_trans==NOTRANS ? Teuchos::NO_TRANS : Teuchos::TRANS        // TRANSA
		,Teuchos::NO_TRANS                                           // TRANSB
		,Y_local.subDim()                                            // M
		,Y_local.numSubCols()                                        // N
		,M_trans==NOTRANS ? M_local.numSubCols() : M_local.subDim()  // K
		,alpha                                                       // ALPHA
		,const_cast<Scalar*>(M_local.values())                       // A
		,M_local.leadingDim()                                        // LDA
		,const_cast<Scalar*>(X_local.values())                       // B
		,X_local.leadingDim()                                        // LDB
		,beta                                                        // BETA
		,Y_local.values()                                            // C
		,Y_local.leadingDim()                                        // LDC
		);
#ifdef TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nSerialMultiVectorBase<Scalar>::apply(...): Time for GEMM = " << timer.totalElapsedTime() << " seconds\n";
#endif

#ifdef TSFCORE_SERIAL_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nSerialMultiVectorBase<Scalar>::apply(...): Total time = " << timerTotal.totalElapsedTime() << " seconds\n";
#endif

}

// Overridden from LinearOp

template<class Scalar>
void SerialMultiVectorBase<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	this->apply(M_trans,static_cast<const MultiVector<Scalar>&>(x),static_cast<MultiVector<Scalar>*>(y),alpha,beta);
}

template<class Scalar>
void SerialMultiVectorBase<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	this->EuclideanLinearOpBase<Scalar>::apply(M_trans,X,Y,alpha,beta);
}

// Overridden from MultiVector

template<class Scalar>
void SerialMultiVectorBase<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &pri_op
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
	) const
{
#ifdef _DEBUG
	// ToDo: Validate input!
	TEST_FOR_EXCEPTION(
		in_applyOp_, std::invalid_argument
		,"SerialMultiVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
		"clear sign that one of the methods getSubMultiVector(...), freeSubMultiVector(...) or commitSubMultiVector(...) "
		"was not implemented properly!"
		);
	apply_op_validate_input(
		"SerialMultiVectorBase<Scalar>::applyOp(...)", *domain(), *range()
		,pri_op,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
		,reduct_objs,pri_first_ele_in,pri_sub_dim_in,pri_global_offset_in
		,sec_first_ele_in,sec_sub_dim_in
		);
#endif
	in_applyOp_ = true;
	apply_op_serial(
		*domain(),*range()
		,pri_op,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs
		,reduct_objs,pri_first_ele_in,pri_sub_dim_in,pri_global_offset_in
		,sec_first_ele_in,sec_sub_dim_in
		);
	in_applyOp_ = false;
}

template<class Scalar>
void SerialMultiVectorBase<Scalar>::getSubMultiVector(
	const Range1D                       &rowRng_in
	,const Range1D                      &colRng_in
	,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
	) const
{
	const Range1D rowRng = validateRowRange(rowRng_in);
	const Range1D colRng = validateColRange(colRng_in);
	const Scalar *localValues = NULL; int leadingDim = 0;
	this->getData(&localValues,&leadingDim);
	sub_mv->initialize(
		rowRng.lbound()-1                             // globalOffset
		,rowRng.size()                                // subDim
		,colRng.lbound()-1                            // colOffset
		,colRng.size()                                // numSubCols
		,localValues
		+(rowRng.lbound()-1)
		+(colRng.lbound()-1)*leadingDim               // values
		,leadingDim                                   // leadingDim
		);
}

template<class Scalar>
void SerialMultiVectorBase<Scalar>::freeSubMultiVector(
	RTOpPack::SubMultiVectorT<Scalar>* sub_mv
	) const
{
	freeData( sub_mv->values() );
	sub_mv->set_uninitialized();
}

template<class Scalar>
void SerialMultiVectorBase<Scalar>::getSubMultiVector(
	const Range1D                                &rowRng_in
	,const Range1D                               &colRng_in
	,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
	)
{
	const Range1D rowRng = validateRowRange(rowRng_in);
	const Range1D colRng = validateColRange(colRng_in);
	Scalar *localValues = NULL; int leadingDim = 0;
	this->getData(&localValues,&leadingDim);
	sub_mv->initialize(
		rowRng.lbound()-1                             // globalOffset
		,rowRng.size()                                // subDim
		,colRng.lbound()-1                            // colOffset
		,colRng.size()                                // numSubCols
		,localValues
		+(rowRng.lbound()-1)
		+(colRng.lbound()-1)*leadingDim               // values
		,leadingDim                                   // leadingDim
		);
}

template<class Scalar>
void SerialMultiVectorBase<Scalar>::commitSubMultiVector(
	RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv
	)
{
	commitData( sub_mv->values() );
	sub_mv->set_uninitialized();
}

// protected

template<class Scalar>
void SerialMultiVectorBase<Scalar>::updateSpace()
{
	if(numRows_ == 0) {
		const VectorSpace<Scalar> *range = this->range().get();
		if(range) {
			numRows_    = range->dim();
			numCols_    = this->domain()->dim();
		}
		else {
			numRows_    = 0;
			numCols_    = 0;
		}
	}
}

template<class Scalar>
Range1D SerialMultiVectorBase<Scalar>::validateRowRange( const Range1D &rowRng_in ) const
{
	const Range1D rowRng = RangePack::full_range(rowRng_in,1,numRows_);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		rowRng.lbound() < 1 || numRows_ < rowRng.ubound(), std::invalid_argument
		,"SerialMultiVectorBase<Scalar>::validateRowRange(rowRng): Error, the range rowRng = ["
		<<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not "
		"in the range [1,"<<numRows_<<"]!"
		);
#endif
	return rowRng;
}

template<class Scalar>
Range1D SerialMultiVectorBase<Scalar>::validateColRange( const Range1D &colRng_in ) const
{
	const Range1D colRng = RangePack::full_range(colRng_in,1,numCols_);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		colRng.lbound() < 1 || numCols_ < colRng.ubound(), std::invalid_argument
		,"SerialMultiVectorBase<Scalar>::validateColRange(colRng): Error, the range colRng = ["
		<<colRng.lbound()<<","<<colRng.ubound()<<"] is not "
		"in the range [1,"<<numCols_<<"]!"
		);
#endif
	return colRng;
}

} // end namespace TSFCore

#endif // TSFCORE_SERIAL_MULTI_VECTOR_BASE_HPP
