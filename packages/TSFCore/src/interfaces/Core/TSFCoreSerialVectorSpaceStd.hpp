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

// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorSpaceStd.hpp

// Define to use MultiVectorCols instead of SerialMultiVector
//#define TSFCORE_VECTOR_SPACE_USE_MULTI_VECTOR_COLS

#ifndef TSFCORE_SERIAL_VECTOR_SPACE_STD_HPP
#define TSFCORE_SERIAL_VECTOR_SPACE_STD_HPP

#include "TSFCoreSerialVectorSpaceStdDecl.hpp"
#include "TSFCoreSerialVectorSpaceBase.hpp"
#include "TSFCoreSerialVectorStd.hpp"
#ifndef TSFCORE_VECTOR_SPACE_USE_MULTI_VECTOR_COLS
#include "TSFCoreSerialMultiVectorStd.hpp"
#endif
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

template<class Scalar>
SerialVectorSpaceStd<Scalar>::SerialVectorSpaceStd( int dim )
{
	initialize(dim);
}

template<class Scalar>
void SerialVectorSpaceStd<Scalar>::initialize( int dim )
{
	dim_ = dim;
}

// Overridden from VectorSpace

template<class Scalar>
Index SerialVectorSpaceStd<Scalar>::dim() const
{
	return dim_;
}

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
SerialVectorSpaceStd<Scalar>::createMember() const
{
	return Teuchos::rcp(new SerialVectorStd<Scalar>(Teuchos::rcp(this,false)));
}

template<class Scalar>
Teuchos::RefCountPtr< MultiVector<Scalar> >
SerialVectorSpaceStd<Scalar>::createMembers(int numMembers) const
{
#ifndef TSFCORE_VECTOR_SPACE_USE_MULTI_VECTOR_COLS
	return Teuchos::rcp(
		new SerialMultiVectorStd<Scalar>(
			Teuchos::rcp(this,false)
			,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(this->smallVecSpcFcty()->createVecSpc(numMembers),true)
			)
		);
#else
	return VectorSpace<Scalar>::createMembers(numMembers);
#endif
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
SerialVectorSpaceStd<Scalar>::clone() const
{
	return Teuchos::rcp(new SerialVectorSpaceStd<Scalar>(*this));
}

} // end namespace TSFCore

#endif // TSFCORE_SERIAL_VECTOR_SPACE_STD_HPP
