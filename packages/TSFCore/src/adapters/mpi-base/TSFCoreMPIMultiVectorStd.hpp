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
// TSFCoreMPIMultiVectorStd.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_STD_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_STD_HPP

// Define to make some verbose output
//#define TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT

#include "TSFCoreMPIMultiVectorStdDecl.hpp"
#include "TSFCoreMPIVectorStd.hpp"
//#include "TSFCoreVectorMultiVector.hpp"

namespace TSFCore {

// Constructors/initializers/accessors

template<class Scalar>
MPIMultiVectorStd<Scalar>::MPIMultiVectorStd()
  :leadingDim_(0)
{}

template<class Scalar>
MPIMultiVectorStd<Scalar>::MPIMultiVectorStd(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >          &mpiRangeSpace
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domainSpace
  ,const Teuchos::RefCountPtr<Scalar>                                    &localValues
  ,const Index                                                           leadingDim
  )
{
  initialize(mpiRangeSpace,domainSpace,localValues,leadingDim);
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::initialize(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >          &mpiRangeSpace
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domainSpace
  ,const Teuchos::RefCountPtr<Scalar>                                    &localValues
  ,const Index                                                           leadingDim
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(mpiRangeSpace.get()==NULL);
  TEST_FOR_EXCEPT(domainSpace.get()==NULL);
  TEST_FOR_EXCEPT(localValues.get()==NULL);
  TEST_FOR_EXCEPT(leadingDim < mpiRangeSpace->localSubDim());
#endif
  mpiRangeSpace_ = mpiRangeSpace;
  domainSpace_   = domainSpace;
  localValues_   = localValues;
  leadingDim_    = leadingDim;
  updateMpiSpace();
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >          *mpiRangeSpace
  ,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  *domainSpace
  ,Teuchos::RefCountPtr<Scalar>                                    *localValues
  ,Index                                                           *leadingDim
  )
{
  if(mpiRangeSpace) *mpiRangeSpace = mpiRangeSpace_;
  if(domainSpace)   *domainSpace   = domainSpace_;
  if(localValues)   *localValues   = localValues_;
  if(leadingDim)    *leadingDim    = leadingDim_;

  mpiRangeSpace_  = Teuchos::null;
  domainSpace_    = Teuchos::null;
  localValues_    = Teuchos::null;
  leadingDim_     = 0;

  updateMpiSpace();
}

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
MPIMultiVectorStd<Scalar>::domainScalarProdVecSpc() const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::domainScalarProdVecSpc() const called!\n";
#endif
  return domainSpace_;
}

// Overridden from MultiVector

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
MPIMultiVectorStd<Scalar>::col(Index j)
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::col() called!\n";
#endif
#ifdef _DEBUG
	TEST_FOR_EXCEPT( j < 1 || this->domain()->dim() < j );
#endif
  return Teuchos::rcp(
    new MPIVectorStd<Scalar>(
      mpiRangeSpace_
      ,Teuchos::rcp( (&*localValues_) + (j-1)*leadingDim_, false )
      ,1
      )
    );
  //return Teuchos::rcp(new VectorMultiVector<Scalar>(subView(Range1D(j,j))));
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MPIMultiVectorStd<Scalar>::subView( const Range1D& col_rng_in )
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::subView() called!\n";
#endif
	const Range1D colRng = validateColRange(col_rng_in);
  return Teuchos::rcp(
    new MPIMultiVectorStd<Scalar>(
      mpiRangeSpace_
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
				mpiRangeSpace_->smallVecSpcFcty()->createVecSpc(colRng.size())
				,true)
      ,Teuchos::rcp( (&*localValues_) + (colRng.lbound()-1)*leadingDim_, false )
      ,leadingDim_
      )
    );
}

// Overridden from MPIMultiVectorBase

template<class Scalar>
Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >
MPIMultiVectorStd<Scalar>::mpiSpace() const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::mpiSpace() const called!\n";
#endif
  return mpiRangeSpace_;
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::getLocalData( Scalar **localValues, Index *leadingDim )
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::getLocalData() called!\n";
#endif
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues==NULL );
	TEST_FOR_EXCEPT( leadingDim==NULL );
#endif
  *localValues = &*localValues_;
  *leadingDim  = leadingDim_;
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::commitLocalData( Scalar *localValues )
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::commitLocalData() called!\n";
#endif
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues!=&*localValues_ );
#endif
	// Nothing to commit!
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::getLocalData( const Scalar **localValues, Index *leadingDim ) const
{
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::getLocalData() called!\n";
#endif
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues==NULL );
	TEST_FOR_EXCEPT( leadingDim==NULL );
#endif
  *localValues = &*localValues_;
  *leadingDim  = leadingDim_;
}

template<class Scalar>
void MPIMultiVectorStd<Scalar>::freeLocalData( const Scalar *localValues ) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( localValues!=&*localValues_ );
#endif
#ifdef TSFCORE_MPI_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nMPIMultiVectorStd<Scalar>::commitLocalData() called!\n";
#endif
	// Nothing to free
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_STD_HPP
