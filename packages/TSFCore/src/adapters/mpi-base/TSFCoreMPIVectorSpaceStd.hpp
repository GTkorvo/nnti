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
// TSFCoreMPIVectorSpaceStd.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_STD_HPP
#define TSFCORE_MPI_VECTOR_SPACE_STD_HPP

#include "TSFCoreMPIVectorSpaceStdDecl.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"
#include "TSFCoreMPIMultiVectorStd.hpp"
#include "TSFCoreMPIVectorStd.hpp"
//#include "TSFCoreVectorMultiVector.hpp"

namespace TSFCore {

template<class Scalar>
MPIVectorSpaceStd<Scalar>::MPIVectorSpaceStd()
  :mpiComm_(MPI_COMM_NULL),localSubDim_(0),globalDim_(0),numProc_(0),procRank_(0)
{
	updateState();
}

template<class Scalar>
MPIVectorSpaceStd<Scalar>::MPIVectorSpaceStd( MPI_Comm mpiComm, const Index localSubDim, const Index globalDim )
  :mpiComm_(MPI_COMM_NULL),localSubDim_(0),globalDim_(0),numProc_(0),procRank_(0)
{
  initialize(mpiComm,localSubDim,globalDim);
}

template<class Scalar>
void MPIVectorSpaceStd<Scalar>::initialize( MPI_Comm mpiComm, const Index localSubDim, const Index globalDim )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( !( localSubDim > 0 ) );
#endif
  mpiComm_     = mpiComm;
  localSubDim_ = localSubDim;
#ifdef RTOp_USE_MPI
  if( mpiComm != MPI_COMM_NULL ) {
		MPI_Comm_size( mpiComm_, &numProc_  );
		MPI_Comm_rank( mpiComm_, &procRank_ );
    if(globalDim < 0) {
			MPI_Allreduce(
				&localSubDim_                           // sendbuf
				,&globalDim_                            // recvbuf
				,1                                      // count
				,Teuchos::RawMPITraits<Index>::type()   // datatype
				,MPI_SUM                                // op
				,mpiComm                                // comm
				);
    }
    else {
#ifdef _DEBUG
      TEST_FOR_EXCEPT( globalDim == 0 );
#endif
      globalDim_ = globalDim; // Danger! we are taking the client's word for this!
    }
	}
	else {
#endif // RTOp_USE_MPI
		numProc_  = 1;
		procRank_ = 0;
    globalDim_ = localSubDim_;
#ifdef RTOp_USE_MPI
	}
#endif
	updateState();
}

template<class Scalar>
void MPIVectorSpaceStd<Scalar>::uninitialize( MPI_Comm *mpiComm, Index *localSubDim, Index *globalDim )
{
  if(mpiComm)     *mpiComm      = mpiComm_;
  if(localSubDim) *localSubDim  = localSubDim_;
  if(globalDim)   *globalDim    = globalDim_;

  mpiComm_      = MPI_COMM_NULL;
  localSubDim_  = 0;
  globalDim_    = 0;
}

// Overridden from VectorSpece

template<class Scalar>
Index MPIVectorSpaceStd<Scalar>::dim() const
{
	return globalDim_;
}

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
MPIVectorSpaceStd<Scalar>::createMember() const
{
  return Teuchos::rcp(
    new MPIVectorStd<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp( new Scalar[localSubDim_], Teuchos::DeallocArrayDelete<Scalar>(), true )
      ,1
      )
    );
  //return Teuchos::rcp(new VectorMultiVector<Scalar>(createMembers(1)));
}

template<class Scalar>
Teuchos::RefCountPtr< MultiVector<Scalar> >
MPIVectorSpaceStd<Scalar>::createMembers(int numMembers) const
{
  return Teuchos::rcp(
    new MPIMultiVectorStd<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(smallVecSpcFcty()->createVecSpc(numMembers),true)
      ,Teuchos::rcp( new Scalar[localSubDim_*numMembers], Teuchos::DeallocArrayDelete<Scalar>(), true )
      ,localSubDim_
      )
    );
}

template<class Scalar>
Teuchos::RefCountPtr<Vector<Scalar> >
MPIVectorSpaceStd<Scalar>::createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( localSubDim_ != raw_v.subDim() );
#endif
  return Teuchos::rcp(
    new MPIVectorStd<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp( raw_v.values(), false )
      ,raw_v.stride()
      )
    );
}

template<class Scalar>
Teuchos::RefCountPtr<const Vector<Scalar> >
MPIVectorSpaceStd<Scalar>::createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( localSubDim_ != raw_v.subDim() );
#endif
  return Teuchos::rcp(
    new MPIVectorStd<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp( const_cast<Scalar*>(raw_v.values()), false )
      ,raw_v.stride()
      )
    );
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVector<Scalar> >
MPIVectorSpaceStd<Scalar>::createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( localSubDim_ != raw_mv.subDim() );
#endif
  return Teuchos::rcp(
    new MPIMultiVectorStd<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(smallVecSpcFcty()->createVecSpc(raw_mv.numSubCols()),true)
      ,Teuchos::rcp( raw_mv.values(), false )
      ,raw_mv.leadingDim()
      )
    );
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVector<Scalar> >
MPIVectorSpaceStd<Scalar>::createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( localSubDim_ != raw_mv.subDim() );
#endif
  return Teuchos::rcp(
    new MPIMultiVectorStd<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(smallVecSpcFcty()->createVecSpc(raw_mv.numSubCols()),true)
      ,Teuchos::rcp( const_cast<Scalar*>(raw_mv.values()), false )
      ,raw_mv.leadingDim()
      )
    );
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
MPIVectorSpaceStd<Scalar>::clone() const
{
	return Teuchos::rcp(new MPIVectorSpaceStd<Scalar>(mpiComm_,localSubDim_,globalDim_));
}

// Overridden from MPIVectorSpaceBase

template<class Scalar>
MPI_Comm MPIVectorSpaceStd<Scalar>::mpiComm() const
{
	return mpiComm_;
}

template<class Scalar>
Index MPIVectorSpaceStd<Scalar>::localSubDim() const
{
	return localSubDim_;
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_SPACE_STD_HPP
