// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIMultiVectorBase.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_HPP

#include <vector>

#include "TSFCoreMPIMultiVectorBaseDecl.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"
#include "RTOpCppToMPI.hpp"
#include "WorkspacePack.hpp"
#include "dynamic_cast_verbose.hpp"
#include "Teuchos_Time.hpp"
#ifdef RTOp_USE_MPI
#  include "Teuchos_RawMPITraits.hpp"
#endif

// Define to see some timing output!
//#define TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES

namespace TSFCore {

template<class Scalar>
MPIMultiVectorBase<Scalar>::MPIMultiVectorBase()
	:in_applyOp_(false)
	,globalDim_(0)
	,localOffset_(-1)
	,localSubDim_(0)
	,numCols_(0)
{}

// Overridden form OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
MPIMultiVectorBase<Scalar>::range() const
{
	return mpiSpace();
}

// Overridden from LinearOp

template<class Scalar>
void MPIMultiVectorBase<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	MultiVectorCols<Scalar>
		X(Teuchos::rcp(const_cast<Vector<Scalar>*>(&x),false)),
		Y(Teuchos::rcp(y,false));
	apply(M_trans,X,&Y,alpha,beta);
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();


#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	Teuchos::Time timerTotal("dummy",true);
	Teuchos::Time timer("dummy");
#endif

	//
	// This function performs one of two operations.
	//
	// The first operation (M_trans == NOTRANS) is:
	//
	//     Y = beta * Y + alpha * M * X
	//
	// where Y and M have compatible (distrubuted?) range vector
	// spaces and X is a locally replicated serial mulit-vector.  This
	// operation does not require any global communication.
	//
	// The second operation (M_trans == TRANS) is:
	//
	//     Y = beta * Y + alpha * M' * X
	//
	// where M and X have compatible (distrubuted?) range vector spaces
	// and Y is a locally replicated serial multi-vector.  This operation
	// requires a glocal reduction.
	//

	//
	// Get spaces and validate compatibility
	//

	const VectorSpace<Scalar>
		&Y_range = *Y->range(),
		&M_range = *this->range(),
		&X_range = *X.range();

	// Get the MPIVectorSpace
	const MPIVectorSpaceBase<Scalar> &mpiSpc = *this->mpiSpace();

	// Get the MPI communicator
	MPI_Comm mpiComm = mpiSpc.mpiComm();
#ifdef _DEBUG
//	std::cout << "MPIMultiVectorBase<Scalar>::apply(...): mpiComm = " << mpiComm << std::endl;
	TEST_FOR_EXCEPTION(
		( globalDim_ > localSubDim_ ) && mpiComm == MPI_COMM_NULL, std::logic_error
		,"MPIMultiVectorBase<Scalar>::apply(...MultiVector<Scalar>...): Error!"
		);
	// ToDo: Write a good general validation function that I can call tha will replace
	// all of these TEST_FOR_EXCEPTION(...) uses

	TEST_FOR_EXCEPTION(
		M_trans==NOTRANS && !mpiSpc.isCompatible(Y_range), Exceptions::IncompatibleVectorSpaces
		,"MPIMultiVectorBase<Scalar>::apply(...MultiVector<Scalar>...): Error!"
		);
	TEST_FOR_EXCEPTION(
		M_trans==TRANS && !mpiSpc.isCompatible(X_range), Exceptions::IncompatibleVectorSpaces
		,"MPIMultiVectorBase<Scalar>::apply(...MultiVector<Scalar>...): Error!"
		);
#endif

	//
	// Get explicit (local) views of Y, M and X
	//

#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.start();
#endif
		
	ExplicitMutableMultiVectorView<Scalar>
		Y_local(
			*Y
			,M_trans==NOTRANS ? Range1D(localOffset_+1,localOffset_+localSubDim_) : Range1D()
			,Range1D()
			);
	ExplicitMultiVectorView<Scalar>
		M_local(
			*this
			,Range1D(localOffset_+1,localOffset_+localSubDim_)
			,Range1D()
			);
	ExplicitMultiVectorView<Scalar>
		X_local(
			X
			,M_trans==NOTRANS ? Range1D() : Range1D(localOffset_+1,localOffset_+localSubDim_)
			,Range1D()
			);
#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nMPIMultiVectorBase<Scalar>::apply(...): Time for getting view = " << timer.totalElapsedTime() << " seconds\n";
#endif
#ifdef _DEBUG		
	TEST_FOR_EXCEPTION(
		M_trans==NOTRANS && ( M_local.numSubCols() != X_local.subDim() || X_local.numSubCols() != Y_local.numSubCols() )
		, Exceptions::IncompatibleVectorSpaces
		,"MPIMultiVectorBase<Scalar>::apply(...MultiVector<Scalar>...): Error!"
		);
	TEST_FOR_EXCEPTION(
		M_trans==TRANS && ( M_local.subDim() != X_local.subDim() || X_local.numSubCols() != Y_local.numSubCols() )
		, Exceptions::IncompatibleVectorSpaces
		,"MPIMultiVectorBase<Scalar>::apply(...MultiVector<Scalar>...): Error!"
		);
#endif

	//
	// If nonlocal (i.e. M_trans==TRANS) then create temorary storage
	// for:
	//
	//     Y_local_tmp = alpha * M(local) * X(local)                 : on nonroot processes
	//
	// or
	//
	//     Y_local_tmp = beta*Y_local + alpha * M(local) * X(local)  : on root process (localOffset_==0)
	// 
	// and set
	//
	//     localBeta = ( localOffset_ == 0 ? beta : 0.0 )
	//
	// Above, we choose localBeta such that we will only perform
	// Y_local = beta * Y_local + ...  on one process (the root
	// process where localOffset_==0x).  Then, when we add up Y_local
	// on all of the processors and we will get the correct result.
	//
	// If strictly local (i.e. M_trans == NOTRANS) then set:
	//
	//      Y_local_tmp = Y_local
	//      localBeta = beta
	//

#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.start();
#endif
		
	wsp::Workspace<Scalar> Y_local_tmp_store(wss,Y_local.subDim()*Y_local.numSubCols());
	RTOpPack::MutableSubMultiVectorT<Scalar> Y_local_tmp;
	Scalar localBeta;
	if( M_trans == TRANS && globalDim_ > localSubDim_ ) {
		// Nonlocal
		Y_local_tmp.initialize(
			0, Y_local.subDim()
			,0, Y_local.numSubCols()
			,&Y_local_tmp_store[0], Y_local.subDim() // leadingDim == subDim (columns are adjacent)
			);
		if( localOffset_ == 0 ) {
			// Root process: Must copy Y_local into Y_local_tmp
			for( int j = 0; j < Y_local.numSubCols(); ++j ) {
				Scalar *Y_local_j = Y_local.values() + Y_local.leadingDim()*j;
				std::copy( Y_local_j, Y_local_j + Y_local.subDim(), Y_local_tmp.values() + Y_local_tmp.leadingDim()*j );
			}
			localBeta = beta;
		}
		else {
			// Not the root process
			localBeta = 0.0;
		}
	}
	else {
		// Local
		Y_local_tmp = Y_local.smv(); // Shallow copy only!
		localBeta = beta;
	}

#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nMPIMultiVectorBase<Scalar>::apply(...): Time for setting up Y_local_tmp and localBeta = " << timer.totalElapsedTime() << " seconds\n";
#endif
		
	//
	// Perform the local multiplication:
	//
	//     Y(local) = localBeta * Y(local) + alpha * op(M(local)) * X(local)
	//
	// or in BLAS lingo:
	//
	//     C        = beta      * C        + alpha * op(A)        * op(B)
	//

#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
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
		,localBeta                                                   // BETA
		,Y_local_tmp.values()                                        // C
		,Y_local_tmp.leadingDim()                                    // LDC
		);
#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nMPIMultiVectorBase<Scalar>::apply(...): Time for GEMM = " << timer.totalElapsedTime() << " seconds\n";
#endif

#ifdef RTOp_USE_MPI
		
	//
	// Perform the global reduction of Y_local_tmp back into Y_local
	//
		
	if( M_trans==TRANS && globalDim_ > localSubDim_ ) {
		// Contiguous buffer for local send
		const Scalar *Y_local_buff = Y_local_tmp.values();
		// Contiguous buffer for final reduction
		wsp::Workspace<Scalar> Y_local_final_buff(wss,Y_local.subDim()*Y_local.numSubCols());
		// Perform the reduction
		MPI_Allreduce(
			Y_local_tmp.values()                     // sendbuff
			,&Y_local_final_buff[0]                  // recvbuff
			,Y_local_final_buff.size()               // count
			,Teuchos::RawMPITraits<Scalar>::type()   // datatype
			,MPI_SUM                                 // op
			,mpiComm                                 // comm
			);
		// Load Y_local_final_buff back into Y_local
		const Scalar *Y_local_final_buff_ptr = &Y_local_final_buff[0];
		for( int j = 0; j < Y_local.numSubCols(); ++j ) {
			Scalar *Y_local_ptr = Y_local.values() + Y_local.leadingDim()*j;
			for( int i = 0; i < Y_local.subDim(); ++i ) {
				(*Y_local_ptr++) = (*Y_local_final_buff_ptr++);
			}
		}
	}

#endif

	// When you get here the view Y_local will be committed back to Y
	// in the distructor to Y_local

#ifdef TSFCORE_MPI_MULTI_VECTOR_BASE_PRINT_TIMES
	timer.stop();
	std::cout << "\nMPIMultiVectorBase<Scalar>::apply(...): Total time = " << timerTotal.totalElapsedTime() << " seconds\n";
#endif

}

// Overridden from MultiVector

template<class Scalar>
void MPIMultiVectorBase<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &pri_op
	,const size_t                   num_multi_vecs
	,const MultiVector<Scalar>*     multi_vecs[]
	,const size_t                   num_targ_multi_vecs
	,MultiVector<Scalar>*           targ_multi_vecs[]
	,RTOp_ReductTarget              reduct_objs[]
	,const Index                    pri_first_ele_in
	,const Index                    pri_sub_dim_in
	,const Index                    pri_global_offset_in
	,const Index                    sec_first_ele_in
	,const Index                    sec_sub_dim_in
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	const Index numCols = this->domain()->dim();
	const MPIVectorSpaceBase<Scalar> &mpiSpc = *mpiSpace();
#ifdef _DEBUG
	// ToDo: Validate input!
	TEST_FOR_EXCEPTION(
		in_applyOp_, std::invalid_argument
		,"MPIMultiVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
		"clear sign that one of the methods getSubMultiVector(...), freeSubMultiVector(...) or commitSubMultiVector(...) "
		"was not implemented properly!"
		);
#endif
	// Flag that we are in applyOp()
	in_applyOp_ = true;
	// First see if this is a locally replicated vector in which case
	// we treat this as a local operation only.
	const bool locallyReplicated = (localSubDim_ == globalDim_);
	// Get the overlap in the current process with the input logical sub-vector
	// from (first_ele_in,sub_dim_in,global_offset_in)
	RTOp_index_type  overlap_first_local_ele  = 0;
	RTOp_index_type  overalap_local_sub_dim   = 0;
	RTOp_index_type  overlap_global_offset    = 0;
	RTOp_parallel_calc_overlap(
		globalDim_, localSubDim_, localOffset_, pri_first_ele_in, pri_sub_dim_in, pri_global_offset_in
		,&overlap_first_local_ele, &overalap_local_sub_dim, &overlap_global_offset
		);
	const Range1D
		local_rng = (
			overlap_first_local_ele!=0
			? Range1D( localOffset_ + overlap_first_local_ele, localOffset_ + overlap_first_local_ele + overalap_local_sub_dim - 1 )
			: Range1D::Invalid
			),
		col_rng(
			sec_first_ele_in
			,sec_sub_dim_in ? sec_first_ele_in + sec_sub_dim_in - 1 : numCols
			);
	// Create sub-vector views of all of the *participating* local data
	wsp::Workspace<RTOpPack::SubMultiVectorT<Scalar> > sub_multi_vecs(wss,num_multi_vecs);
	wsp::Workspace<RTOpPack::MutableSubMultiVectorT<Scalar> > targ_sub_multi_vecs(wss,num_targ_multi_vecs);
	if( overlap_first_local_ele != 0 ) {
		if(1){for(int k = 0; k < num_multi_vecs; ++k ) {
			multi_vecs[k]->getSubMultiVector( local_rng, col_rng, &sub_multi_vecs[k] );
			sub_multi_vecs[k].setGlobalOffset( overlap_global_offset );
		}}
		if(1){for(int k = 0; k < num_targ_multi_vecs; ++k ) {
			targ_multi_vecs[k]->getSubMultiVector( local_rng, col_rng, &targ_sub_multi_vecs[k] );
			targ_sub_multi_vecs[k].setGlobalOffset( overlap_global_offset );
		}}
	}
	// Apply the RTOp operator object (all processors must participate)
	RTOpPack::MPI_apply_op(
		locallyReplicated ? MPI_COMM_NULL : mpiSpc.mpiComm()                               // comm
		,pri_op                                                                            // op
		,-1                                                                                // root_rank (perform an all-reduce)
		,col_rng.size()                                                                    // num_cols
		,num_multi_vecs                                                                    // num_multi_vecs
		,num_multi_vecs && overlap_first_local_ele ? &sub_multi_vecs[0] : NULL             // sub_multi_vecs
		,num_targ_multi_vecs                                                               // num_targ_multi_vecs
		,num_targ_multi_vecs && overlap_first_local_ele ? &targ_sub_multi_vecs[0] : NULL   // targ_sub_multi_vecs
		,reduct_objs                                                                       // reduct_objs
		);
	// Free and commit the local data
	if(1){for(int k = 0; k < num_multi_vecs; ++k ) {
		sub_multi_vecs[k].setGlobalOffset(local_rng.lbound()-1);
		multi_vecs[k]->freeSubMultiVector( &sub_multi_vecs[k] );
	}}
	if(1){for(int k = 0; k < num_targ_multi_vecs; ++k ) {
		targ_sub_multi_vecs[k].setGlobalOffset(local_rng.lbound()-1);
		targ_multi_vecs[k]->commitSubMultiVector( &targ_sub_multi_vecs[k] );
	}}
	// Flag that we are leaving applyOp()
	in_applyOp_ = false;
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::getSubMultiVector(
	const Range1D                       &rowRng_in
	,const Range1D                      &colRng_in
	,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
	) const
{
	const Range1D rowRng = validateRowRange(rowRng_in);
	const Range1D colRng = validateColRange(colRng_in);
	if( rowRng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rowRng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		MultiVector<Scalar>::getSubMultiVector(rowRng_in,colRng_in,sub_mv);
		return;
	}
	// rng consists of all local data so get it!
	const Scalar *localValues = NULL;
	int leadingDim = 0;
	this->getLocalData(&localValues,&leadingDim);
	sub_mv->initialize(
		rowRng.lbound()-1                             // globalOffset
		,rowRng.size()                                // subDim
		,colRng.lbound()-1                            // colOffset
		,colRng.size()                                // numSubCols
		,localValues
		+(rowRng.lbound()-localOffset_-1)
		+(colRng.lbound()-1)*leadingDim               // values
		,leadingDim                                   // leadingDim
		);
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::freeSubMultiVector(
	RTOpPack::SubMultiVectorT<Scalar>* sub_mv
	) const
{
	if( sub_mv->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_mv->globalOffset()+sub_mv->subDim() ) {
		// Let the default implementation handle it!
		MultiVector<Scalar>::freeSubMultiVector(sub_mv);
		return;
	}
	freeLocalData( sub_mv->values() );
	sub_mv->set_uninitialized();
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::getSubMultiVector(
	const Range1D                                &rowRng_in
	,const Range1D                               &colRng_in
	,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
	)
{
	const Range1D rowRng = validateRowRange(rowRng_in);
	const Range1D colRng = validateColRange(colRng_in);
	if( rowRng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rowRng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		MultiVector<Scalar>::getSubMultiVector(rowRng_in,colRng_in,sub_mv);
		return;
	}
	// rng consists of all local data so get it!
	Scalar *localValues = NULL;
	int leadingDim = 0;
	this->getLocalData(&localValues,&leadingDim);
	sub_mv->initialize(
		rowRng.lbound()-1                             // globalOffset
		,rowRng.size()                                // subDim
		,colRng.lbound()-1                            // colOffset
		,colRng.size()                                // numSubCols
		,localValues
		+(rowRng.lbound()-localOffset_-1)
		+(colRng.lbound()-1)*leadingDim               // values
		,leadingDim                                   // leadingDim
		);
}

template<class Scalar>
void MPIMultiVectorBase<Scalar>::commitSubMultiVector(
	RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv
	)
{
	if( sub_mv->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_mv->globalOffset()+sub_mv->subDim() ) {
		// Let the default implementation handle it!
		MultiVector<Scalar>::commitSubMultiVector(sub_mv);
		return;
	}
	commitLocalData( sub_mv->values() );
	sub_mv->set_uninitialized();
}

// protected

template<class Scalar>
void MPIMultiVectorBase<Scalar>::updateMpiSpace()
{
	if(globalDim_ == 0) {
		const MPIVectorSpaceBase<Scalar> *mpiSpace = this->mpiSpace().get();
		if(mpiSpace) {
			globalDim_    = mpiSpace->dim();
			localOffset_  = mpiSpace->localOffset();
			localSubDim_  = mpiSpace->localSubDim();
			numCols_      = this->domain()->dim();
		}
		else {
			globalDim_    = 0;
			localOffset_  = -1;
			localSubDim_  = 0;
			numCols_      = 0;
		}
	}
}

template<class Scalar>
Range1D MPIMultiVectorBase<Scalar>::validateRowRange( const Range1D &rowRng_in ) const
{
	const Range1D rowRng = RangePack::full_range(rowRng_in,1,globalDim_);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		rowRng.lbound() < 1 || globalDim_ < rowRng.ubound(), std::invalid_argument
		,"MPIMultiVectorBase<Scalar>::validateRowRange(rowRng): Error, the range rowRng = ["
		<<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not "
		"in the range [1,"<<globalDim_<<"]!"
		);
#endif
	return rowRng;
}

template<class Scalar>
Range1D MPIMultiVectorBase<Scalar>::validateColRange( const Range1D &colRng_in ) const
{
	const Range1D colRng = RangePack::full_range(colRng_in,1,numCols_);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		colRng.lbound() < 1 || numCols_ < colRng.ubound(), std::invalid_argument
		,"MPIMultiVectorBase<Scalar>::validateColRange(colRng): Error, the range colRng = ["
		<<colRng.lbound()<<","<<colRng.ubound()<<"] is not "
		"in the range [1,"<<numCols_<<"]!"
		);
#endif
	return colRng;
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_HPP
