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
// TSFCoreMPIMultiVectorBaseDecl.hpp

#ifndef TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
#define TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP

#include "TSFCoreMultiVectorDecl.hpp"
#include "TSFCoreEuclideanLinearOpBaseDecl.hpp"
#include "Teuchos_BLAS.hpp"

namespace TSFCore {

///
template<class Scalar> class MPIVectorSpaceBase;

///
/** Base class for MPI-based SPMD multi-vectors.
 *
 * By inheriting from this base class, multi-vector implementations
 * allow their multi-vector objects to be seamlessly combined with
 * other MPI-based multi-vector objects (of different concrete types)
 * in <tt>applyOp()</tt> and <tt>apply()</tt>.  A big part of this
 * protocol is that every multi-vector object can expose an
 * <tt>MPIVectorSpaceBase</tt> object through the virtual function
 * <tt>mpiSpace()</tt>.
 *
 * This base class contains an implementation of <tt>applyOp()</tt>
 * that relies on implementations of the methods (<tt>const</tt>)
 * <tt>getSubMultiVector()</tt>, <tt>freeSubMultiVector()</tt>,
 * (non-<tt>const</tt>) <tt>getSubMultiVector()</tt> and
 * <tt>commitSubMultiVector()</tt> (which all have default
 * implementations in this subclass).  In essence, this implemenation
 * will only call the <tt>getSubMultiVector()</tt> methods using a
 * range of (global) indexes for elements that exist on the local
 * processor.  As long as the number of local elements on each
 * processor is fairly large, the virtual function call overhead will
 * be minimal and this will result in a near optimal implementation.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must overide only five functions:
 * <tt>mpiSpace()</tt>, <tt>getLocalData(const Scalar**,Index*)</tt>,
 * <tt>freeLocalData(const Scalar**,Index*)</tt>,
 * <tt>getLocalData(Scalar**,Index*)</tt>,
 * <tt>commitLocalData(Scalar**,Index*)</tt>.  Note that overidding
 * the <tt>mpiSpace()</tt> function requires implementing or using a
 * pre-implemented concrete <tt>MPIVectorSpace</tt> object.
 *
 * If the <tt>getSubMultiVector()</tt> methods are ever called with
 * index ranges outside of those of the local processor, then the
 * default implementations in <tt>MultiVector</tt> of all of the
 * methods (<tt>const</tt>) <tt>MultiVector::getSubMultiVector()</tt>,
 * <tt>MultiVector::freeSubMultiVector()</tt>, (non-<tt>const</tt>)
 * <tt>MultiVector::getSubMultiVector()</tt> and
 * <tt>MultiVector::commitSubMultiVector()</tt> are called in instead.
 * Alternatively, a subclass could provide more specialized
 * implemenations of these methods (for more efficient gather/scatter
 * operations) if desired but this should not be needed for most use
 * cases.
 *
 * It is interesting to note that in the above use case that the
 * explicit subvector access methods call on its default
 * implementation defined in <tt>MultiVector</tt> (which calls on
 * <tt>applyOp()</tt>) and the operator will be properly applied since
 * the version of <tt>applyOp()</tt> implemented in this class will
 * only request local vector data and hence there will only be two
 * levels of recussion for any call to an explicit subvector access
 * method.  This is a truly elegant result.
 *
 * As described in the documentation for <tt>MPIVectorSpaceBase</tt>,
 * it is possible that at runtime it may be discovered that the
 * mapping of vector data to processors does not fall under this
 * design in which case the method <tt>applyOp()</tt> should be
 * overridden to handle this which will of course remove the
 * possibility of interoperability with other MPI-based vector
 * objects.  As long as ghost data is not included this should never
 * be an issue.
 *
 * Note that multi-vector subclass derived from this base class must
 * only be directly used in SPMD mode for this to work properly.
 *
 * \ingroup TSFCore_adapters_MPI_support_grp
 */
template<class Scalar>
class MPIMultiVectorBase
	: virtual public MultiVector<Scalar>
	, virtual public EuclideanLinearOpBase<Scalar>
{
public:

	///
	using EuclideanLinearOpBase<Scalar>::apply;
	///
	using MultiVector<Scalar>::applyOp;

	/** @name  Constructors / initializers / accessors */
	//@{

	///
	MPIMultiVectorBase();

	//@}

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI-based vector space object for the range of <tt>*this</tt> multi-vectorr.
	 */
	virtual Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

	///
	/** Returns a non-<tt>const</tt> pointer to a fortran-style view of the local multi-vector data.
	 *
	 * @param  localValues [out] On output <tt>*localValues</tt> will point to 
	 *                     the first element in the first colum of the local multi-vector
	 *                     stored as a column-major dense Fortran-style matrix.
	 * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
	 *                     of the Fortran-style local multi-vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>localValues!=NULL</tt>
	 * <li> <tt>leadingDim!=NULL</tt>
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>*localValues!=NULL</tt>
	 * <li> <tt>*leadingDim!=0</tt>
	 * </ul>
	 *
	 * The function <tT>commitLocalData()</tt> must be called to
	 * commit changes to the data.
	 */
	virtual void getLocalData( Scalar **localValues, Index *leadingDim ) = 0;

	///
	/** Commit view of local data that was gotten from <tt>getLocalData()</tt>.
	 *
	 * @param  localValues [in/out] On input <tt>localValues</tt> must be the pointer set
	 *                     by <tt>getLocalData()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>localValues!=NULL</tt>
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>*this</tt> will be updated to the entires in <tt>*localValues</tt>.
	 * </ul>
	 */
	virtual void commitLocalData( Scalar *localValues ) = 0;

	///
	/** Returns a <tt>const</tt>  pointer to a fortran-style view of the local multi-vector data.
	 *
	 * @param  localValues [out] On output <tt>*localValues</tt> will point to 
	 *                     the first element in the first colum of the local multi-vector
	 *                     stored as a column-major dense Fortran-style matrix.
	 * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
	 *                     of the Fortran-style local multi-vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>localValues!=NULL</tt>
	 * <li> <tt>leadingDim!=NULL</tt>
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>*localValues!=NULL</tt>
	 * <li> <tt>*leadingDim!=0</tt>
	 * </ul>
	 */
	virtual void getLocalData( const Scalar **localValues, Index *leadingDim ) const = 0;

	///
	/** Free view of local data that was gotten from <tt>getLocalData()</tt>.
	 *
	 * @param  localValues [in/out] On input <tt>localValues</tt> must be the pointer set
	 *                     by <tt>getLocalData()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>localValues!=NULL</tt>
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>*this</tt> will be updated to the entires in <tt>*localValues</tt>.
	 * </ul>
	 */
	virtual void freeLocalData( const Scalar *localValues ) const = 0;

	//@}

	/** @name Overridden from EuclideanLinearOpBase */
	//@{

	/// Returns <tt>mpiSpace</tt>.
	Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;

	///
	/** Wraps the <tt>Vector</tt> objects in <tt>MultiVector</tt> objects then calls
	 * the <tt>MultiVector</tt> version of <tt>euclideanApply()</tt>
	 */
	void euclideanApply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;

	///
	/** Uses GEMM(...) and MPI_Allreduce(...) to implement.
	 *
	 * ToDo: Finish documentation!
	 */
	void euclideanApply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;

	//@}

	/** @name Overridden from LinearOp */
	//@{

	///
	/** Calls <tt>EuclideanLinearOp::apply()</tt> to disambiguate <tt>apply()</tt>
	 */
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;

	///
	/** Calls <tt>EuclideanLinearOp::apply()</tt> to disambiguate <tt>apply()</tt>
	 */
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;

	//@}

	/** @name Overridden from MultiVector */
	//@{
	///
	void applyOp(
		const RTOpPack::RTOpT<Scalar>   &primary_op
		,const int                      num_multi_vecs
		,const MultiVector<Scalar>*     multi_vecs[]
		,const int                      num_targ_multi_vecs
		,MultiVector<Scalar>*           targ_multi_vecs[]
		,RTOpPack::ReductTarget*        reduct_objs[]
		,const Index                    primary_first_ele
		,const Index                    primary_sub_dim
		,const Index                    primary_global_offset
		,const Index                    secondary_first_ele
		,const Index                    secondary_sub_dim
		) const;
	///
	void getSubMultiVector(
		const Range1D                       &rowRng
		,const Range1D                      &colRng
		,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
		) const;
	///
	void freeSubMultiVector( RTOpPack::SubMultiVectorT<Scalar>* sub_mv ) const;
	///
	void getSubMultiVector(
		const Range1D                                &rowRng
		,const Range1D                               &colRng
		,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
		);
	///
	void commitSubMultiVector( RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv );
	//@}

protected:

	///
	/** Subclasses should call whenever the structure of the VectorSpace changes.
	 *
	 * This function can be overridden by subclasses but this
	 * particualar function implementation must be called from within
	 * any override.
	 */
	virtual void updateMpiSpace();

	///
	/** Validate and resize the row range.
	 *
	 * This function throws an exception if the input range is invalid
	 */
	Range1D validateRowRange( const Range1D& rowRng ) const;

	///
	/** Validate and resize the column range.
	 *
	 * This function throws an exception if the input range is invalid
	 */
	Range1D validateColRange( const Range1D& rowCol ) const;
	
private:
	
	// ///////////////////////////////////////
	// Private data members
	
	mutable bool in_applyOp_;

	mutable Teuchos::BLAS<int,Scalar> blas_;

	// cached
	Index  globalDim_;
	Index  localOffset_;
	Index  localSubDim_;
	Index  numCols_;

	
}; // end class MPIMultiVectorBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_MULTI_VECTOR_BASE_DECL_HPP
