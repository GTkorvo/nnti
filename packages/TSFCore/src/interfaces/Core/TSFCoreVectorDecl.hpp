// /////////////////////////////////////////////////////////////
// TSFCoreVectorDecl.hpp

#ifndef TSFCORE_VECTOR_DECL_HPP
#define TSFCORE_VECTOR_DECL_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreTypes.hpp"
#include "RTOpCpp.hpp"
#include "RTOp_SparseSubVector.h"

namespace TSFCore {

///
/** Abstract interface for finite-dimensional dense vectors.
 *
 * This interface contains a mimimal set of operations.  The main
 * feature of this interface is the operation <tt>applyOp()</tt>.
 * Every standard (i.e. BLAS) and nearly every non-standard
 * element-wise operation that can be performed on a set of vectors
 * can be performed efficiently through reduction/transformation
 * operators.  More standard vector operations could be included in
 * this interface and allow for specialized implementations but in
 * general, assuming the sub-vectors are large enough, such
 * implementations would not be significantly faster than those
 * implemented through reduction/transformation operators.  There are
 * some operations however that can not always be efficiently
 * implemented with reduction/transforamtion operators and a few of
 * these important methods are included in this interface.  The
 * <tt>applyOp()</tt> method allows to client to specify a sub-set of
 * the vector elements to include in reduction/transformation
 * operation.  This greatly increases the generality of this vector
 * interface as vector objects can be used as sub objects in larger
 * composite vectors and sub-views of a vector can be created.
 *
 * <b>Heads Up !!!!!!</b>
 * There already exists reduction/transformation operator-based
 * implementatations of several standard vector operations and
 * some convenience functions that wrap these operators and call
 * <tt>applyOp()</tt> can be found \ref TSFCore_VectorStdOps_grp "here".
 *
 * This interface also allows a client to extract a sub-set of
 * elements in an explicit form non-mutable <tt>RTOpPack::SubVectorT<Scalar></tt>
 * non-mutable <tt>RTOpPack::MutableSubVectorT<Scalar></tt> objects using the
 * methods <tt>getSubVector()</tt>.  In general, this is very bad
 * thing to do and should be avoided at all costs.  However, there are
 * some situations where this is needed and therefore it is supported.
 * The default implementation of these methods use
 * reduction/transformation operators with <tt>applyOp()</tt> in order
 * to extract and set the needed elements.
 *
 * In addition to being able to extract an explicit non-mutable view
 * of some (small?) sub-set of elements, this interface allows a
 * client to either extract a explicit mutable sub-views using
 * <tt>getSubVector()</tt> or to set sub-vectors using
 * <tt>setSubVector()</tt>.
 *
 * <b>Notes for subclass developers</b>
 *
 * In order to create a concreate subclass of this interface, only two
 * methods must be overridden: <tt>space()</tt> and
 * <tt>applyOp()</tt>.  Overridding the <tt>space()</tt> method
 * requires defining a concreate <tt>VectorSpace</tt> class (which has
 * only three pure virtual methods).
 *
 * A subclass should only bother overidding the explicit sub-vector
 * access methods <tt>getSubVector()</tt>, <tt>freeSubVector()</tt>,
 * <tt>commitSubVector()</tt> and <tt>setSubVector()</tt> if profiling
 * shows that these methods are taking too much time.  In most use
 * cases, these default implementations will not impact the overall
 * runtime for a numerical algorithm in any significant way.
 */
template<class Scalar>
class Vector {
public:

	///
	virtual ~Vector() {}

	/** @name Pure virtual methods (must be overridden by subclass) */
	//@{

	///
	/** Returns a smart pointer to the vector space that this vector belongs to.
	 *
	 * A return value of <tt>space().get()==NULL</tt> is a flag that
	 * <tt>*this</tt> is not fully initialized.
	 */
	virtual MemMngPack::ref_count_ptr< const VectorSpace<Scalar> > space() const = 0;

	///
	/** Apply a reduction/transformation,operation over a set of vectors:
	 * <tt>op(op(v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 *
	 * The vector <tt>this</tt> that this method is called on is
	 * assumed to be one of the vectors in
	 * <tt>v[0]...v[nv-1],z[0]...z[nz-1]</tt>.  This method is not
	 * called directly by the client but is instead
	 * <tt>TSFCore::applyOp()</tt>.
	 *
	 * See the documentation for the method <tt>TSFCore::applyOp()</tt>
	 * for a description of what this method does.
	 */
	virtual void applyOp(
		const RTOpPack::RTOpT<Scalar>    &op
		,const size_t                    num_vecs
		,const Vector<Scalar>*           vecs[]
		,const size_t                    num_targ_vecs
		,Vector<Scalar>*                 targ_vecs[]
		,RTOp_ReductTarget               reduct_obj
		,const Index                     first_ele
		,const Index                     sub_dim
		,const Index                     global_offset
		) const = 0;

	//@}

	/** @name Explicit sub-vector access */
	//@{

	///
	/** Get a non-mutable explicit view of a sub-vector.
	 *
	 * @param  rng      [in] The range of the elements to extract the sub-vector view.
	 * @param  sub_vec  [in/out] View of the sub-vector.  Prior to the
	 *                  first call <tt>RTOp_sub_vector_null(sub_vec)</tt> must
	 *                  have been called for the correct behavior.  Technically
	 *                  <tt>*sub_vec</tt> owns the memory but this memory can be freed
	 *                  only by calling <tt>this->freeSubVector(sub_vec)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() <= this->dim()) == true</tt>
	 *      (<tt>throw std::out_of_range</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*sub_vec</tt> contains an explicit non-mutable view to the elements
	 *      in the range <tt>full_range(rng,1,this->dim())</tt>
	 * </ul>
	 *
	 * This is only a transient view of a sub-vector that is to be
	 * immediately used and then released with a call to
	 * <tt>freeSubVector()</tt>.
	 *
	 * Note that calling this operation might require some dynamic
	 * memory allocations and temporary memory.  Therefore, it is
	 * critical that <tt>this->freeSubSector(sub_vec)</tt> is called
	 * to clean up memory and avoid memory leaks after the sub-vector
	 * is used.
	 *
	 * If <tt>this->getSubVector(...,sub_vec)</tt> was previously
	 * called on <tt>sub_vec</tt> then it may be possible to reuse
	 * this memory if it is sufficiently sized.  The user is
	 * encouraged to make multiple calls to
	 * <tt>this->getSubVector(...,sub_vec)</tt> before
	 * <tt>this->freeSubSector(sub_vec)</tt> to finally clean up all
	 * of the memory.  Of course the same <tt>sub_vec</tt> object must
	 * be passed to the same vector object for this to work correctly.
	 *
	 * This method has a default implementation based on a vector
	 * reduction operator class (see
	 * <tt>RTOp_ROp_get_sub_vector.h</tt>) and calls
	 * <tt>applyOp()</tt>.  Note that the footprint of the reduction
	 * object (both internal and external state) will be
	 * O(<tt>rng.size()</tt>).  For serial applications this is faily
	 * reasonable and will not be a major performance penalty.  For
	 * parallel applications, however, this is a terrible
	 * implementation and must be overridden if <tt>rng.size()</tt> is
	 * large at all.  Although, this method should not even be used in
	 * case where the vector is very large.  If a subclass does
	 * override this method, it must also override
	 * <tt>freeSubSector()</tt> which has a default implementation
	 * which is a companion to this method's default implementation.
	 */
	virtual void getSubVector( const Range1D& rng, RTOpPack::SubVectorT<Scalar>* sub_vec ) const;

	///
	/** Free an explicit view of a sub-vector.
	 *
	 * @param  sub_vec
	 *				[in/out] The memory refered to by <tt>sub_vec->values</tt>
	 *				will be released if it was allocated and <tt>*sub_vec</tt>
	 *              will be zeroed out using <tt>RTOp_sub_vector_null(sub_vec)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>sub_vec</tt> must have been passed through a call to 
	 *      <tt>this->getSubVector(...,sub_vec)</tt>
	 * </ul>
 	 *
	 * Postconditions:<ul>
	 * <li> <tt>sub_vec->values == NULL</tt>, <tt>sub_vec->sub_dim == 0</tt> 
	 * </ul>
	 *
	 * The sub-vector view must have been allocated by this->getSubVector() first.
	 *
	 * This method has a default implementation which is a companion
	 * to the default implementation for <tt>getSubVector()</tt>.  If
	 * <tt>getSubVector()</tt> is overridden by a subclass then this
	 * method must be overridden also!
	 */
	virtual void freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const;

	///
	/** Get a mutable explicit view of a sub-vector.
	 *
	 * @param  rng      [in] The range of the elements to extract the sub-vector view.
	 * @param  sub_vec  [in/out] Mutable view of the sub-vector.  Prior to the
	 *                  first call <tt>RTOp_mutable_sub_vector_null(sub_vec)</tt> must
	 *                  have been called for the correct behavior.  Technically
	 *                  <tt>*sub_vec</tt> owns the memory but this memory can be freed
	 *                  only by calling <tt>this->commitSubVector(sub_vec)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() <= this->dim()) == true</tt>
	 *      (throw <tt>std::out_of_range</tt>)
	 * </ul>
 	 *
	 * Postconditions:<ul>
	 * <li> <tt>*sub_vec</tt> contains an explicit mutable view to the elements
	 *      in the range <tt>full_range(rng,1,this->dim())</tt>
	 * </ul>
	 *
	 * This is only a transient view of a sub-vector that is to be
	 * immediately used and then committed back with a call to
	 * <tt>commitSubVector()</tt>.
	 *
	 * Note that calling this operation might require some internal
	 * allocations and temporary memory.  Therefore, it is critical
	 * that <tt>this->commitSubVector(sub_vec)</tt> is called to
	 * commit the changed entires and clean up memory and avoid memory
	 * leaks after the sub-vector is modified.
	 *
	 * If <tt>this->getSubVector(...,sub_vec)</tt> was previously
	 * called on <tt>sub_vec</tt> then it may be possible to reuse
	 * this memory if it is sufficiently sized.  The user is
	 * encouraged to make multiple calls to
	 * <tt>this->getSubVector(...,sub_vec)</tt> before
	 * <tt>this->commitSubVector(sub_vec)</tt> to finally clean up all
	 * of the memory.  Of course the same <tt>sub_vec</tt> object must
	 * be passed to the same vector object for this to work correctly.
	 *
	 * Changes to the underlying sub-vector are not guarrenteed to
	 * become permanent until <tt>this->getSubVector(...,sub_vec)</tt>
	 * is called agian, or <tt>this->commitSubVector(sub_vec)</tt> is
	 * called.
	 *
	 * This method has a default implementation based on a vector
	 * reduction operator class (see
	 * <tt>RTOp_ROp_get_sub_vector.h</tt>) and calls
	 * <tt>applyOp()</tt>.  Note that the footprint of the reduction
	 * object (both internal and external state) will be
	 * O(<tt>rng.size()</tt>).  For serial applications this is faily
	 * reasonable and will not be a major performance penalty.  For
	 * parallel applications, this will be a terrible thing to do and
	 * must be overridden if <tt>rng.size()</tt> is large at all.  If
	 * a subclass does override this method, it must also override
	 * <tt>commitSubVector()</tt> which has a default implementation
	 * which is a companion to this method's default implementation.
	 *
	 */
	virtual void getSubVector( const Range1D& rng, RTOpPack::MutableSubVectorT<Scalar>* sub_vec );

	///
	/** Commit changes for a mutable explicit view of a sub-vector.
	 *
	 * @param sub_vec
	 *				[in/out] The data in <tt>sub_vec->values</tt> will be written
	 *              back internal storage and the memory refered to by
	 *              <tt>sub_vec->values</tt> will be released if it was allocated
	 *				and <tt>*sub_vec</tt> will be zeroed out using
	 *				<tt>RTOp_mutable_sub_vector_null(sub_vec)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>sub_vec</tt> must have been passed through a call to 
	 *      <tt>this->getSubVector(...,sub_vec)</tt>
	 * </ul>
 	 *
	 * Postconditions:<ul>
	 * <li> <tt>sub_vec->values == NULL</tt>, <tt>sub_vec->sub_dim == 0</tt> 
	 * </ul>
	 *
	 * The sub-vector view must have been allocated by
	 * <tt>this->getSubVector()</tt> first.
	 *
	 * This method has a default implementation which is a companion
	 * to the default implementation for <tt>getSubVector()</tt>.  If
	 * <tt>getSubVector()</tt> is overridden by a subclass then this
	 * method must be overridden also!
	 */
	virtual void commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec );

	///
	/** Set a specific sub-vector.
	 *
	 * After this function returns, the corresponding elements in
	 * <tt>this</tt> vector object will be set equal to those in the
	 * input vector.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>sub_vec.global_offset + sub_vec.sub_dim <= this->dim()</tt>
	 *      (<tt>throw std::out_of_range</tt>)
	 * </ul>
 	 *
	 * Postconditions:<ul>
	 * <li> All of the elements in the range
	 *      <tt>[sub_vec.global_offset+1,sub_vec.global_offset+sub_vec.sub_dim]</tt>
	 *      in <tt>*this</tt> are set to 0.0 except for those that have that
	 *      have entries in <tt>sub_vec</tt> which are set to the values spacificed
	 *      by <tt>(*this)(sub_vec.global_offset+vec.local_offset+sub_vec.indices[sub_vec.indices_stride*(k-1)])
	 *      = vec.values[vec.value_stride*(k-1)]</tt>, for <tt>k = 1..sub_vec.sub_nz</tt>
	 * </ul>
	 *
	 * The default implementation of this operation uses a transformation operator class
	 * (see RTOp_TOp_setSubVector.h) and calls <tt>applyOp()</tt>.  Be forewarned
	 * however, that the operator objects state data (both internal and external) will be 
	 * O(<tt>sub_vec.sub_nz</tt>).  For serial applications, this is entirely adequate.  For parallel
	 * applications this will be very bad!
	 *
	 * @param  sub_vec  [in] Represents the elements in the subvector to be set.
	 */
	virtual void setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec );

	//@}

#ifdef DOXYGEN_COMPILE
	VectorSpace<Scalar>*      space; // doxygen only!
#endif

}; // end class Vector

///
/** Apply a reduction/transformation,operation over a set of vectors:
 * <tt>op(op(v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)</tt>.
 *
 * @param  op	[in] Reduction/transformation operator to apply over each sub-vector
 *				and assemble the intermediate targets into <tt>reduct_obj</tt> (if
 *              <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
 * @param  num_vecs
 *				[in] Number of nonmutable vectors in <tt>vecs[]</tt>.
 *              If <tt>vecs==NULL</tt> then this argument is ignored but should be set to zero.
 * @param  vecs
 *				[in] Array (length <tt>num_vecs</tt>) of a set of pointers to
 *				nonmutable vectors to include in the operation.
 *				The order of these vectors is significant to <tt>op</tt>.
 * @param  num_targ_vecs
 *				[in] Number of mutable vectors in <tt>targ_vecs[]</tt>.
 *              If <tt>targ_vecs==NULL</tt>	then this argument is ignored but should be set to zero.
 * @param  targ_vecs
 *				[in] Array (length <tt>num_targ_vecs</tt>) of a set of pointers to
 *				mutable vectors to include in the operation.
 *				The order of these vectors is significant to <tt>op</tt>.
 *				If <tt>targ_vecs==NULL</tt> then <tt>op</tt> is called with no mutable vectors.
 * @param  reduct_obj
 *				[in/out] Target object of the reduction operation.
 *				This object must have been created by the <tt>op.reduct_obj_create_raw(&reduct_obj)</tt>
 *              function first.  The reduction operation will be added to <tt>(*reduct_obj)</tt> if
 *              <tt>(*reduct_obj)</tt> has already been through a reduction.  By allowing the info in
 *              <tt>(*reduct_obj)</tt> to be added to the reduction over all of these vectors, the reduction
 *              operation can be accumulated over a set of abstract vectors	which can be useful for implementing
 *              composite vectors for instance.  If <tt>op.get_reduct_type_num_entries(...)</tt> returns
 *              <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
 *              <tt>reduct_obj</tt> should be set to <tt>RTOp_REDUCT_OBJ_NULL</tt> and no reduction will be performed.
 * @param  first_ele
 *				[in] (default = 1) The index of the first element in <tt>this</tt> to be included.
 * @param  sub_dim
 *              [in] (default = 0) The number of elements in these vectors to include in the reduction/transformation
 *              operation.  The value of <tt>sub_dim == 0</tt> means to include all available elements.
 * @param  global_offset
 *				[in] (default = 0) The offset applied to the included vector elements.
 *
 *
 * Preconditions:<ul>
 * <li> [<tt>num_vecs > 0</tt>] <tt>vecs[k]->space()->isCompatible(*this->space()) == true</tt>,
 *      for <tt>k = 0...num_vecs-1</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
 * <li> [<tt>num_targ_vecs > 0</tt>] <tt>targ_vecs[k]->space()->isCompatible(*this->space()) == true</tt>,
 *      for <tt>k = 0...num_targ_vecs-1</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
 * <li> [<tt>num_targ_vecs > 0</tt>] The vectors pointed to by <tt>targ_vecs[k]</tt>, for
 *      <tt>k = 0...num_targ_vecs-1</tt> must not alias each other or any of the vectors
 *      <tt>vecs[k]</tt>, for <tt>k = 0...num_vecs-1</tt>.  <b>You have be warned!!!!</b>
 * <li> <tt>1 <= first_ele <= this->dim()</tt> (throw <tt>std::out_of_range</tt>)
 * <li> <tt>global_offset >= 0</tt> (throw <tt>std::invalid_argument</tt>)
 * <li> <tt>sub_dim - (first_ele - 1) <= this->dim()</tt> (throw <tt>std::length_error</tt>).
 * </ul>
 *
 * Postconditions:<ul>
 * <li> The vectors in <tt>targ_vecs[]</tt> may be modified as determined by the definition of <tt>op</tt>.
 * <li> [<tt>reduct_obj!=RTOp_REDUCT_OBJ_NULL</tt>] The reduction object <tt>reduct_obj</tt> contains
 *      the combined reduction from the input state of <tt>reduct_obj</tt> and the reductions that
 *      where accumulated during this this method invokation.
 * </ul>
 *
 * The logical vector passed to the
 * <tt>op\ref RTOpPack::RTOp::apply_op ".apply_op(...)"</tt>
 * method is: \verbatim

 v(k + global_offset) = this->get_ele(first_ele + k - 1)
 , for k = 1 ... sub_dim
 \endverbatim
 *
 * where <tt>v</tt> represents any one of the input or input/output
 * vectors.  The situation where <tt>first_ele == 1</tt> and
 * <tt>global_offset > 1</tt> corresponds to the case where the
 * vectors represent consitituent vectors in a larger aggregrate
 * vector.  The situation where <tt>first_ele > 1</tt> and
 * <tt>global_offset == 0</tt> is for when a sub-view of the vectors
 * are being treated as full vectors.  Other combinations of these
 * arguments are also possible.
 */
template<class Scalar>
inline
void applyOp(
	const RTOpPack::RTOpT<Scalar>   &op
	,const size_t                   num_vecs
	,const Vector<Scalar>*          vecs[]
	,const size_t                   num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOp_ReductTarget              reduct_obj
	,const Index                    first_ele     = 1
	,const Index                    sub_dim       = 0
	,const Index                    global_offset = 0
	)
{
	if(num_vecs)
		vecs[0]->applyOp(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
	else if (num_targ_vecs)
		targ_vecs[0]->applyOp(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
}


} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_DECL_HPP
