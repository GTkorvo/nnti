// /////////////////////////////////////////////////////////////////////////////////////
// IterQuantityAccessContiguousDecl.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//
// Change Log:
//	11/18/99:
//		* last_updated() Added
//		* set_not_updated(k) Added
//		* resize(n) Added
//		* Lazy initialization implemented.
//  12/20/01:
//      * AbstractFactory added to handle memory management

#ifndef ITER_QUANITY_ACCESS_CONTIGUOUS_DECL_H
#define ITER_QUANITY_ACCESS_CONTIGUOUS_DECL_H

#include <vector>
#include <limits>

#include "IterQuantityAccess.hpp"
#include "AbstractFactoryStd.hpp"

namespace IterationPack {

// ToDo: Use an implementation subclass for the operations to avoid code blot.

///
/** Iteration Quanities subclass for contiguous iterations.
  *
  * This class implements the IterQuantityAccess interface for
  * the case where storage is provided for consecutive iterations.
  * This class allows the number of iterations to be variable on
  * construction.  For 3 storage locations, this class could
  * provide memory for the intervals [k-5, k-4, k-3], or [k-1, k, k+1] or
  * [k+1, k+2, k+3] etc.  The particular interval being represented is
  * dependent on the sequence in which <tt>set_k(offset)</tt> and <tt>next_iteration()</tt>
  * are called.  The only rule is that backward interation is not allowed.
  * For example, if the [k, k+1] interval is being represented then
  * <tt>set_k(-1)</tt> will throw an exception but <tt>set_k(+5)</tt> would not.
  * Backward memory is keep as long as possible and is only changed by <tt>set_k()</tt>.
  * For example, if the range [k-1, k] is being represented and <tt>set_k(+1)</tt> is called,
  * the interval [k, k+1] is now represented and the value for the kth iteration
  * is still keep but the value for the k-1 iteration is lost.  This
  * could have been determined a priori by calling <tt>will_loose_mem()</tt>. 
  * 
  * The default constructor, copy constructor and assignment operators are not
  * not allowed (i.e. declared private).
  */
template<class T_info>
class IterQuantityAccessContiguous : public IterQuantityAccess<T_info> {
public:

	///
	typedef MemMngPack::ref_count_ptr<
		const MemMngPack::AbstractFactory<T_info> >          abstract_factory_ptr_t;
	///
	typedef MemMngPack::AbstractFactoryStd<T_info,T_info>    abstract_factory_std_t;

	/** @name Constructors/initalizers */
	//@{

	///
	/** Construct storage for <tt>num_quantities</tt> with the name <tt>name</tt> given an abstract factory.
	 *
	 * After construction <tt>this->set_k(offset)</tt> can be called for any
	 * <tt>offset</tt> in the range of legal integers.
	 *
	 * Preconditions: <ul>
	 * <li> <tt>num_quantities > 0</tt> (throw std::length_error)
	 * </ul>
	 *
	 * If \c abstract_factory.get() == NULL then the client had better call \c this->set_factory()
	 * with an non-NULL factory before any attempt is made to call \c get_k() or \c set_k() or an
	 * exception will be thrown.
	 */
	IterQuantityAccessContiguous(
		int                              num_quantities
		,const std::string&              name
#ifdef _MIPS_CXX // MipsPro 7.3.1.1 tries to instantiate the default type even when one is specified?
		,const abstract_factory_ptr_t&   abstract_factory
#else
		,const abstract_factory_ptr_t&   abstract_factory  = MemMngPack::rcp(new abstract_factory_std_t())
#endif
		);

	///
	/** Set the abstract factory to use to allocate storate.
	 *
	 * Postconditions:<ul>
	 * <li> \c this will become uninitialized and current memory will be wipped out.
	 * </ul>
	 *
	 * If \c abstract_factory.get() == NULL then the client had better call \c this->set_factory()
	 * again later with a non-NULL factory before any attempt is made to call \c get_k() or
	 * \c set_k() or an exception will be thrown.
	 */
	void set_factory( const abstract_factory_ptr_t& abstract_factory );

	///
	/** Resize the number of contiguous storage locations.
	 *
	 * Postconditions:<ul>
	 * <li> \c this will become uninitialized and current memory will be wipped out.
	 * </ul>
	 */
	void resize( int num_quantities );

	///
	~IterQuantityAccessContiguous();

	//@}

	/** @name Access */
	//@{

	/// Return the number of continous storage locations
	int num_quantities() const;

	//@}

	/** @name Overridden from IterQuantity */
	//@{

	///
	IterQuantity* clone() const;
	///
	const char* name() const; 
	///
	bool has_storage_k(int offset) const;
	///
	bool updated_k(int offset) const;
	///
	int last_updated() const;
	///
	void set_not_updated_k(int offset);
	///
	void set_all_not_updated();
	///
	bool will_loose_mem(int offset, int set_offset) const;
	///
	void next_iteration();
	///
	void print_concrete_type( std::ostream& out ) const;

	//@}

	/** @name Overridden from IterQuantityAccess */
	//@{

	///
	T_info& get_k(int offset);
	///
	const T_info& get_k(int offset) const;
	///
	T_info& set_k(int offset);
	///
	T_info& set_k(int set_offset, int get_offset);

	//@}

private:

	// ///////////////////////////////////////////////
	// Private types

	//
	typedef std::vector<bool>                                                       updated_t;
	//
	typedef std::vector<typename abstract_factory_ptr_t::element_type::obj_ptr_t>   store_t;
	//
	typedef std::vector<T_info*>                                                    quantities_t;

	// ///////////////////////////////////////////////
	// Private data members

	// number of contigous iterations memory is reserved for.
	int	num_quantities_;
	// The name of the quantity (useful for debugging)
	std::string name_;
	// The abstract factory used to create the objects themselves
	abstract_factory_ptr_t   abstract_factory_;
	// The highest offset for iteration we are providing storage for.  We are providing
	// storage for iterations:
	//   [ k + max_offset_, k + max_offset_ - 1, ..., k + max_offset_ - num_quanities_ + 1 ].
	// For max_offset_ == 1 and num_quantities_  == 3:  [ k+1, k, k-1 ].
	int	max_offset_;	
	// Flags for if the iteration quanity was updated.
	//    updated_[max_offset - offset]
	//        for offset = max_offset, max_offset - 1, ..., max_offset_ - num_quanities_ + 1
	// returns true if and only if the quantity (k + offset) has been updated.
	updated_t     updated_;
	// Storage vector for the iteration quantities
	store_t       store_;
	// The vector of pointers to the storage quantities
	quantities_t  quantities_;

	// ///////////////////////////////////////////////
	// Private member functions

	// Returns true if storage is initialized and false otherwise.
	bool is_initialized() const;

	// Called to make sure that we are initialized before a
	// nonconst operation is performed.
	void lazy_initialization();

	// Called to release current memory
	void release_mem();

	// Not defined and not to be called
	IterQuantityAccessContiguous();
	IterQuantityAccessContiguous(const IterQuantityAccessContiguous&);
	IterQuantityAccessContiguous& operator=(const IterQuantityAccessContiguous&);

};	// end class IterQuantityAccessContiguous

// /////////////////////////////////////////////////////////////
// Inline members

template <class T_info>
inline
int IterQuantityAccessContiguous<T_info>::num_quantities() const
{	
	return num_quantities_; 
}

}	// end namespace IterationPack

#endif	// ITER_QUANITY_ACCESS_CONTIGUOUS_DECL_H
