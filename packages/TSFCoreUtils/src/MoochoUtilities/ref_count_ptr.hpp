// ///////////////////////////////////////////////////////////////////////
// ref_count_ptr.hpp
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

#ifndef REF_COUNT_PTR_H
#define REF_COUNT_PTR_H

#include "ref_count_ptr_decl.hpp"

// /////////////////////////////////////////////////////////////////////////
// Inline implementations below, not for the client to look at.

namespace MemMngPack {

namespace PrivateUtilityPack {

// Assert that the pointer is not null
void assert_not_null(const void *);

// Node class to keep track of the delete address and
// the reference count for ref_count_ptr<...>
class ref_count_ptr_node {
public:
	//
	ref_count_ptr_node(bool has_ownership)
		: count_(1), has_ownership_(has_ownership)
	{}
	//
	virtual ~ref_count_ptr_node() {}
	//
	int count() const {
		return count_;	
	}
	//
	int incr_count() {
		return ++count_;
	}
	//
	int deincr_count() {
		return --count_;
	}
	//
	void has_ownership(bool has_ownership) {
		has_ownership_ = has_ownership;
	}
	//
	bool has_ownership() const {
		return has_ownership_;
	}

private:
	int         count_;
	bool        has_ownership_;
	// not defined and not to be called
	ref_count_ptr_node();
	ref_count_ptr_node(const ref_count_ptr_node&);
	ref_count_ptr_node& operator=(const ref_count_ptr_node&);
};	// end class ref_count_ptr_node;

// Implementation class for actually deleting the object if has_ownership() == true.
template<class T, class Dealloc_T>
class ref_count_ptr_node_tmpl : public ref_count_ptr_node {
public:

	//
	ref_count_ptr_node_tmpl(T* p, Dealloc_T dealloc, bool has_ownership)
		: ref_count_ptr_node(has_ownership), ptr_(p), dealloc_(dealloc)
	{}
	//
	const Dealloc_T& dealloc() const { return dealloc_; }
	//
	~ref_count_ptr_node_tmpl() {
		if( has_ownership() )
			dealloc_.free(ptr_);
	}

private:

	T           *ptr_;
	Dealloc_T   dealloc_;
	// not defined and not to be called
	ref_count_ptr_node_tmpl();
	ref_count_ptr_node_tmpl(const ref_count_ptr_node_tmpl&);
	ref_count_ptr_node_tmpl& operator=(const ref_count_ptr_node_tmpl&);

}; // end class ref_count_ptr_node_tmpl<T>

}	// end namespace PrivateUtilityPack 

// /////////////////////////////////////////////////////////////////////////////////
// Inline member functions for ref_count_ptr<...>.

template<class T>
inline
ref_count_ptr<T>::ref_count_ptr( ENull )
	: ptr_(NULL)
	, node_(NULL)
{}

template<class T>
inline
ref_count_ptr<T>::ref_count_ptr( T* p, bool has_ownership )
	: ptr_(p)
	, node_( p ? new PrivateUtilityPack::ref_count_ptr_node_tmpl<T,DeallocDelete<T> >(p,DeallocDelete<T>(),has_ownership) : NULL )
{}

#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
template<class T>
REF_COUNT_PTR_INLINE
template<class Dealloc_T>
ref_count_ptr<T>::ref_count_ptr( T* p, Dealloc_T dealloc, bool has_ownership )
	: ptr_(p)
	, node_( p ? new PrivateUtilityPack::ref_count_ptr_node_tmpl<T,Dealloc_T>(p,dealloc,has_ownership) : NULL )
{}
#endif

#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
template<class T>
REF_COUNT_PTR_INLINE
template <class T2>
ref_count_ptr<T>::ref_count_ptr(const ref_count_ptr<T2>& r_ptr)
	: ptr_(const_cast<T2*>(r_ptr.get()))                 // will not compile if T1 is not an ancestor of T2
	, node_(const_cast<node_t*>(r_ptr.access_node()))
{
	if(node_) node_->incr_count();
}
#endif

template<class T>
REF_COUNT_PTR_INLINE
ref_count_ptr<T>::ref_count_ptr(const ref_count_ptr<T>& r_ptr)
	: ptr_(r_ptr.ptr_), node_(r_ptr.node_)
{
	if(node_) node_->incr_count();
}

template<class T>
REF_COUNT_PTR_INLINE
ref_count_ptr<T>::~ref_count_ptr()
{
	if(node_ && node_->deincr_count() == 0 ) delete node_;
}

template<class T>
REF_COUNT_PTR_INLINE
ref_count_ptr<T>& ref_count_ptr<T>::operator=(const ref_count_ptr<T>& r_ptr) {
	if(node_) {
		if( r_ptr.node_ == node_ )
			return *this; // Assignment to self!
		if( !node_->deincr_count() ) {
			delete node_;
		}
	}
	ptr_  = r_ptr.ptr_;
	node_ = r_ptr.node_;
	if(node_) node_->incr_count();
	return *this;
}

template<class T>
inline
T* ref_count_ptr<T>::operator->() const {
	assert_not_null();
	return ptr_;
}

template<class T>
inline
T& ref_count_ptr<T>::operator*() const {
	assert_not_null();
	return *ptr_;
}

template<class T>
inline
T* ref_count_ptr<T>::get() const {
	return ptr_;
}

template<class T>
REF_COUNT_PTR_INLINE
T* ref_count_ptr<T>::release() {
	if(node_)
		node_->has_ownership(false);
	return ptr_;
}

template<class T>
REF_COUNT_PTR_INLINE
int ref_count_ptr<T>::count() const {
	if(node_)
		return node_->count();
	return 0;
}

template<class T>
REF_COUNT_PTR_INLINE
void ref_count_ptr<T>::set_has_ownership() {
	if(node_)
		node_->has_ownership(true);
}

template<class T>
REF_COUNT_PTR_INLINE
bool ref_count_ptr<T>::has_ownership() const {
	if(node_)
		return node_->has_ownership();
	return false;
}

template<class T>
REF_COUNT_PTR_INLINE
bool ref_count_ptr<T>::shares_resource(const ref_count_ptr<T>& r_ptr) const {
	return node_ == r_ptr.node_;
}

// private

template<class T>
inline
void ref_count_ptr<T>::assert_not_null() const {
	PrivateUtilityPack::assert_not_null(ptr_);
}

// very bad public functions

template<class T>
inline
ref_count_ptr<T>::ref_count_ptr( T* p, node_t* node)
	: ptr_(p), node_(node)
{
	if(node_) node_->incr_count();
}

template<class T>
inline
T*& ref_count_ptr<T>::access_ptr()
{	return ptr_; }

template<class T>
inline
typename ref_count_ptr<T>::node_t*& ref_count_ptr<T>::access_node()
{	return node_; }

template<class T>
inline
typename ref_count_ptr<T>::node_t* ref_count_ptr<T>::access_node() const
{	return node_; }

}	// end namespace MemMngPack

// /////////////////////////////////////////////////////////////////////////////////
// Inline member functions for conversions for ref_count_ptr<...>.

template<class T>
inline
MemMngPack::ref_count_ptr<T>
MemMngPack::rcp( T* p, bool owns_mem )
{
	return ref_count_ptr<T>(p,owns_mem);
}

#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
template<class T, class Dealloc_T>
inline
MemMngPack::ref_count_ptr<T>
MemMngPack::rcp( T* p, Dealloc_T dealloc, bool owns_mem )
{
	return ref_count_ptr<T>(p,dealloc,owns_mem);
}
#endif

template<class T2, class T1>
REF_COUNT_PTR_INLINE
MemMngPack::ref_count_ptr<T2>
MemMngPack::rcp_implicit_cast(const ref_count_ptr<T1>& p1)
{
	T2 *check = p1.get();	// Make the compiler check if the conversion is legal
	ref_count_ptr<T2> p2;
	if(p1.access_node()) {
		p2.access_ptr()  = check;
		p2.access_node() = const_cast<ref_count_ptr<T1>&>(p1).access_node();
		p2.access_node()->incr_count();
	}
	return p2;
}

template<class T2, class T1>
REF_COUNT_PTR_INLINE
MemMngPack::ref_count_ptr<T2>
MemMngPack::rcp_static_cast(const ref_count_ptr<T1>& p1)
{
	T2 *check = static_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
	ref_count_ptr<T2> p2;
	if(p1.access_node()) {
		p2.access_ptr()  = check;
		p2.access_node() = const_cast<ref_count_ptr<T1>&>(p1).access_node();
		p2.access_node()->incr_count();
	}
	return p2;
}

template<class T2, class T1>
REF_COUNT_PTR_INLINE
MemMngPack::ref_count_ptr<T2>
MemMngPack::rcp_const_cast(const ref_count_ptr<T1>& p1)
{
	T2 *check = const_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
	ref_count_ptr<T2> p2;
	if(p1.access_node()) {
		p2.access_ptr()  = check;
		p2.access_node() = const_cast<ref_count_ptr<T1>&>(p1).access_node();
		p2.access_node()->incr_count();
	}
	return p2;
}

template<class T2, class T1>
REF_COUNT_PTR_INLINE
MemMngPack::ref_count_ptr<T2>
MemMngPack::rcp_dynamic_cast(const ref_count_ptr<T1>& p1)
{
	ref_count_ptr<T2> p2; // NULL by default
	if( p1.get() ) {
		T2 *check = dynamic_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
		if(check) {
			p2.access_ptr()  = check;
			p2.access_node() = const_cast<ref_count_ptr<T1>&>(p1).access_node();
			p2.access_node()->incr_count();
		}
	}
	return p2;
}

#endif	// REF_COUNT_PTR_H
