// /////////////////////////////////////////////////////////
// WorkspacePack.hpp
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

#ifndef WORK_SPACE_PACK_H
#define WORK_SPACE_PACK_H

#include "ref_count_ptr.hpp"
#include "ThrowException.hpp"

namespace WorkspacePack {

class WorkspaceStore;
class RawWorkspace;

///
/** \defgroup WorkspacePack_grp Set of utilities for allocating temporary workspace.
 * \ingroup Misc_grp
 *
 * The goal of this set of utilities is to allow the user to create
 * arrays of uninitialized or default initialized objects as automatic
 * variables on the stack to be used for temporary workspace.
 *
 */
//@{

///
/** Pointer to global workspace object.
 *
 * This is a smart pointer to a workspace object can be set at any time
 * and will serve as the default workspace.  This object can serve as a
 * single workspace that can be used by all of the functions in 
 * an entire process thread for all of its temporary workspace memory needs.
 * By default this pointer is set to NULL and it is up to some entity
 * to set this pointer to a valid object.  If the application is to be
 * threaded, then sharing a single WorkspacePack::WorkspaceStore object between
 * threads will result in incorrect behavior and could potentially crash the
 * program in some cases.
 */
extern MemMngPack::ref_count_ptr<WorkspaceStore> default_workspace_store;

///
/** Encapulsation object for raw temporary workspace that has been allocated.
 * These objects can only be created on the stack and should not be included
 * as the member of any other classes.
 */
class RawWorkspace {
public:
	///
	friend class WorkspaceStore;
	///
	/** Allocate num_bytes bytes of temporary workspace.
	 * When this object is created if <tt>workspace_store != NULL</tt> the <tt>workspace_store</tt> object
	 * will be used to get the raw memory.  If <tt>workspace_store == NULL || </tt>
	 * <tt>workspace_store->num_bytes_remaining() < num_bytes</tt> then this memory
	 * will have to be dynamically allocated.
	 *
	 * Preconditons:<ul>
	 * <li> <tt>num_bytes >= 0</tt> (throw <tt>std::invalid_arguemnt)
	 * </ul>
	 *
	 * Postconditons:<ul>
	 * <li> <tt>this-></tt>num_bytes() == <tt>num_bytes</tt>
	 * <li> [<tt>num_bytes > 0</tt>] <tt>this-></tt>workspace_ptr() <tt>+ i</tt> for <tt>i = 0,..num_bytes-1</tt>
	 *      points to valid raw ininitialized allocated memory.
	 * <li> [<tt>num_bytes == 0</tt>] <tt>this-></tt>workspace_ptr() == NULL</tt>
	 * </ul>
	 *
	 * @param  workspace_store  [in] Pointer to the workspace object to get the memory from.
	 *                          This can be <tt>NULL</tt> in which case <tt>new T[]</tt> and 
	 *                          <tt>delete []</tt> will be used instead.
	 * @param  num_bytes        [in] The number of bytes to allocate.
	 */
	RawWorkspace(WorkspaceStore* workspace_store, size_t num_bytes);
	/// Deallocate workspace.
	~RawWorkspace();
	/// Return the number of bytes of raw workspace.
	size_t num_bytes() const;
	/// Give a raw pointer to the beginning of the workspace.
	char* workspace_ptr();
	///
	const char* workspace_ptr() const;
private:
	WorkspaceStore   *workspace_store_;
	char             *workspace_begin_;
	char             *workspace_end_;
	bool             owns_memory_;  // If true then the pointed to memory was allocated with
	                                // new so we need to call delete on it when we are destroyed.
	// not defined and not to be called
	RawWorkspace();
	RawWorkspace(const RawWorkspace&);
	RawWorkspace& operator=(const RawWorkspace&);
	static void* operator new(size_t);
	static void operator delete(void*);
}; // end class RawWorkspace

///
/** Templated class for workspace creation.
 *
 * Objects of this type are what should be created by the user
 * instead of RawWorkspace objects since this class will properly
 * initialize memory using placement new and allows typed operator[]
 * access to the array elements.  The default constructor,
 * copy constructor and assignment operations are not allowed and
 * objects can not be allocated with new.
 *
 * It is important to note that the constructors and destructors will
 * only be called if <tt>call_constructors=true</tt> (the default) is
 * passed to the consructor Workspace().  For build-in types that do not
 * need constructors and destructors called, the client should pass in
 * <tt>call_constructors=false</tt> .  Otherwise we would have to call
 * constructors and destructors on all of the memory and that could
 * considerably slow things down.
 *
 * With simple built-in data types (i.e. <tt>call_constructors=false</tt>)
 * the cost of creating and destroying
 * one of these objects should be O(1) independent of how much data
 * is requested.  This is true as long as no dynamic memory has to
 * be allocated (this is determined the object <tt>workspace_store</tt>
 * passed to the constructor Workspace()).
 */
template<class T>
class Workspace {
public:
	///
	/** Allocates a num_elements array of temporary objects.
	 *
	 * @param  workspace_store  [in] Pointer to the workspace object to get the memory from.
	 *                          This can be <tt>NULL</tt> in which case <tt>new T[]</tt> and 
	 *                          <tt>delete []</tt> will be used instead.
	 * @param  num_elements     [in] The number of bytes to allocate.
	 * @param  call_consructors [in] If <tt>true</tt> then constructors and destructors will be
	 *                          called on the allocated memory.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>num_element >= 0</tt> (throw <tt>std::invalid_argument)
	 * </ul>
	 *
	 * Postconditons:<ul>
	 * <li> <tt>this-></tt>size() == <tt>num_elements</tt>
	 * <li> [<tt>num_elements > 0</tt>] <tt>this-></tt>operator[i], for <tt>i = 0,..num_elements-1</tt>
	 *      points to valid allocated object of type <tt>T</tt>.
	 * <li> [num_elements > 0 && call_constructors==true</tt>] <tt>this-></tt>operator[i],
	 *      for <tt>i = 0,..num_elements-1</tt> was allocated as
	 *      <tt>new (&this->operator[i]) T()</tt>.
	 * </ul>
	 *
	 * When this object is created the <tt>workspace_store</tt> object
	 * will be used to get the raw memory if <tt>workspace_store != NULL</tt>.
	 * If <tt>workspace_store == NULL || workspace_store->num_bytes_remaining()</tt>
	 * <tt> < sizeof(T)*num_elements</tt> then this memory will have to be dynamically allocated.
	 * The memory is default initialized (or uninitialized) using placement new.  The
	 * constructor will only be called with placement new if <tt>call_constructor == ture</tt>.
	 * Otherwise, the memory will be left uninitlaized.  This is okay for integral types
	 * like <tt>double</tt> and <tt>int</tt> but not okay for class types like <tt>std::string</tt> etc.
	 */
	Workspace(WorkspaceStore* workspace_store, size_t num_elements, bool call_constructors = true);
	///
	/** The destructor on the elements will only be called if <tt>call_constructors == true</tt> was
	 * passed to the constructor.
	 */
	~Workspace();
	/// Return the number of elements in the array.
	size_t size() const;
	///
	/** Non-const zero based element access.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>0 <= i && i < size()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 */
	T& operator[](size_t i);
	///
	/** Const zero based element access.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>0 <= i && i < size()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 */
	const T& operator[](size_t i) const;
private:
	RawWorkspace  raw_workspace_;
	bool          call_constructors_;
	// not defined and not to be called
	Workspace();
	Workspace(const RawWorkspace&);
	Workspace& operator=(const RawWorkspace&);
	static void* operator new(size_t);
	static void operator delete(void*);
}; // end class Workspace

///
/** Workspace encapsulation class.
 *
 * Base class for objects that allocate a huge block of memory
 * at once and then allow RawWorkspace (an hense Workspace<T>) objects to be created
 * that make use of this memory in a stack-like fasion.  The classes WorkspaceStore
 * and RawWorkspace work closely together and are useless on their own.
 *
 * Through this interface, a client can not initialize or resize the size of the
 * available workspace and can not directly instantiate objects of this type.
 * Instead it must create a derived WorkspaceStoreInitializeable object defined later.
 */
class WorkspaceStore {
public:
	///
	friend class RawWorkspace;
	///
	~WorkspaceStore();
	///
	/** Return the total number of bytes that where initially allocated.
	 */
	size_t num_bytes_total() const;
	///
	/** Return the number of bytes remaining currently.
	 */
	size_t num_bytes_remaining() const;
	///
	/** Return the number of static memory allocations granted thus far.
	 * This is the number of memory allocations requested by the creation
	 * of RawWorkspace objects where there was sufficient preallocated memory
	 * to satisfy the request.
	 */
	int num_static_allocations() const;
	///
	/** Return the number of dynamic memory allocations granted thus far.
	 * This is the number of memory allocations requested by the creation
	 * of RawWorkspace objects where there was not sufficient preallocated memory
	 * to satisfy the request and dynamic memory had to be created.
	 */
	int num_dyn_allocations() const;
protected:
	///
	WorkspaceStore(size_t num_bytes);
	///
	void protected_initialize(size_t num_bytes);
private:
	char    *workspace_begin_; // Points to the beginning of raw allocated workspace.
	                           // If NULL then no workspace has been allocated yet.
	char    *workspace_end_;   // Points to one past the last byte of allocated workspace.
	                           // workspace_end_ >= workspace_begin_
	char    *curr_ws_ptr_;     // Points to the first available byte of workspace.
	                           // workspace_begin_ <= curr_ws_ptr_ <= workspace_end_
    int     num_static_allocations_; // Number of workspace allocation using already
	                           // allocated memory.
	int     num_dyn_allocations_; // Number of workspace allocations using dynamic
                               // memory because the current workspace store was
	                           // overridden
	// Not definted and not to be called
	WorkspaceStore(const WorkspaceStore&);
	WorkspaceStore& operator=(const WorkspaceStore&);
}; // end class WorkspaceStore

///
/** WorkspaceStore class that can be used to actually reinitialize memory.
 *
 * The client can create concrete instances of this type and initalize
 * the memory used.  The client should call <tt>initialize(num_bytes)</tt> to set the number
 * of bytes to allocate where <tt>num_bytes</tt> should be large enough to satisfy all but
 * the largests of memory request needs.
 */
class WorkspaceStoreInitializeable
	: public WorkspaceStore
{
public:
	///
	/** Default constructs to no memory set and will dynamically
	 * allocate all memory requested.
	 */
	WorkspaceStoreInitializeable(size_t num_bytes = 0);
	///
	/** Set the size block of memory to be given as workspace.
	 *
	 * If there are any instantiated RawWorkspace objects then this
	 * function willl throw an exception.  It must be called before
	 * any RawWorkspace objects are created.
	 */
	void initialize(size_t num_bytes);
}; // end class WorkspaceStoreInitializeable

//@}

// /////////////////////////////////////
// Inline members for Workspace<T>

template<class T>
//inline
Workspace<T>::Workspace(WorkspaceStore* workspace_store, size_t num_elements, bool call_constructors)
	: raw_workspace_(workspace_store,sizeof(T)*num_elements), call_constructors_(call_constructors)
{
	if(call_constructors_) {
		char* raw_ptr = raw_workspace_.workspace_ptr();
		for( size_t k = 0; k < num_elements; ++k, raw_ptr += sizeof(T) )
			::new (raw_ptr) T(); // placement new
	}
}

template<class T>
//inline
Workspace<T>::~Workspace()
{
	if(call_constructors_) {
		const size_t num_elements = this->size();
		char* raw_ptr = raw_workspace_.workspace_ptr();
		for( size_t k = 0; k < num_elements; ++k, raw_ptr += sizeof(T) )
			reinterpret_cast<T*>(raw_ptr)->~T();
	}
}

template<class T>
inline
size_t Workspace<T>::size() const
{
	return raw_workspace_.num_bytes() / sizeof(T);
}

template<class T>
inline
T& Workspace<T>::operator[](size_t i)
{
#ifdef _DEBUG
	THROW_EXCEPTION( !( 0 <= i && i < this->size() ), std::invalid_argument, "Workspace<T>::operator[](i): Error!" );
#endif	
	return reinterpret_cast<T*>(raw_workspace_.workspace_ptr())[i];
}

template<class T>
inline
const T& Workspace<T>::operator[](size_t i) const
{
	return const_cast<Workspace<T>*>(this)->operator[](i);
}

#ifdef _PG_CXX // Should not have to define this but pgCC is complaining!
template<class T>
inline
void* Workspace<T>::operator new(size_t)
{
	assert(0);
	return NULL;
}
#endif

// should not have to define this but the gcc-2.95.2 compiler is complaining!
template<class T>
inline
void Workspace<T>::operator delete(void*)
{
	assert(0);
}

// /////////////////////////////////////
// Inline members for WorkspaceStore

inline
size_t WorkspaceStore::num_bytes_total() const
{
	return workspace_end_ - workspace_begin_;
}

inline
size_t WorkspaceStore::num_bytes_remaining() const
{
	return workspace_end_ - curr_ws_ptr_;
}

inline
int WorkspaceStore::num_static_allocations() const
{
	return num_static_allocations_;
}

inline
int WorkspaceStore::num_dyn_allocations() const
{
	return num_dyn_allocations_;
}

// /////////////////////////////////////////////////
// Inline members for WorkspaceStoreInitializeable

inline
WorkspaceStoreInitializeable::WorkspaceStoreInitializeable(size_t num_bytes)
	: WorkspaceStore(num_bytes)
{}

inline
void WorkspaceStoreInitializeable::initialize(size_t num_bytes)
{
	protected_initialize(num_bytes);
}

// /////////////////////////////////////
// Inline members for RawWorkspace

inline
size_t RawWorkspace::num_bytes() const
{
	return workspace_end_ - workspace_begin_;
}

inline
char* RawWorkspace::workspace_ptr()
{
	return workspace_begin_;
}

inline
const char* RawWorkspace::workspace_ptr() const
{
	return workspace_begin_;
}

// should not have to define this but the gcc-2.95.2 compiler is complaining!
inline
void RawWorkspace::operator delete(void*)
{
	assert(0);
}

} // end namespace WorkspacePack

#endif // WORK_SPACE_PACK_H
