// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_REFCOUNTPTR_DECL_H
#define TEUCHOS_REFCOUNTPTR_DECL_H

/*! \file Teuchos_RefCountPtrDecl.hpp
    \brief Reference-counted pointer class and non-member templated function
	definitions
*/

#include "Teuchos_any.hpp"

#ifdef REFCOUNTPTR_INLINE_FUNCS
#define REFCOUNTPTR_INLINE inline
#else
#define REFCOUNTPTR_INLINE
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace MemMngPack {} // ToDo: Take out latter!
#endif

/** \class Teuchos::DeallocDelete
    \brief Policy class for deallocator that uses <tt>delete</tt> to delete a
    pointer which is used by <tt>RefCountPtr</tt>.
 */

/** \class Teuchos::RefCountPtr
    \brief Templated class for reference counted smart pointers (see \ref RefCountPtr_stuff "description").
 *
 * This is a class for smart reference counted pointer objects
 * that deletes an object (if the object is owned) that was allocated
 * by <tt>new</tt> after all refereces to it have been removed.
 *
 * To see how to use this class see \ref RefCountPtr_stuff "description".
 */

namespace Teuchos {

namespace PrivateUtilityPack {
	class RefCountPtr_node;
}

/** \defgroup RefCountPtr_stuff Reference counting smart pointer class for automatic garbage collection.
  
For a carefully written discussion about what this class is and basic
details on how to use it see the
<A HREF="http://software.sandia.gov/Trilinos/RefCountPtrBeginnersGuideSAND.pdf">beginners guide</A>.

<b>Quickstart for <tt>RefCountPtr</tt></b>
 
Here we present a short, but fairly comprehensive, quick-start for the
use of <tt>RefCountPtr<></tt>.  The use cases described here
should cover the overwhelming majority of the use instances of
<tt>RefCountPtr<></tt> in a typical program.

The following class hierarchy will be used in the C++ examples given
below.

\code

class A { public: virtual ~A(){} virtual void f(){} };
class B1 : virtual public A {};
class B2 : virtual public A {};
class C : virtual public B1, virtual public B2 {};

class D {};
class E : public D {};

\endcode

All of the following code examples used in this quickstart are assumed
to be in the namespace <tt>Teuchos</tt> or have appropriate <tt>using
Teuchos::...</tt> declarations.  This removes the need to explicitly
use <tt>Teuchos::</tt> to qualify classes, functions and other
declarations from the <tt>Teuchos</tt> namespace.  Note that some of
the runtime checks are denoted as "debug runtime checked" which means
that checking will only be performed in a debug build (that is one
where the macro _DEBUG is defined at compile time).

<ol>

<li> <b>Creation of <tt>RefCountPtr<></tt> objects</b>

<ol>

<li> <b>Creating a <tt>RefCountPtr<></tt> object using <tt>new</tt></b>

\code
RefCountPtr<C> c_ptr = rcp(new C);
\endcode

<li> <b>Creating a <tt>RefCountPtr<></tt> object to an array allocated using <tt>new[n]</tt></b>

\code

RefCountPtr<C> c_ptr = rcp(new C[n],DeallocArrayDelete<C>(),true);
\endcode

<li> <b>Initializing a <tt>RefCountPtr<></tt> object to NULL</b>

\code
RefCountPtr<C> c_ptr;
\endcode
or
\code
RefCountPtr<C> c_ptr = null;
\endcode

<li> <b>Initializing a <tt>RefCountPtr<></tt> object to an object
       \underline{not} allocated with <tt>new</tt></b>

\code
C              c;
RefCountPtr<C> c_ptr = rcp(&c,false);
\endcode


<li> <b>Representing constantness and non-constantness</b>

<ol>

<li> <b>Non-constant pointer to non-constant object</b>
\code
RefCountPtr<C> c_ptr;
\endcode

<li> <b>Constant pointer to non-constant object</b>
\code
const RefCountPtr<C> c_ptr;
\endcode

<li> <b>Non-Constant pointer to constant object</b>
\code
RefCountPtr<const C> c_ptr;
\endcode

<li> <b>Constant pointer to constant object</b>
\code
const RefCountPtr<const C> c_ptr;
\endcode

</ol>

<li> <b>Copy constructor (implicit casting)</b>

\code
RefCountPtr<C>       c_ptr  = rcp(new C); // No cast
RefCountPtr<A>       a_ptr  = c_ptr;      // Cast to base class
RefCountPtr<const A> ca_ptr = a_ptr;      // Cast from non-const to const
\endcode

</ol>

<li> <b>Reinitialization of <tt>RefCountPtr<></tt> objects (using assignment operator)</b>

<ol>

<li> <b>Resetting from a raw pointer</b>

\code
RefCountPtr<A> a_ptr;
a_ptr = rcp(new C());
\endcode

<li> <b>Resetting to null</b>

\code
RefCountPtr<A> a_ptr = rcp(new C());
a_ptr = null; // The C object will be deleted here
\endcode

<li> <b>Assigning from a <tt>RefCountPtr<></tt> object</b>

\code
RefCountPtr<A> a_ptr1;
RefCountPtr<A> a_ptr2 = rcp(new C());
a_ptr1 = a_ptr2; // Now a_ptr1 and a_ptr2 point to same C object
\endcode

</ol>

<li> <b>Accessing the reference-counted object</b>

<ol>

<li> <b>Access to object reference (debug runtime checked)</b> : <tt>RefCountPtr::operator*()</tt> 

\code
C &c_ref = *c_ptr;
\endcode

<li> <b>Access to object pointer (unchecked, may return <tt>NULL</tt>)</b> : <tt>RefCountPtr::get()</tt>

\code
C *c_rptr = c_ptr.get();
\endcode

<li> <b>Access to object pointer (debug runtime checked, will not return <tt>NULL</tt>)</b> : <tt>RefCountPtr::operator*()</tt>

\code
C *c_rptr = &*c_ptr;
\endcode

<li> <b>Access of object's member (debug runtime checked)</b> : <tt>RefCountPtr::operator->()</tt>

\code
c_ptr->f();
\endcode

<li> <b>Testing for null</b>

\code
if( !a_ptr.get() ) std::cout << "a_ptr is null!\n";
\endcode

or

\code
if( a_ptr == null ) std::cout << "a_ptr is null!\n";
\endcode

or

\code
if( is_null(a_ptr) ) std::cout << "a_ptr is null!\n";
\endcode

</ol>

<li> <b>Casting</b>

<ol>

<li> <b>Implicit casting (see copy constructor above)</b>

<li> <b>Casting away <tt>const</tt></b> : <tt>rcp_const_cast()</tt>

\code
RefCountPtr<const A>  ca_ptr = rcp(new C);
RefCountPtr<A>        a_ptr  = rcp_const_cast<A>(ca_ptr); // cast away const!
\endcode

<li> <b>Static cast (no runtime check)</b> : <tt>rcp_static_cast()</tt>

\code
RefCountPtr<D>     d_ptr = rcp(new E);
RefCountPtr<E>     e_ptr = rcp_static_cast<E>(d_ptr); // Unchecked, unsafe?
\endcode

<li> <b>Dynamic cast (runtime checked, failed cast allowed)</b> : <tt>rcp_dynamic_cast()</tt>

\code
RefCountPtr<A>     a_ptr  = rcp(new C);
RefCountPtr<B1>    b1_ptr = rcp_dynamic_cast<B1>(a_ptr);  // Checked, safe!
RefCountPtr<B2>    b2_ptr = rcp_dynamic_cast<B2>(b1_ptr); // Checked, safe!
RefCountPtr<C>     c_ptr  = rcp_dynamic_cast<C>(b2_ptr);  // Checked, safe!
\endcode

<li> <b>Dynamic cast (runtime checked, failed cast not allowed)</b> : <tt>rcp_dynamic_cast()</tt>

\code
RefCountPtr<A>     a_ptr1  = rcp(new C);
RefCountPtr<A>     a_ptr2  = rcp(new A);
RefCountPtr<B1>    b1_ptr1 = rcp_dynamic_cast<B1>(a_ptr1,true);  // Success!
RefCountPtr<B1>    b1_ptr2 = rcp_dynamic_cast<B1>(a_ptr2,true);  // Throw std::bad_cast!
\endcode

</ol>

<li> <b>Managing extra data</b>

<ol>

<li> <b>Adding extra data (post destruction of extra data)</b> : <tt>set_extra_data()</tt>

\code
set_extra_data(rcp(new B1),"A:B1",&a_ptr);
\endcode

<li> <b>Adding extra data (pre destruction of extra data)</b> : <tt>get_extra_data()</tt>

\code
set_extra_data(rcp(new B1),"A:B1",&a_ptr,PRE_DESTORY);
\endcode

<li> <b>Retrieving extra data</b> : <tt>get_extra_data()</tt>

\code
get_extra_data<RefCountPtr<B1> >(a_ptr,"A:B1")->f();
\endcode

<li> <b>Resetting extra data</b> : <tt>get_extra_data()</tt>

\code
get_extra_data<RefCountPtr<B1> >(a_ptr,"A:B1") = rcp(new C);
\endcode

</ol>

</ol>

 */
//@{

/// Used to initialize a <tt>RefCountPtr</tt> object to NULL using an implicit conversion!
enum ENull { null };

/// Used to specify a pre or post destruction of extra data
enum EPrePostDestruction { PRE_DESTROY, POST_DESTROY };

/// Deallocator for <tt>new</tt> which calls <tt>delete</tt>
template<class T>
class DeallocDelete
{
public:
	/// Gives the type (required)
	typedef T ptr_t;
	/// Deallocates a pointer <tt>ptr</tt> using <tt>delete ptr</tt> (required).
	void free( T* ptr ) { if(ptr) delete ptr; }
};

/// Deallocator for <tt>new []</tt> which calls <tt>delete []</tt>
template<class T>
class DeallocArrayDelete
{
public:
	/// Gives the type (required)
	typedef T ptr_t;
	/// Deallocates a pointer <tt>ptr</tt> using <tt>delete [] ptr</tt> (required).
	void free( T* ptr ) { if(ptr) delete [] ptr; }
};

/** \brief . */
template<class T>
class RefCountPtr {
public:
	/** \brief . */
	typedef T	element_type;
	/** \brief Initialize <tt>RefCountPtr<T></tt> to NULL.
	 *
	 * This allows clients to write code like:
	 \code
	 RefCountPtr<int> p = null;
	 \endcode
	 or
	 \code
	 RefCountPtr<int> p;
	 \endcode
	 * and construct to <tt>NULL</tt>
	 */
	RefCountPtr( ENull null_arg = null );
	/** \brief Initialize from another <tt>RefCountPtr<T></tt> object.
	 *
	 * After construction, <tt>this</tt> and <tt>r_ptr</tt> will
	 * reference the same object.
	 *
	 * This form of the copy constructor is required even though the
	 * below more general templated version is sufficient since some
	 * compilers will generate this function automatically which will
	 * give an incorrect implementation.
	 *
	 * Postconditons:<ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	RefCountPtr(const RefCountPtr<T>& r_ptr);
	/** \brief Initialize from another <tt>RefCountPtr<T2></tt> object (implicit conversion only).
	 *
	 * This function allows the implicit conversion of smart pointer objects just
	 * like with raw C++ pointers.  Note that this function will only compile
	 * if the statement <tt>T1 *ptr = r_ptr.get()</tt> will compile.
	 *
	 * Postconditons: <ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	template<class T2>
	RefCountPtr(const RefCountPtr<T2>& r_ptr);
	/** \brief Removes a reference to a dynamically allocated object and possibly deletes
	 * the object if owned.
	 *
	 * Peforms <tt>delete ...</tt> if <tt>this->has_ownership() == true</tt> and
	 * <tt>this->count() == 1</tt>.  If <tt>this->count() == 1</tt> but <tt>this->has_ownership() == false</tt>
	 * then the object is not deleted.
	 * If <tt>this->count() > 1</tt> then then internal reference count
	 * shared by all the other related <tt>RefCountPtr<...></tt> objects for this shared
	 * object is deincremented by one.  If <tt>this->get() == NULL</tt> then nothing happens.
	 */
	~RefCountPtr();
	/** \brief Copy the pointer to the referenced object and increment the reference count.
	 *
	 * If <tt>this->has_ownership() == true</tt> and <tt>this->count() == 1</tt> before this operation
	 * is called, then the object pointed to by <tt>this->get()</tt> will be deleted (using <tt>delete</tt>)
	 * prior to binding to the pointer (possibly <tt>NULL</tt>) pointed to in <tt>r_ptr</tt>.
	 * Assignment to self (i.e. <tt>this->get() == r_ptr.get()</tt>) is harmless and this
	 * function does nothing.
	 *
	 * Postconditons:
	 * <ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	RefCountPtr<T>& operator=(const RefCountPtr<T>& r_ptr);
	/** \brief Pointer (<tt>-></tt>) access to members of underlying object.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>this->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
	 * </ul>
	 */
	T* operator->() const;
	/** \brief Dereference the underlying object.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>this->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
	 * </ul>
	 */
	T& operator*() const;
    /** \brief Get the raw C++ pointer to the underlying object.
	 */
	T* get() const;
	/** \brief Release the ownership of the underlying dynamically allocated object.
	 *
	 * After this function is called then the client is responsible for calling
	 * delete on the returned pointer no matter how many <tt>ref_count_prt<T></tt> objects
	 * have a reference to it.  If <tt>this-></tt>get() <tt>== NULL</tt>, then this call is
	 * meaningless.
	 *
	 * Note that this function does not have the exact same semantics as does
	 * <tt>auto_ptr<T>::release()</tt>.  In <tt>auto_ptr<T>::release()</tt>, <tt>this</tt>
	 * is set to <tt>NULL</tt> while here in RefCountPtr<T>:: release() only an ownership flag is set
	 * and <tt>this</tt> still points to the same object.  It would be difficult to duplicate
	 * the behavior of <tt>auto_ptr<T>::release()</tt> for this class.
	 *
	 * Postconditions:
	 * <ul>
	 * <li> <tt>this->has_ownership() == false</tt>
	 * </ul>
	 *
	 * @return Returns the value of <tt>this->get()</tt>
	 */
	T* release();
	/** \brief Return the number of <tt>RefCountPtr<></tt> objects that have a reference
	 * to the underlying pointer that is being shared.
	 *
	 * @return  If <tt>this->get() == NULL</tt> then this function returns 0.
	 * Otherwise, this function returns <tt>> 0</tt>.
	 */
	int count() const;
	/** \brief Give <tt>this</tt> and other <tt>RefCountPtr<></tt> objects ownership 
	 * of the referenced object <tt>this->get()</tt>.
	 *
	 * See ~RefCountPtr() above.  This function
	 * does nothing if <tt>this->get() == NULL</tt>.
	 *
	 * Postconditions:
	 * <ul>
	 * <li> If <tt>this->get() == NULL</tt> then
	 *   <ul>
	 *   <li> <tt>this->has_ownership() == false</tt> (always!).
	 *   </ul>
	 * <li> else
	 *   <ul>
	 *   <li> <tt>this->has_ownership() == true</tt>
	 *   </ul>
	 * </ul>
	 */
	void set_has_ownership();
	/** \brief Returns true if <tt>this</tt> has ownership of object pointed to by <tt>this->get()</tt> in order to delete it.
	 *
	 * See ~RefCountPtr() above.
	 *
	 * @return If this->get() <tt>== NULL</tt> then this function always returns <tt>false</tt>.
	 * Otherwise the value returned from this function depends on which function was
	 * called most recently, if any; set_has_ownership() (<tt>true</tt>)
	 * or release() (<tt>false</tt>).
	 */
	bool has_ownership() const;
	/** \brief Returns true if the smart pointers share the same underlying reference-counted object.
	 *
	 * This method does more than just check if <tt>this->get() == r_ptr.get()</tt>.
	 * It also checks to see if the underlying reference counting machinary is the
	 * same.
	 */
	bool shares_resource(const RefCountPtr<T>& r_ptr) const;
	/** \brief Throws <tt>std::logic_error</tt> if <tt>this->get()==NULL</tt>, otherwise returns reference to <tt>*this</tt>. */
	const RefCountPtr<T>& assert_not_null() const;

private:

	// //////////////////////////////////////
	// Private types

	typedef PrivateUtilityPack::RefCountPtr_node			node_t;

	// //////////////////////////////////////////////////////////////
	// Private data members

	T       *ptr_;  // NULL if this pointer is null
	node_t	*node_;	// NULL if this pointer is null

public:
#ifndef DOXYGEN_COMPILE
	// These constructors should be private but I have not had good luck making
	// this portable (i.e. using friendship etc.) in the past
	RefCountPtr( T* p, bool has_ownership );
	template<class Dealloc_T>
	RefCountPtr( T* p, Dealloc_T dealloc, bool has_ownership );
	// This is a very bad breach of encapsulation that is needed since MS VC++ 5.0 will
	// not allow me to declare template functions as friends.
	RefCountPtr( T* p, node_t* node);
	T*&           access_ptr();
	node_t*&      access_node();
	node_t*       access_node() const;
#endif

};	// end class RefCountPtr<...>

/** \brief Create a <tt>RefCountPtr</tt> object properly typed.
 *
 * @param  p  [in] Pointer to an object to be reference counted.
 * @param owns_mem
 *            [in] If <tt>owns_mem==true</tt>  then <tt>delete p</tt>
 *            will be called when the last reference to this object
 *            is removed.  If <tt>owns_mem==false</tt> then nothing
 *            will happen to delete the the object pointed to by
 *            <tt>p</tt> when the last reference is removed.
 *
 * Preconditions:<ul>
 * <li> If <tt>owns_mem==true</tt> then <tt>p</tt> must have been
 *      created by calling <tt>new</tt> to create the object since
 *      <tt>delete p</tt> will be called eventually.
 * </ul>
 *
 * If the pointer <tt>p</tt> did not come from <tt>new</tt> then
 * either the client should use the version of <tt>rcp()</tt> that
 * that uses a deallocator policy object or should pass in 
 * <tt>owns_mem = false</tt>.
 */
template<class T>
RefCountPtr<T> rcp( T* p, bool owns_mem
#ifndef __sun
	= true
#endif
	);
#ifdef __sun // RAB: 20040303: Sun needs to fix their compiler
template<class T> inline RefCountPtr<T> rcp( T* p ) { return rcp(p,true); }
#endif

/** \brief Initialize from a raw pointer with a deallocation policy.
 *
 * @param  p       [in] Raw C++ pointer that \c this will represent.
 * @param  dealloc [in] Deallocator policy object (copied by value) that defines
 *                 a function <tt>void Dealloc_T::free(T* p)</tt> that will
 *                 free the underlying object.
 * @param  owns_mem
 *                 [in] If true then <tt>return</tt> is allowed to delete
 *                 the underlying pointer by calling <tt>dealloc.free(p)</tt>.
 *                 when all references have been removed.
 *
 * Postconditions:<ul>
 * <li> The function <tt>void Dealloc_T::free(T* p)</tt> exists.
 * </ul>
 *
 * Postconditions:<ul>
 * <li> <tt>return.get() == p</tt>
 * <li> If <tt>p == NULL</tt> then
 *   <ul>
 *   <li> <tt>return.count() == 0</tt>
 *   <li> <tt>return.has_ownership() == false</tt>
 *   </ul>
 * <li> else
 *   <ul>
 *   <li> <tt>return.count() == 1</tt>
 *   <li> <tt>return.has_ownership() == owns_mem</tt>
 *   </ul>
 * </ul>
 *
 * By default, <tt>return</tt> has ownership to delete the object
 * pointed to by <tt>p</tt> when <tt>return</tt> is deleted (see
 * <tt>~RefCountPtr())</tt>.  It is vitually important that if
 * <tt>owns_mem == true</tt> that the address <tt>p</tt> that is
 * passed in is the same address that was returned by <tt>new</tt>.
 * With multiple inheritance this is not always the case.  See the
 * above discussion.  This class is templated to accept a deallocator
 * object that will free the pointer.  The other functions use a
 * default deallocator of type <tt>DeallocDelete</tt> which has a method
 * <tt>DeallocDelete::free()</tt> which just calls <tt>delete p</tt>.
 */
template<class T, class Dealloc_T>
RefCountPtr<T> rcp( T* p, Dealloc_T dealloc, bool owns_mem );

/** \brief Returns true if <tt>p.get()==NULL</tt>.
 */
template<class T>
bool is_null( const RefCountPtr<T> &p );

/** \brief Returns true if <tt>p.get()==NULL</tt>.
 */
template<class T>
bool operator==( const RefCountPtr<T> &p, ENull );

/** \brief Returns true if <tt>p.get()!=NULL</tt>.
 */
template<class T>
bool operator!=( const RefCountPtr<T> &p, ENull );

/** \brief Return true if two <tt>RefCountPtr</tt> objects point to the same
 * referenced-counted object.
 */
template<class T1, class T2>
bool operator==( const RefCountPtr<T1> &p1, const RefCountPtr<T2> &p2 );

/** \brief Return true if two <tt>RefCountPtr</tt> objects do not point to the
 * same referenced-counted object.
 */
template<class T1, class T2>
bool operator!=( const RefCountPtr<T1> &p1, const RefCountPtr<T2> &p2 );

/** \brief Implicit cast of underlying <tt>RefCountPtr</tt> type from <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * The function will compile only if (<tt>T2* p2 = p1.get();</tt>) compiles.
 *
 * This is to be used for conversions up an inheritance hierarchy and from non-const to
 * const and any other standard implicit pointer conversions allowed by C++.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_implicit_cast(const RefCountPtr<T1>& p1);

/** \brief Static cast of underlying <tt>RefCountPtr</tt> type from <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * The function will compile only if (<tt>static_cast<T2*>(p1.get());</tt>) compiles.
 *
 * This can safely be used for conversion down an inheritance hierarchy
 * with polymorphic types only if <tt>dynamic_cast<T2>(p1.get()) == static_cast<T2>(p1.get())</tt>.
 * If not then you have to use <tt>rcp_dynamic_cast<tt><T2>(p1)</tt>.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_static_cast(const RefCountPtr<T1>& p1);

/** \brief Constant cast of underlying <tt>RefCountPtr</tt> type from <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * This function will compile only if (<tt>const_cast<T2*>(p1.get());</tt>) compiles.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_const_cast(const RefCountPtr<T1>& p1);

/** \brief Dynamic cast of underlying <tt>RefCountPtr</tt> type from <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * @param  p1             [in] The smart pointer casting from
 * @param  throw_on_fail  [in] If <tt>true</tt> then if the cast fails (for <tt>p1.get()!=NULL) then
 *                        a <tt>std::bad_cast</tt> exception is thrown with a very informative
 *                        error message.
 *
 * Postconditions:<ul>
 * <li> If <tt>( p1.get()!=NULL && throw_on_fail==true && dynamic_cast<T2*>(p1.get())==NULL ) == true</tt>
 *      then an <tt>std::bad_cast</tt> exception is thrown with a very informative error message.
 * <li> If <tt>( p1.get()!=NULL && dynamic_cast<T2*>(p1.get())!=NULL ) == true</tt>
 *      then <tt>return.get() == dynamic_cast<T2*>(p1.get())</tt>.
 * <li> If <tt>( p1.get()!=NULL && throw_on_fail==false && dynamic_cast<T2*>(p1.get())==NULL ) == true</tt>
 *      then <tt>return.get() == NULL</tt>.
 * <li> If <tt>( p1.get()==NULL ) == true</tt>
 *      then <tt>return.get() == NULL</tt>.
 * </ul>
 *
 * This function will compile only if (<tt>dynamic_cast<T2*>(p1.get());</tt>) compiles.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_dynamic_cast(
	const RefCountPtr<T1>& p1
	,bool throw_on_fail
#ifndef __sun
	= false
#endif
	);
#ifdef __sun // RAB: 20041019: Sun needs to fix their compiler
template<class T2, class T1> inline RefCountPtr<T2> rcp_dynamic_cast( const RefCountPtr<T1>& p1 )
{ return rcp_dynamic_cast<T2>(p1,false); }
#endif

/** \brief Set extra data associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  extra_data
 *               [in] Data object that will be set (copied)
 * @param  name  [in] The name given to the extra data.  The value of
 *               <tt>name</tt> together with the data type <tt>T1</tt> of the
 *               extra data must be unique from any other such data or
 *               the other data will be overwritten.
 * @param  p     [out] On output, will be updated with the input <tt>extra_data</tt>
 * @param  force_unique
 *               [in] Determines if this type and name pair must be unique
 *               in which case if an object with this same type and name
 *               already exists, then an exception will be thrown.
 *               The default is <tt>true</tt> for safety.
 * @param  destroy_when
 *               [in] Determines when <tt>extra_data</tt> will be destoryed
 *               in relation to the underlying reference-counted object.
 *               If <tt>destroy_when==PRE_DESTROY</tt> then <tt>extra_data</tt>
 *               will be deleted before the underlying reference-counted object.
 *               If <tt>destroy_when==POST_DESTROY</tt> (the default) then <tt>extra_data</tt>
 *               will be deleted after the underlying reference-counted object.
 *
 * If there is a call to this function with the same type of extra
 * data <tt>T1</tt> and same arguments <tt>p</tt> and <tt>name</tt>
 * has already been made, then the current piece of extra data already
 * set will be overwritten with <tt>extra_data</tt>.  However, if the
 * type of the extra data <tt>T1</tt> is different, then the extra
 * data can be added and not overwrite existing extra data.  This
 * means that extra data is keyed on both the type and name.  This
 * helps to minimize the chance that clients will unexpectedly
 * overwrite data by accident.
 *
 * When the last <tt>RefcountPtr</tt> object is removed and the
 * reference-count node is deleted, then objects are deleted in the following
 * order: (1) All of the extra data that where added with
 * <tt>destroy_when==PRE_DESTROY</tt> are first, (2) then the underlying
 * reference-counted object is deleted, and (3) the rest of the extra data
 * that was added with <tt>destroy_when==PRE_DESTROY</tt> is then deleted.
 * The order in which the objects are destroyed is not guaranteed.  Therefore,
 * clients should be careful not to add extra data that has deletion
 * dependancies (instead consider using nested RefCountPtr objects as extra
 * data which will guarantee the order of deletion).
 *
 * Preconditions:<ul>
 * <li> <tt>p->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> If this function has already been called with the same template
 *      type <tt>T1</tt> for <tt>extra_data</tt> and the same string <tt>name</tt>
 *      and <tt>force_unique==true</tt>, then an <tt>std::invalid_argument</tt>
 *      exception will be thrown.
 * </ul>
 *
 * Note, this function is made a non-member function to be consistent
 * with the non-member <tt>get_extra_data()</tt> functions.
 */
template<class T1, class T2>
void set_extra_data( const T1 &extra_data, const std::string& name, RefCountPtr<T2> *p, bool force_unique
#ifndef __sun
                     = true
#endif
                     ,EPrePostDestruction destroy_when
#ifndef __sun
                     = POST_DESTROY
#endif
	);
#ifdef __sun
template<class T1, class T2>
inline void set_extra_data( const T1 &extra_data, const std::string& name, RefCountPtr<T2> *p )
{ set_extra_data( extra_data, name, p, true, POST_DESTROY ); }
template<class T1, class T2>
inline void set_extra_data( const T1 &extra_data, const std::string& name, RefCountPtr<T2> *p, bool force_unique )
{ set_extra_data( extra_data, name, p, force_unique, POST_DESTROY ); }
#endif

/** \brief Get a non-const reference to extra data associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  p    [in] Smart pointer object that extra data is being extraced from.
 * @param  name [in] Name of the extra data.
 *
 * @return Returns a non-const reference to the extra_data object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> <tt>name</tt> and <tt>T1</tt> must have been used in a previous
 *      call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>).
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 */
template<class T1, class T2>
T1& get_extra_data( RefCountPtr<T2>& p, const std::string& name );

/** \brief Get a const reference to extra data associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  p    [in] Smart pointer object that extra data is being extraced from.
 * @param  name [in] Name of the extra data.
 *
 * @return Returns a const reference to the extra_data object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> <tt>name</tt> and <tt>T1</tt> must have been used in a previous
 *      call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>).
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * Also note that this const version is a false sense of security
 * since a client can always copy a const <tt>RefCountPtr</tt> object
 * into a non-const object and then use the non-const version to
 * change the data.  However, its presence will help to avoid some
 * types of accidental changes to this extra data.
 */
template<class T1, class T2>
const T1& get_extra_data( const RefCountPtr<T2>& p, const std::string& name );

/** \brief Get a non-const pointer to extra data (if it exists) associated
 * with a <tt>RefCountPtr</tt> object.
 *
 * @param  p    [in] Smart pointer object that extra data is being extraced from.
 * @param  name [in] Name of the extra data.
 *
 * @return Returns a non-const pointer to the extra_data object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> If <tt>name</tt> and <tt>T1</tt> have been used in a previous
 *      call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>)
        then <tt>return !=NULL</tt> and otherwise <tt>return == NULL</tt>.
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 */
template<class T1, class T2>
T1* get_optional_extra_data( RefCountPtr<T2>& p, const std::string& name );

/** \brief Get a const pointer to extra data (if it exists) associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  p    [in] Smart pointer object that extra data is being extraced from.
 * @param  name [in] Name of the extra data.
 *
 * @return Returns a const pointer to the extra_data object if it exists.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> If <tt>name</tt> and <tt>T1</tt> have been used in a previous
 *      call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>)
        then <tt>return !=NULL</tt> and otherwise <tt>return == NULL</tt>.
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * Also note that this const version is a false sense of security
 * since a client can always copy a const <tt>RefCountPtr</tt> object
 * into a non-const object and then use the non-const version to
 * change the data.  However, its presence will help to avoid some
 * types of accidental changes to this extra data.
 */
template<class T1, class T2>
const T1* get_optional_extra_data( const RefCountPtr<T2>& p, const std::string& name );

/** \brief Return a non-<tt>const</tt> reference to the underlying deallocator object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> The deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 *      (throws <tt>std::logic_error</tt>)
 * </ul>
 *
 */
template<class Dealloc_T, class T>
Dealloc_T& get_dealloc( RefCountPtr<T>& p );

/** \brief Return a <tt>const</tt> reference to the underlying deallocator object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> The deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 *      (throws <tt>std::logic_error</tt>)
 * </ul>
 *
 * Note that the <tt>const</tt> version of this function provides only
 * a very ineffective attempt to avoid accidental changes to the
 * deallocation object.  A client can always just create a new
 * non-<tt>const</tt> <tt>RefCountPtr<T></tt> object from any
 * <tt>const</tt> <tt>RefCountPtr<T></tt> object and then call the
 * non-<tt>const</tt> version of this function.
 */
template<class Dealloc_T, class T>
const Dealloc_T& get_dealloc( const RefCountPtr<T>& p );

//@}

} // end namespace Teuchos

#endif	// TEUCHOS_REFCOUNTPTR_DECL_H
