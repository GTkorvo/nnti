// ///////////////////////////////
// TestRefCountPtrMain.cpp
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

#include <assert.h>

#include <iostream>

#include "Teuchos_RefCountPtr.hpp"

// Return constants from class functions
const int
    A_g_return  = 1,
	A_f_return  = 2,
	B1_g_return = 3,
	B1_f_return = 4,
	B2_g_return = 5,
	B2_f_return = 6,
	C_g_return  = 7,
	C_f_return  = 8,
	D_g_return  = 9,
	D_f_return  = 10,
	E_g_return  = 11,
	E_f_return  = 12;

// ////////////////////////////////////////////////
//
// Polymorphic multiple inheritance example
//
//            -----
//           |  A  |
//            -----
//             /|\
//              | 
//         ------------
//        |            |
//      -----        ------
//     |  B1 |      |  B2  |
//      -----        ------
//       /|\          /|\
//        |            |
//         ------------
//              |
//            -----
//           |  C  |
//            -----
//

class A
{
	int A_g_, A_f_;
public:
	A() : A_g_(A_g_return), A_f_(A_f_return) {}
	virtual ~A() {}
	virtual int A_g() { return A_g_; }
	virtual int A_f() const { return A_f_; }
};

class B1 : virtual public A
{
	int B1_g_, B1_f_;
public:
	B1() : B1_g_(B1_g_return), B1_f_(B1_f_return) {}
	virtual int B1_g() { return B1_g_; }
	virtual int B1_f() const { return B1_f_; }
};

class B2 : virtual public A
{
	int B2_g_, B2_f_;
public:
	B2() : B2_g_(B2_g_return), B2_f_(B2_f_return) {}
	virtual int B2_g() { return B2_g_; }
	virtual int B2_f() const { return B2_f_; }
};

class C : virtual public B1, virtual public B2
{
	int C_g_, C_f_;
public:
	C() : C_g_(C_g_return), C_f_(C_f_return) {}
	virtual int C_g() { return C_g_; }
	virtual int C_f() const { return C_f_; }
	
};

// ////////////////////////////////////////////////
//
// Non-polymophic classes hiearchy examlpe
//
//            -----
//           |  D  |
//            -----
//             /|\
//              | 
//            -----
//           |  E  |
//            -----
//

class D 
{
	int D_g_, D_f_;
public:
	D() : D_g_(D_g_return), D_f_(D_f_return) {}
	int D_g() { return D_g_; }
	int D_f() const { return D_f_; }
};

class E : public D
{
	int E_g_, E_f_;
public:
	E() : E_g_(E_g_return), E_f_(E_f_return) {}
	int E_g() { return E_g_; }
	int E_f() const { return E_f_; }
};

// //////////////////////////

// 
// Uncomment these macros to see example errors
//

//#define SHOW_COMPILE_TIME_ERRORS
//#define SHOW_RUN_TIME_ERROR_1
//#define SHOW_RUN_TIME_ERROR_2
//#define SHOW_RUN_TIME_ERROR_3
//#define SHOW_RUN_TIME_ERROR_4
//#define SHOW_MEMORY_LEAK_1

//
// This program prints minimal output to standard error
//
int main() {

	namespace rcp = MemMngPack;
	using rcp::RefCountPtr;
	using rcp::rcp_implicit_cast;
	using rcp::rcp_const_cast;
	using rcp::rcp_static_cast;
	using rcp::rcp_dynamic_cast;
	
	try {

	// Create some smart pointers

	RefCountPtr<A>       a_ptr1  = rcp_implicit_cast<A>(rcp::rcp(new C));
	assert( a_ptr1.get() );
	assert( a_ptr1.count()  == 1 );
#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
	RefCountPtr<D>       d_ptr1  = rcp::rcp(new E);
#else
	RefCountPtr<D>       d_ptr1  = rcp_implicit_cast<D>(rcp::rcp(new E));
#endif
	assert( d_ptr1.get() );
	assert( d_ptr1.count()  == 1 );

	{

		// Create some more smart points (no new memory!)

		const RefCountPtr<const A> ca_ptr1 = rcp::rcp_const_cast<const A>(a_ptr1); 
		assert( a_ptr1.count()  == 2 );
		assert( ca_ptr1.get() );
		assert( ca_ptr1.count() == 2 );
		const RefCountPtr<const D> cd_ptr1 = rcp::rcp_const_cast<const D>(d_ptr1);
		assert( d_ptr1.count()  == 2 );
		assert( cd_ptr1.get() );
		assert( cd_ptr1.count() == 2 );

#ifdef SHOW_RUN_TIME_ERROR_1
		// Conversion using get() is a no no!  When a_ptr2 is deleted so will the allocated
		// object and then a_ptr1 will be corrupted after this block ends!
		const RefCountPtr<A> a_ptr2 = a_ptr1.get();
#endif

		// Test assignment functions

		a_ptr1 = rcp::rcp_const_cast<A>(ca_ptr1); // Should be okay, assignment to self

#ifdef SHOW_COMPILE_TIME_ERRORS
		ca_ptr1 = ca_ptr1; // Should not compile since ca_ptr1 is declared constant
		ca_ptr1->A_g();    // Should not compile since A_g() is a non-const member function
#endif

		// Test function calls through operaor->(...)

		assert( a_ptr1->A_g()  == A_g_return );
		assert( a_ptr1->A_f()  == A_f_return );
		assert( ca_ptr1->A_f() == A_f_return );
		assert( d_ptr1->D_g()  == D_g_return );
		assert( d_ptr1->D_f()  == D_f_return );
		assert( cd_ptr1->D_f() == D_f_return );
		
		// Test funciton calls through operator*(...)

		assert( (*a_ptr1).A_g()  == A_g_return );
		assert( (*a_ptr1).A_f()  == A_f_return );
		assert( (*ca_ptr1).A_f() == A_f_return );
		assert( (*d_ptr1).D_g()  == D_g_return );
		assert( (*d_ptr1).D_f()  == D_f_return );
		assert( (*cd_ptr1).D_f() == D_f_return );

		// Test dynamic and static conversions

		// Cast down the inheritance hiearchy (const A -> const B1)
		const RefCountPtr<const B1> cb1_ptr1 = rcp::rcp_dynamic_cast<const B1>(ca_ptr1);
		assert( cb1_ptr1.get() );
		assert( cb1_ptr1.count() == 3 );
		assert( ca_ptr1.count()  == 3 );
		assert( a_ptr1.count()   == 3 );

		// Cast up the inheritance hiearchy (const B1 -> const A)
		assert( rcp::rcp_implicit_cast<const A>(cb1_ptr1)->A_f()  == A_f_return );
#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
		assert( RefCountPtr<const A>(cb1_ptr1)->A_f()           == A_f_return );
#endif
		// Implicit cast from const to non-const (A -> const A)
		assert( rcp::rcp_implicit_cast<const A>(a_ptr1)->A_f()    == A_f_return );
#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
		assert( RefCountPtr<const A>(a_ptr1)->A_f()             == A_f_return );
#endif
		// Cast away constantness (const B1 -> B1)
		assert( rcp::rcp_const_cast<B1>(cb1_ptr1)->B1_g()         == B1_g_return );
		// Cast across the inheritance hiearchy (const B1 -> const B2)
		assert( rcp::rcp_dynamic_cast<const B2>(cb1_ptr1)->B2_f() == B2_f_return );
		// Cast down the inheritance hiearchy (const B1 -> const C)
		assert( rcp::rcp_dynamic_cast<const C>(cb1_ptr1)->C_f()   == C_f_return );

		// Cast away constantness (const C -> C)
		const RefCountPtr<C>
			c_ptr1 = rcp::rcp_const_cast<C>(rcp::rcp_dynamic_cast<const C>(ca_ptr1));
		assert( c_ptr1.get() );
		assert( c_ptr1.count()   == 4 );
		assert( ca_ptr1.count()  == 4 );
		assert( a_ptr1.count()   == 4 );

		// Cast down the inheritance hiearchy using static_cast<...> (const D -> const E)
		const RefCountPtr<const E>
			ce_ptr1 = rcp::rcp_static_cast<const E>(cd_ptr1); // This is not checked at runtime!
		assert( ce_ptr1.get() );
		assert( ce_ptr1.count()  == 3 );
		assert( cd_ptr1.count()  == 3 );
		assert( d_ptr1.count()   == 3 );

		// Cast up the inheritance hiearchy (const E -> const D)
		assert( rcp::rcp_implicit_cast<const D>(ce_ptr1)->D_f()   == D_f_return ); 
		// Cast away constantness (const E -> E)
		assert( rcp::rcp_const_cast<E>(ce_ptr1)->E_g()            == E_g_return );
		assert( ce_ptr1->D_f()                                    == D_f_return );

#ifdef SHOW_COMPILE_TIME_ERRORS
		// Try to cast down inheritance hiearchy using dynamic_cast<...> (const D -> const E)
		rcp::rcp_dynamic_cast<const E>( cd_ptr1 )->E_f();  // This should not compile since D and E are not polymophic
#endif

		try {
			// Try to cast form one interface to another that is not supported (B2 -> B1).
			// The RefCountPtr<B1> returned from rcp_dynamic_cast<...> should be null!
			// Note that RefCountPtr<...>::optertor->() should throw an exception but even
			// so no memory leak occurs.  If you don't believe me then step through with a
			// debugger and see for yourself.
			assert( rcp::rcp_dynamic_cast<B1>( rcp::rcp(new B2) )->B1_g() == B1_g_return );
			return -1; // Should not be executed!
		}
		catch( const std::logic_error )
		{}

		// Manually clean up some memory

		delete d_ptr1.release();  // Now d_ptr1.get() no longer points to a valid object but okay
		                          // as long as no other access to this object is attempted! (see below)
#ifdef SHOW_RUN_TIME_ERROR_2
		assert( d_ptr1->D_g() == D_g_return ); // Should cause a segmentation fault since d_ptr.get() was deleted!
#endif

#ifdef SHOW_MEMORY_LEAK_1
		a_ptr1.release(); // If we release but do not delete manually then this is a memory leak!
#endif

		// Here at the end of the block, all of the other smart pointers are deleted!
	}
	// Check that all of the other references where removed but these
	assert( a_ptr1.count() == 1 );
	assert( d_ptr1.count() == 1 );

	// Assign some other dynamically created objects.
	
	a_ptr1 = rcp::rcp(new A);  // In each case the current dynamically allocated object is deleted ...
	a_ptr1 = rcp_implicit_cast<A>(rcp::rcp(new B1)); // before the new reference is set.
	a_ptr1 = rcp_implicit_cast<A>(rcp::rcp(new B2)); // ""
	a_ptr1 = rcp_implicit_cast<A>(rcp::rcp(new C));  // ""
	d_ptr1 = rcp::rcp(new D);                        // ""
#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
	d_ptr1 = rcp::rcp(new E);                        // ""
#else
	d_ptr1 = rcp_implicit_cast<D>(rcp::rcp(new E));
#endif

	// Assign pointers to some automatic objects that do not need deleted.
	// We can do this but we need to remove ownership of the pointer
	// from the smart pointer objects so that they do not try to
	// delete them.  If we forget then delete will be called on these
	// pointers and will cause a runtime error.

	C c; // Automatic object what will be deleted by compiler at end of block
#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
	a_ptr1 = rcp::rcp(&c);
#else
	a_ptr1 = rcp_implicit_cast<A>(rcp::rcp(&c));
#endif
#ifndef SHOW_RUN_TIME_ERROR_3
	// Release ownership so that a_ptr1 will not try to delete &c when a_ptr1 goes out of scope
	a_ptr1.release();
#endif	

	E e; // Automatic object what will be deleted by compiler at end of block
#ifdef REF_COUNT_PTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
	d_ptr1 = rcp::rcp(&e);
#else
	d_ptr1 = rcp_implicit_cast<D>(rcp::rcp(&e));
#endif
#ifndef SHOW_RUN_TIME_ERROR_4
	// Release ownership so that d_ptr1 will not try to delete &e when a_ptr1 goes out of scope
	d_ptr1.release();
#endif	

	// Set pointers to null to force releasing any owned memory
	a_ptr1 = rcp::null;
	d_ptr1 = rcp::null;

	std::cerr << "RefCountPtr<...> seems to checks out!\n";

	} // end try
	catch( const std::exception &excpt ) {
		std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
		return -1;
	}

	return 0;

}
