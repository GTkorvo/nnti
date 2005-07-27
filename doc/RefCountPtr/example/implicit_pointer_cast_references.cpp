class A {};
class C : public A{};

// Typedefs
typedef A*        ptr_A;
typedef const A*  ptr_const_A;

// Takes pointer to A (by value)
void foo1( ptr_A a_ptr );

// Takes pointer to const A (by value)
void foo2( ptr_const_A a_ptr );

// Takes pointer to A (by const reference) 
void foo3( const ptr_A &a_ptr );

// Takes pointer to A (by non-const reference) 
void foo4( ptr_A &a_ptr );

void boo1()
{
  C* c = new C;
  A* a = c;
  foo1(c); // Okay, implicit cast to base class
  foo2(a); // Okay, implicit cast to const
  foo2(c); // Okay, implicit cast to base class and const
  foo3(c); // Okay, implicit case to base class
  foo4(a); // Okay, no cast
  foo4(c); // Error, can not cast from (C&*) to (A&*)!
}
