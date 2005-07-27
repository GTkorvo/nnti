#include <iostream>

class A {};
class C : public A{};

typedef A*        ptr_A;
typedef const A*  ptr_const_A;

ptr_A       foo9(const ptr_A       &a) { std::cout << "Called non-const A version!\n"; return a; }
ptr_const_A foo9(const ptr_const_A &a) { std::cout << "Called const A version!\n"; return a; }

ptr_A boo5() {
  C* c = new C;
  return foo9(c);                       // Error, call is ambiguous!
  ptr_A a = c;
  return foo9(a);                       // Okay, calls foo9(ptr_A)
}

int main()
{
  boo5();
  return 0;
}
