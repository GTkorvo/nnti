  //
  // Header-like declarations
  //
  
  #include <iostream>
  #include <cstdlib>
  
  namespace NamespaceA {
  
  template<class T>
  class A {
  public:
    explicit A( const T& a ) : a_(a) {}
    void print( std::ostream &os ) const { os << "\na="<<a_<<"\n"; }
  private:
    T a_;
  };
  
  } // namespace NamespaceA
  
  // Add a using declaration to inject 'A' into another namespace
  namespace NamespaceB { using NamespaceA::A; }
  
  // Now use the A class without namespace qualification in NamespaceB
  namespace NamespaceB {
  
  A<double> foo( std::ostream &os, const A<int> &aa );
  
  } // namespace NamespaceB
  
  // Create another A class in the global namespace.  With care, we should
  // not have any problems with this and our code should not be affected by
  // the presence of this class.
  template<class T>
  class A {
  public:
    explicit A( const T& a ) : a_(a) 
      { std::cerr << "\nOh no, called ::A::A(...)!\n"; std::exit(1); }
    void print( std::ostream &os ) { os << "\na="<<a_<<"\n"; }
  private:
    T a_;
  };
  
  // See what happens when you define another class A in NamespaceB which
  // conflicts with the using declaration!  This should not be allowed and
  // should be caught by the compiler!
  
  #ifdef SHOW_DUPLICATE_CLASS_A
  
  namespace NamespaceB {
  
  template<class T>
  class A {
  public:
    explicit A( const T& a ) : a_(a) 
      { std::cerr << "\nOh no, called ::A::A(...)!\n"; exit(1); }
    void print( std::ostream &os ) { os << "\na="<<a_<<"\n"; }
  private:
    T a_;
  };
  
  } // namespace NamespaceB
  
  #endif // SHOW_DUPLICATE_CLASS_A
 
  
  //
  // Implementations
  //
  
  // Define function in NamespaceB without namespace qualification for class A
  NamespaceB::A<double>
  NamespaceB::foo( std::ostream &os, const A<int> &aa )
  {
    A<double> ab(2.0);
    aa.print(std::cout);
    ab.print(std::cout);
    return ab;
  }
  // NOTE: Above, we need explicit namespace qualification for the return
  // type 'NamespaceB::A<double>' since we use namespace qualification to
  // define nonmember functions (see Thyra coding guidelines).  Without this
  // namespace qualification, the global class '::A' would be assumed and
  // you would get a compilation error.  However, within the function, which
  // is in the scope of NamespaceB, we don't need namespace qualifications!
  
  
  //
  // User's code.  This code does not typically live in a namespace (or is
  // in another unrelated namespace).  Here, some explicit namespace
  // qualification and using declarations will be required to avoid
  // ambiguities.
  //
  
  int main()
  {
  
  #if defined(SHOW_MISSING_USING_DECL)
    // Here, no using declaration is provided.  This will result in the
    // global class '::A' being used below which will result in a compiler
    // error when the NamespaceB::foo(...) function is called.  This is a
    // feature!
  #elif defined(SHOW_ERRONEOUS_USING_DIRECTIVE)
    // Here we try to just inject all of the names from NamespaceA into the
    // local scope.  However, this will result in the names 'NamespaceA::A'
    // and '::A' being equally visible which will result in a compiler error
    // when the first unqualified 'A' object is created below!
    using namespace NamespaceA;
  #else
    // Inject the class name 'A' into the local scope and will override any
    // (sloppy) names polluting the global namespace.  This will cause the
    // global '::A' class to be hidden (which is good!).
    using NamespaceA::A;
  #endif
  
    A<int> aa(5);
    A<double> ab = NamespaceB::foo(std::cout,aa);
    ab.print(std::cout);
    
    return 0;
    
  }
