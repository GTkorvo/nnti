#include "example_get_args.hpp"

// Abstract interfaces
class UtilityBase {
public:
  virtual void f() const = 0;
};
class UtilityBaseFactory {
public:
  virtual UtilityBase* createUtility() const = 0;
};

// Concrete implementations
class UtilityA : public UtilityBase {
public:
  void f() const { std::cout<<"\nUtilityA::f() called, this="<<this<<"\n"; }
};
class UtilityB : public UtilityBase {
public:
  void f() const { std::cout<<"\nUtilityB::f() called, this="<<this<<"\n"; }
};
class UtilityAFactory : public UtilityBaseFactory {
public:
  UtilityBase* createUtility() const { return new UtilityA(); }
};
class UtilityBFactory : public UtilityBaseFactory {
public:
  UtilityBase* createUtility() const { return new UtilityB(); }
};

// Client classes
class ClientA {
public:
  void f( const UtilityBase &utility ) const { utility.f(); }
};
class ClientB {
  UtilityBase *utility_;
public:
  ClientB() : utility_(0) {}
  ~ClientB() { delete utility_; }
  void initialize( UtilityBase *utility ) { utility_ = utility; }
  void g( const ClientA &a ) { a.f(*utility_); }
};
class ClientC {
  const UtilityBaseFactory *utilityFactory_;
  UtilityBase              *utility_;
  bool                     shareUtility_;
public:
  ClientC( const UtilityBaseFactory *utilityFactory, bool shareUtility )
    :utilityFactory_(utilityFactory)
    ,utility_(utilityFactory->createUtility())
    ,shareUtility_(shareUtility) {}
  ~ClientC() { delete utilityFactory_; delete utility_; }
  void h( ClientB *b ) {
    if( shareUtility_ ) b->initialize(utility_);
    else                b->initialize(utilityFactory_->createUtility());
  }
};

// Main program
int main( int argc, char* argv[] )
{
  // Read options from the commandline
  bool useA, shareUtility;
  example_get_args(argc,argv,&useA,&shareUtility);
  // Create factory
  UtilityBaseFactory *utilityFactory = 0;
  if(useA) utilityFactory = new UtilityAFactory();
  else     utilityFactory = new UtilityBFactory();
  // Create clients
  ClientA a;
  ClientB b1, b2;
  ClientC c(utilityFactory,shareUtility);
  // Do some stuff
  c.h(&b1);
  c.h(&b2);
  b1.g(a);
  b2.g(a);
  // Cleanup memory
  delete utilityFactory;
}
