#include "Teuchos_RefCountPtr.hpp"
#include "example_get_args.hpp"

// Inject symbols for RefCountPtr so we don't need Teuchos:: qualification
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// Abstract interfaces
class UtilityBase {
public:
  virtual void f() const = 0;
};
class UtilityBaseFactory {
public:
  virtual RefCountPtr<UtilityBase> createUtility() const = 0;
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
  RefCountPtr<UtilityBase> createUtility() const { return rcp(new UtilityA()); }
};
class UtilityBFactory : public UtilityBaseFactory {
public:
  RefCountPtr<UtilityBase> createUtility() const { return rcp(new UtilityB()); }
};

// Client classes
class ClientA {
public:
  void f( const UtilityBase &utility ) const { utility.f(); }
};
class ClientB {
  RefCountPtr<UtilityBase> utility_;
public:
  void initialize(const RefCountPtr<UtilityBase> &utility) { utility_=utility; }
  void g( const ClientA &a ) { a.f(*utility_); }
};
class ClientC {
  RefCountPtr<UtilityBaseFactory> utilityFactory_;
  RefCountPtr<UtilityBase>        utility_;
  bool                            shareUtility_;
public:
  ClientC( const RefCountPtr<UtilityBaseFactory> &utilityFactory, bool shareUtility )
    :utilityFactory_(utilityFactory)
    ,utility_(utilityFactory->createUtility())
    ,shareUtility_(shareUtility) {}
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
  RefCountPtr<UtilityBaseFactory> utilityFactory;
  if(useA) utilityFactory = rcp(new UtilityAFactory());
  else     utilityFactory = rcp(new UtilityBFactory());
  // Create clients
  ClientA a;
  ClientB b1, b2;
  ClientC c(utilityFactory,shareUtility);
  // Do some stuff
  c.h(&b1);
  c.h(&b2);
  b1.g(a);
  b2.g(a);
}
