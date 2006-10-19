#ifndef SIMPLEAPP_EPETRALAGMULT_H
#define SIMPLEAPP_EPETRALAGMULT_H

#include "Aristos_EpetraVector.hpp"
#include "Aristos_LagMult.hpp"
#include "SimpleApp_SimpleEpetraDataPool.hpp"

namespace SimpleApp {

class SimpleEpetraLagMult : public Aristos::LagMult {
private:

  Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat_;
  
public:

  SimpleEpetraLagMult( Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat );
  
  void getValue( const Aristos::Vector &x, const Aristos::Vector &g, Aristos::Vector &l, double &tol) const;

};

} // namespace SimpleApp

#endif
