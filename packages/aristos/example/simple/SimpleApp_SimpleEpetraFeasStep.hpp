#ifndef SIMPLEAPP_EPETRAFEASSTEP_H
#define SIMPLEAPP_EPETRAFEASSTEP_H

#include "Aristos_EpetraVector.hpp"
#include "Aristos_FeasStep.hpp"
#include "SimpleApp_SimpleEpetraDataPool.hpp"

namespace SimpleApp {

class SimpleEpetraFeasStep : public Aristos::FeasStep {
private:

  Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat_;
  
public:

  SimpleEpetraFeasStep( Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat );
  
  void getValue( const Aristos::Vector &x, const Aristos::Vector &c, Aristos::Vector &fs, double &tol) const;

};

} // namespace SimpleApp

#endif
