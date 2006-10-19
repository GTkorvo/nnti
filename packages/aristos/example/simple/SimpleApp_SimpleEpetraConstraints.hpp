#ifndef SIMPLEAPP_EPETRACONSTRAINTS_H
#define SIMPLEAPP_EPETRACONSTRAINTS_H

#include "Aristos_EpetraVector.hpp"
#include "Aristos_Constraints.hpp"
#include "SimpleApp_SimpleEpetraDataPool.hpp"

namespace SimpleApp {

class SimpleEpetraConstraints : public Aristos::Constraints {
private:

  Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat_;
  
public:

  SimpleEpetraConstraints( Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat );
  
  void getValue( const Aristos::Vector &x, Aristos::Vector &c ) const;

  void applyJacobian( bool Trans, const Aristos::Vector &x, const Aristos::Vector &v, Aristos::Vector &Jv) const;

  void applyNullSp( bool Trans, const Aristos::Vector &x, const Aristos::Vector &v, Aristos::Vector &Wv,
                    double * tol) const;

};

} // namespace SimpleApp

#endif
