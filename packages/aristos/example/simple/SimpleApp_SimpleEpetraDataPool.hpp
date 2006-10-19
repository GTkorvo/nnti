#ifndef SIMPLEAPP_DATAPOOL_H
#define SIMPLEAPP_DATAPOOL_H

#include "Aristos_DataPool.hpp"
#include "Aristos_EpetraVector.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

namespace SimpleApp {

class SimpleEpetraDataPool : public Aristos::DataPool {
private:

  Epetra_CrsMatrix jacobian_;
  Epetra_CrsMatrix augmat_;
  Epetra_LinearProblem augsys_;
  
public:

  SimpleEpetraDataPool( Epetra_CrsMatrix &jacobian, Epetra_CrsMatrix &augmat );
  
  void computeAll( const Aristos::Vector &x );

  void solveAugsys( const Epetra_Vector &rhs, Epetra_Vector &result );

  Epetra_CrsMatrix* getJacobian();
  Epetra_CrsMatrix* getAugmat();

};

} // namespace SimpleApp

#endif
