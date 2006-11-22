#ifndef GLDISTAPP_SCHUROP_H
#define GLDISTAPP_SCHUROP_H

#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Object.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RefCountPtr.hpp"

namespace GLdistApp {

class SchurOp : public Epetra_Operator, public Epetra_Object {
private:
  // Volume stiffness matrix.
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> A_;
  // Control/state mass matrix.
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> B_;
  // Jacobian of the nonlinear term.
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> N_;
  
public:
  SchurOp(Teuchos::RefCountPtr<Epetra_FECrsMatrix> A,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> B,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> N)
    : A_(A), B_(B), N_(N) {};
             
  ~SchurOp() {};

  inline int SetUseTranspose(bool UseTranspose) { return -1; }
  
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { return -1; }

  double NormInf() const { return 0.0; }

  // Returns a character string describing the operator
  const char* Label() const {return(Epetra_Object::Label());}

  bool UseTranspose() const {return(false);}

  // Returns false because this class can't compute an Inf-norm.
  bool HasNormInf() const {return(false);}

  const Epetra_Comm& Comm() const {return A_->Comm();}

  const Epetra_Map& OperatorDomainMap() const {return (Epetra_Map &) A_->DomainMap();}

  const Epetra_Map& OperatorRangeMap() const {return (Epetra_Map &) A_->RangeMap();}

};

}  // End namespace GLdistApp.

#endif
