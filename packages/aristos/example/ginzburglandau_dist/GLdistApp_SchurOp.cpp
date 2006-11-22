#include "GLdistApp_SchurOp.hpp"

namespace GLdistApp {

int SchurOp::Apply(const Epetra_MultiVector& x, Epetra_MultiVector& y) const
{
  Epetra_MultiVector statetmp1( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector statetmp2( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector nn( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector an( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector na( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector aa( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector ctrldomain( (Epetra_BlockMap&)B_->DomainMap(), 1);
  Epetra_MultiVector bb( (Epetra_BlockMap&)B_->RangeMap(), 1);

  N_->Multiply(true, x, statetmp1); // N'*x
  A_->Multiply(true, x, statetmp2); // A'*x

  N_->Multiply(false, statetmp1, nn); //N*N'*x
  A_->Multiply(false, statetmp1, an); //A*N'*x
  N_->Multiply(false, statetmp2, na); //N*A'*x
  A_->Multiply(false, statetmp2, aa); // A*A'*x
  
  B_->Multiply(true, x, ctrldomain); // B'*x
  B_->Multiply(false, ctrldomain, bb); // B*B'*x

  y.Update(1.0, nn, 1.0, an, 0.0);
  y.Update(1.0, na, 1.0, aa, 1.0);
  y.Update(1.0, bb, 1.0);

  return 0;
};

}  // End namespace GLdistApp.
