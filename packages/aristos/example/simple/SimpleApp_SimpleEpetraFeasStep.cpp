#include <SimpleApp_SimpleEpetraFeasStep.hpp>
#include "Epetra_Map.h"

namespace SimpleApp {

SimpleEpetraFeasStep::SimpleEpetraFeasStep( Teuchos::RefCountPtr<SimpleEpetraDataPool> dat )
  : dat_(dat)
{}


void SimpleEpetraFeasStep::getValue( const Aristos::Vector &x, const Aristos::Vector &c,
                                 Aristos::Vector &fs, double &tol) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> ec =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(c))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> efs =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(fs)).getVector());

  // Now do arithmetic in Epetra-land.

  Epetra_Vector rhs((dat_->getAugmat())->RowMap());
  
  for (int i=0; i < 5; i++)
    rhs[i] = 0.0;
  for (int i=5; i < 8; i++)
    rhs[i] = (*ec)[i-5];

  Epetra_Vector result((dat_->getAugmat())->RowMap());

  dat_->solveAugsys(rhs, result);

  Epetra_Vector tmp(ec->Map(), false);

  for (int i=0; i < 3; i++)
    tmp[i] = result[i+5];

  Epetra_CrsMatrix* myjac = dat_->getJacobian();

  myjac->Multiply(true, tmp, *efs);

}

} // namespace SimpleApp
