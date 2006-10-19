#include <SimpleApp_SimpleEpetraLagMult.hpp>
#include "Epetra_Map.h"

namespace SimpleApp {

SimpleEpetraLagMult::SimpleEpetraLagMult( Teuchos::RefCountPtr<SimpleEpetraDataPool> dat )
  : dat_(dat)
{}


void SimpleEpetraLagMult::getValue( const Aristos::Vector &x, const Aristos::Vector &g,
                                 Aristos::Vector &l, double &tol) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> eg =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(g))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> el =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(l)).getVector());

  // Now do arithmetic in Epetra-land.

  Epetra_Vector rhs((dat_->getAugmat())->RowMap());
  
  for (int i=0; i < 5; i++)
    rhs[i] = (*eg)[i];
  for (int i=5; i < 8; i++)
    rhs[i] = 0.0;

  Epetra_Vector result((dat_->getAugmat())->RowMap());

  dat_->solveAugsys(rhs, result);

  for (int i=0; i < 3; i++)
    (*el)[i] = result[i+5];

}

} // namespace SimpleApp
