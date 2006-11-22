#include <GLdistApp_GLdistYUEpetraLagMult.hpp>

namespace GLdistApp {

GLdistYUEpetraLagMult::GLdistYUEpetraLagMult(
    Teuchos::RefCountPtr<GLdistYUEpetraDataPool> dat )
  :dat_(dat)
{}


void GLdistYUEpetraLagMult::getValue( const Aristos::Vector &x,
    const Aristos::Vector &g, Aristos::Vector &l, double &tol) const
{

  // Dynamic cast back to Epetra_MultiVectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> egy =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(g))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> egu =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(g))).getUVector();
  Teuchos::RefCountPtr<Epetra_MultiVector> el =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(l)).getYVector());
  
  Teuchos::RefCountPtr<Epetra_MultiVector> egp = Teuchos::rcp(new Epetra_MultiVector(egy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> ey = Teuchos::rcp(new Epetra_MultiVector(egy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> eu = Teuchos::rcp(new Epetra_MultiVector(egu->Map(),1));

  dat_->solveAugsysDyn( egy, egu, egp, ey, eu, el, &tol );

}

} // namespace GLdistApp
