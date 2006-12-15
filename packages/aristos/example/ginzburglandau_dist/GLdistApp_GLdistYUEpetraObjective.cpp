#include "GLdistApp_GLdistYUEpetraObjective.hpp"

namespace GLdistApp {

GLdistYUEpetraObjective::GLdistYUEpetraObjective( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat )
  :dat_(dat)
{}

double GLdistYUEpetraObjective::getValue( const Aristos::Vector &x ) const
{
  Teuchos::RefCountPtr<const Epetra_MultiVector> ey =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> eu =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getUVector();

  double beta = dat_->getbeta();

  Teuchos::RefCountPtr<Epetra_FECrsMatrix> H = dat_->getH();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> R = dat_->getR();

  Teuchos::RefCountPtr<Epetra_FEVector> q = dat_->getq();

  Epetra_MultiVector yq(*ey);
  Epetra_MultiVector Hyq(*ey);
  Epetra_MultiVector Ru(*eu);
  
  yq.Update(-1.0, *q, 1.0);
  H->Multiply(false, yq, Hyq);
  R->Multiply(false, *eu, Ru);

  
  double doty[1];
  double dotu[1];
  
  yq.Dot(Hyq, doty);
  eu->Dot(Ru, dotu);

  return 0.5*doty[0] + 0.5*beta*dotu[0];
}

void GLdistYUEpetraObjective::getGradient( const Aristos::Vector &x, Aristos::Vector &g ) const
{
  // Dynamic cast back to YUEpetra vectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> ey =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> eu =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getUVector();

  Teuchos::RefCountPtr<Epetra_MultiVector> egy =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(g)).getYVector());
  Teuchos::RefCountPtr<Epetra_MultiVector> egu =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(g)).getUVector());

  double beta = dat_->getbeta();

  Teuchos::RefCountPtr<Epetra_FECrsMatrix> H = dat_->getH();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> R = dat_->getR();

  Teuchos::RefCountPtr<Epetra_FEVector> q = dat_->getq();

  Epetra_MultiVector yq(*ey);
  
  yq.Update(-1.0, *q, 1.0);

  H->Multiply(false, yq, *egy);
  R->Multiply(false, *eu, *egu);
  egu->Scale(beta);

}

} // namespace GLsApp
