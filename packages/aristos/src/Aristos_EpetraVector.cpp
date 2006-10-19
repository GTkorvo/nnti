#include "Aristos_EpetraVector.hpp"

namespace Aristos {

EpetraVector::EpetraVector( const Teuchos::RefCountPtr<Epetra_Vector> &epetra_vec )
  :epetra_vec_(epetra_vec)
    {}

// Overridden from Vector

double EpetraVector::innerProd( const Vector &x ) const
{
  double dot[1];
  EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(const_cast <Vector&>(x));
  //EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(x);
  epetra_vec_->Dot( *ex.epetra_vec_, dot );
  return dot[0];
}

void EpetraVector::linComb( const double &alpha, const Vector &x, const double &beta )
{
  EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(const_cast <Vector&>(x));
  //EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(x);
  epetra_vec_->Update( alpha, *ex.epetra_vec_, beta );
}

void EpetraVector::Scale( const double &alpha )
{
  epetra_vec_->Scale( alpha );
}

void EpetraVector::Set( const double &alpha )
{
  epetra_vec_->PutScalar( alpha );
}

void EpetraVector::Set( const double &alpha, const Vector &x )
{
  EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(const_cast <Vector&>(x));
  epetra_vec_->Scale( alpha, *ex.epetra_vec_ );
}

Teuchos::RefCountPtr<Vector> EpetraVector::createVector() const
{
  return Teuchos::rcp( new EpetraVector( Teuchos::rcp(new Epetra_Vector(epetra_vec_->Map(),false)) ));
}

Teuchos::RefCountPtr<const Epetra_Vector> EpetraVector::getVector() const
{
  return epetra_vec_;
}

} // namespace Aristos
