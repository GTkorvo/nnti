#include "Aristos_SledgeVector.hpp"

namespace Aristos {

SledgeVector::SledgeVector( const Teuchos::RefCountPtr<USEME::umDVec> &sledge_vec )
  :sledge_vec_(sledge_vec)
    {}

// Overridden from Vector

double SledgeVector::innerProd( const Vector &x ) const
{
  SledgeVector &ex = Teuchos::dyn_cast<SledgeVector>(const_cast <Vector&>(x));
  return sledge_vec_->d_dot(*ex.sledge_vec_);
}

void SledgeVector::linComb( const double &alpha, const Vector &x, const double &beta )
{
  SledgeVector &ex = Teuchos::dyn_cast<SledgeVector>(const_cast <Vector&>(x));

  (*sledge_vec_) = alpha*(*ex.sledge_vec_) + beta*(*sledge_vec_);
}

void SledgeVector::Scale( const double &alpha )
{
  (*sledge_vec_) *= alpha;
}

void SledgeVector::Set( const double &alpha )
{
  (*sledge_vec_) = alpha;
}

void SledgeVector::Set( const double &alpha, const Vector &x )
{
  SledgeVector &ex = Teuchos::dyn_cast<SledgeVector>(const_cast <Vector&>(x));
  (*sledge_vec_) = alpha * (*ex.sledge_vec_);
}

Teuchos::RefCountPtr<Vector> SledgeVector::createVector() const
{
  return Teuchos::rcp( new SledgeVector( Teuchos::rcp(new USEME::umDVec(sledge_vec_->size(), "vec")) ));
}

Teuchos::RefCountPtr<const USEME::umDVec> SledgeVector::getVector() const
{
  return sledge_vec_;
}

} // namespace Aristos
