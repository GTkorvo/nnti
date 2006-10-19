#ifndef ARISTOS_SLEDGEVECTOR_H
#define ARISTOS_SLEDGEVECTOR_H

#include "Aristos_Vector.hpp"
#include "Sledge++.h"

/** \class Aristos::SledgeVector
    \brief The Aristos::Vector / USEME::umDVec (Sledge) adapter class.

    Holds a pointer to an USEME::umDVec and implements the member functions of the Aristos::Vector class.
*/


namespace Aristos {

class SledgeVector : public Vector {

private:

  Teuchos::RefCountPtr<USEME::umDVec>  sledge_vec_;

public:

  SledgeVector( const Teuchos::RefCountPtr<USEME::umDVec> &sledge_vec );

  /** \name Overridden from Vector */
  //@{

  double innerProd( const Vector &x ) const;

  void linComb( const double &alpha, const Vector &x, const double &beta );

  void Scale( const double &alpha );

  void Set( const double &alpha );

  void Set( const double &alpha, const Vector &x );

  Teuchos::RefCountPtr<Vector> createVector() const;

  //@}

  /** Returns a reference counted pointer to the private data container.
  */
  Teuchos::RefCountPtr<const USEME::umDVec> getVector() const;




}; // class SledgeVector

} // namespace Aristos

#endif
