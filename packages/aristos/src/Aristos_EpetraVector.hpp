#ifndef ARISTOS_EPETRAVECTOR_H
#define ARISTOS_EPETRAVECTOR_H

#include "Aristos_Vector.hpp"
#include "Epetra_Vector.h"


/** \class Aristos::EpetraVector
    \brief The Aristos::Vector / Epetra_Vector adapter class.
      
    Holds a pointer to an Epetra_Vector and implements the member functions of the Aristos::Vector class.
*/


namespace Aristos {

class EpetraVector : public Vector {

private:
  
  Teuchos::RefCountPtr<Epetra_Vector>  epetra_vec_;

public:

  EpetraVector( const Teuchos::RefCountPtr<Epetra_Vector> &epetra_vec );

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
  Teuchos::RefCountPtr<const Epetra_Vector> getVector() const;


}; // class EpetraVector

} // namespace Aristos

#endif
