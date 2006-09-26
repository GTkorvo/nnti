#ifndef RBGEN_UTILS_H
#define RBGEN_UTILS_H

class Epetra_MultiVector;

namespace RBGen {
  
  double BasisAngle( const Epetra_MultiVector& S, const Epetra_MultiVector& T );
  
} // end namespace RBGen

#endif // RBGEN_UTILS_H
