#ifndef MLAPI_DOUBLEVECTOR_UTILS_H
#define MLAPI_DOUBLEVECTOR_UTILS_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "MLAPI_MultiVector.h"

namespace MLAPI {

/*!
\file MLAPI_MultiVector_Utils.h

\brief Utilities for MultiVector's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

//! Creates a new vector, x, such that x = y.
MultiVector Duplicate(const MultiVector& y)
{
  MultiVector x(y.GetVectorSpace());
  x.Update(y);
  return(x);
}

//! Extracts a given vector from input multivector.
MultiVector Extract(const MultiVector& y, const int v)
{
  if ((v < 0) || v > y.GetNumVectors())
    ML_THROW("Wrong input parameter v (" +
             GetString(v) + ")", -1);
      
  MultiVector x(y.GetVectorSpace(), 1);
  for (int i = 0 ; i < x.GetMyLength() ; ++i)
    x(i) = y(i,v);

  return(x);
}

} // namespace MLAPI

#endif // HAVE_ML_MLAPI

#endif
