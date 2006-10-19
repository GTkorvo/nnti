#ifndef ARISTOS_LAGMULT_H
#define ARISTOS_LAGMULT_H

#include "Aristos_Vector.hpp"

/** \class Aristos::LagMult
    \brief Provides the interface for the evaluation of Lagrange multipliers.
*/


namespace Aristos {

class LagMult {
public:

  virtual ~LagMult() {}
  
  /** \brief Returns Lagrange multiplier estimate.

      \param x [in]     - Current SQP iterate vector.
      \param g [in]     - The gradient of the objective function.
      \param l [out]    - The vector of Lagrange multiplier estimates.
      \param tol [in]   - Tolerance for inexact computations.

      \return None.

      \par Detailed Description:

      Interface function that estimates the Lagrange multipliers. To be subclassed and implemented by the user.

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void getValue( const Vector &x, const Vector &g, Vector &l, double &tol) const = 0;

}; // class LagMult

} // namespace Aristos

#endif
