#ifndef ARISTOS_HESSVEC_H
#define ARISTOS_HESSVEC_H

#include "Aristos_Vector.hpp"

/** \class Aristos::HessVec
    \brief Provides the interface for the evaluation of Hessian-of-Lagrangian times vector products.
*/


namespace Aristos {

class HessVec {
public:

  virtual ~HessVec() {}
  
  /** \brief Returns Hessian-of-Lagrangian times vector product.

      \param x [in]     - Current SQP iterate vector.
      \param l [in]     - The vector of Lagrange multipliers.
      \param s [in]     - The vector to which the Hessian-of-Lagrangian operator is applied.
      \param Hs [out]   - Resulting vector.

      \return None.

      \par Detailed Description:

      Interface function that evaluates the vector obtained by applying the Hessian-of-Lagrangian operator
      to a vector. To be subclassed and implemented by the user.

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void getValue( const Vector &x, const Vector &l, const Vector &s, Vector &Hs) const = 0;

}; // class HessVec

} // namespace Aristos

#endif
