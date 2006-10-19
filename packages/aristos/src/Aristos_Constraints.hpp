#ifndef ARISTOS_CONSTRAINTS_H
#define ARISTOS_CONSTRAINTS_H

#include "Aristos_Vector.hpp"

/** \class Aristos::Constraints
    \brief Provides interfaces for the evaluation of quantities related to the constraint function.

    These include:\n
    \li constraint evaluation
    \li application of the Jacobian operator
    \li application of the constraint null space operator 
*/


namespace Aristos {

class Constraints {
public:

  virtual ~Constraints() {}
  
  /** \brief Evaluates the vector of constraints.

      \param x [in]  - Current SQP iterate vector.
      \param c [out] - Value of constraints.

      \return None.

      \par Detailed Description:
      
      Interface function that evaluates the vector of constraints. To
      be subclassed and implemented by the user. 

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void getValue( const Vector &x, Vector &c ) const = 0;



  /** \brief Applies the Jacobian operator to vector.
  
      \param Trans [in] - If Trans==<tt>true</tt>, apply adjoint operator, otherwise apply standard operator.
      \param x [in]     - Current SQP iterate vector.
      \param v [in]     - The vector to which the Jacobian operator is applied.
      \param Jv [out]   - Resulting vector.

      \return None.

      \par Detailed Description:

      Interface function that evaluates the vector obtained by applying the Jacobian operator
      (or the adjoint thereof) to a vector. To be subclassed and implemented by the user. 

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void applyJacobian( bool Trans, const Vector &x, const Vector &v, Vector &Jv) const = 0;



  /** \brief Applies the constraint null space operator to a vector.
  
      \param Trans [in] - If Trans==<tt>true</tt>, apply adjoint operator, otherwise apply standard operator.
      \param x [in]     - Current SQP iterate vector.
      \param v [in]     - The vector to which the constraint null space operator is applied.
      \param Wv [out]   - Resulting vector.
      \param tol [in]   - Tolerance for inexact computations.

      \return None.

      \par Detailed Description:

      Interface function that evaluates the vector obtained by applying the constraint null space operator
      (or the adjoint thereof) to a vector. To be subclassed and implemented by the user. 

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void applyNullSp( bool Trans, const Vector &x, const Vector &v, Vector &Wv, double * tol) const = 0;

}; // class Constraints

} // namespace Aristos

#endif
