#ifndef ARISTOS_FEASSTEP_H
#define ARISTOS_FEASSTEP_H

#include "Aristos_Vector.hpp"

/** \class Aristos::FeasStep
    \brief Provides the interface for the computation of the feasibility step, used in the
           quasi-normal step computation.

    This can be a Newton step, solution of the state equations, etc. It is in general
    related to finding the minimum norm solution of:\n
    \f[
    \min \| c_x(x_k)s + c(x_k) \|
    \f]
    where \f$c_x(x_k)\f$ is the Jacobian operator and \f$c(x_k)\f$ is the value of constraints,
    at the current iterate \f$x_k\f$.
*/


namespace Aristos {

class FeasStep {
public:

  virtual ~FeasStep() {}
  
  /** \brief Returns feasibility step (part of the quasi-normal step computation).
    
      \param x [in]     - Current SQP iterate vector.
      \param c [in]     - Value of the constraints.
      \param fs [out]   - Feasibility step.
      \param tol [in]   - Tolerance for inexact computations.

      \return None.
      
      \par Detailed Description:

      Interface function that evaluates the feasibility step. 
      To be subclassed and implemented by the user.
        
      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void getValue( const Vector &x, const Vector &c, Vector &fs, double &tol) const = 0;

}; // class FeasStep

} // namespace Aristos

#endif
