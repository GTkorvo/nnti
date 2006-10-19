#ifndef ARISTOS_OBJECTIVE_H
#define ARISTOS_OBJECTIVE_H

#include "Aristos_Vector.hpp"

/** \class Aristos::Objective
    \brief Provides interfaces for the evaluation of quantities related to the objective function.

    These include:\n
    \li the value of the objective,
    \li the gradient of the objective.
*/


namespace Aristos {

class Objective {
public:

  virtual ~Objective() {}
  
  /** \brief Returns the value of the objective function.

      \param x [in]  - Current SQP iterate vector.

      \return value of objective.

      \par Detailed Description:

      Interface function that evaluates the objective function. To be subclassed and implemented by the user.

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual double getValue( const Vector &x ) const = 0;



  /** \brief Computes the gradient of the objective function.

      \param x [in]  - Current SQP iterate vector.
      \param g [out] - Gradient of the objective function.

      \return None.

      \par Detailed Description:

      Interface function that evaluates the gradient of the objective function. To
      be subclassed and implemented by the user.

      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax (this is only one example):\n
      <tt>
      Teuchos::RefCountPtr<const YourVectorClass> ex =
        (Teuchos::dyn_cast<Aristos::YourVectorAdapter>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void getGradient( const Vector &x, Vector &g ) const = 0;

}; // class Objective

} // namespace Aristos

#endif
