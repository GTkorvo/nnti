#ifndef ARISTOS_DATAPOOL_H
#define ARISTOS_DATAPOOL_H

#include "Aristos_Vector.hpp"

/** \class Aristos::DataPool
    \brief Provides the interface to a generic data pool.

    The only member function is computeAll, responsible for the computation of all data that can be
    stored and reused within one (or more) SQP iteration. The data is entirely problem specific.
*/


namespace Aristos {

class DataPool {
public:

  virtual ~DataPool() {}
  
  /** \brief Recompute all stored quantities.
      \param x [in]  - Current SQP iterate vector.
    
      \return None. 
    
      \par Detailed Description:

      Interface function that evaluates and stores problem-specific quantities
      that can be reused within one (or more) SQP iteration.
  
      \note The Aristos::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax:\n
      <tt>
      Teuchos::RefCountPtr<const umDVec> ex =
        (Teuchos::dyn_cast<Aristos::SledgeVector>(const_cast<Aristos::Vector&>(x))).getVector();
      </tt>
  */
  virtual void computeAll( const Vector &x ) = 0;

}; // class DataPool

} // namespace Aristos

#endif
