#ifndef SIMPLEAPP_EPETRAOBJECTIVE_H
#define SIMPLEAPP_EPETRAOBJECTIVE_H

#include "Aristos_Objective.hpp"
#include "Aristos_EpetraVector.hpp"

namespace SimpleApp {

// Does inner product
class SimpleEpetraObjective : public Aristos::Objective {
private:

public:

  SimpleEpetraObjective();
  
  double getValue( const Aristos::Vector &x ) const;

  void getGradient( const Aristos::Vector &x, Aristos::Vector &g ) const;

};

} // namespace SimpleApp

#endif
