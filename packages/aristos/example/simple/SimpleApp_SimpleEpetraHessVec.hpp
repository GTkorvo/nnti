#ifndef SIMPLEAPP_EPETRAHESSVEC_H
#define SIMPLEAPP_EPETRAHESSVEC_H

#include "Aristos_EpetraVector.hpp"
#include "Aristos_HessVec.hpp"
#include "SimpleApp_SimpleEpetraDataPool.hpp"

namespace SimpleApp {

class SimpleEpetraHessVec : public Aristos::HessVec {
private:

public:

  SimpleEpetraHessVec();
  
  void getValue( const Aristos::Vector &x, const Aristos::Vector &l, const Aristos::Vector &s,
                 Aristos::Vector &Hs) const;

};

} // namespace SimpleApp

#endif
