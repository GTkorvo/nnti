#ifndef GLDISTAPP_YUEPETRAHESSVEC_H
#define GLDISTAPP_YUEPETRAHESSVEC_H

#include "Aristos_YUEpetraVector.hpp"
#include "Aristos_HessVec.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"

namespace GLdistApp {

/**
    \brief Implements the Aristos::HessVec interface module for the parallel
    Ginzburg-Landau (GLdist) application.
*/
class GLdistYUEpetraHessVec : public Aristos::HessVec {
private:
  Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool> dat_;
public:

  GLdistYUEpetraHessVec( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool> dat );

  void getValue( const Aristos::Vector &x, const Aristos::Vector &l,
      const Aristos::Vector &s, Aristos::Vector &Hs) const;

};

} // namespace GLdistApp

#endif
