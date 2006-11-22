#ifndef GLDISTAPP_YUEPETRAFEASSTEP_H
#define GLDISTAPP_YUEPETRAFEASSTEP_H

#include "Aristos_YUEpetraVector.hpp"
#include "Aristos_FeasStep.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"

namespace GLdistApp {

/**
    \brief Implements the Aristos::FeasStep interface module for the parallel
    Ginzburg-Landau (GLdist) application.
*/
class GLdistYUEpetraFeasStep : 
public Aristos::FeasStep {
private:

  Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat_;

public:

  GLdistYUEpetraFeasStep( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat );

  void getValue( const Aristos::Vector &x, const Aristos::Vector &c,
      Aristos::Vector &fs, double &tol) const;

};

} // namespace GLdistApp

#endif
