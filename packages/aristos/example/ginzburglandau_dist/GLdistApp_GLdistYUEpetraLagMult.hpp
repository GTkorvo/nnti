#ifndef GLDISTAPP_YUEPETRALAGMULT_H
#define GLDISTAPP_YUEPETRALAGMULT_H

#include "Aristos_YUEpetraVector.hpp"
#include "Aristos_LagMult.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"

namespace GLdistApp {

/**
    \brief Implements the Aristos::LagMult interface module for the parallel
    Ginzburg-Landau (GLdist) application.
*/
class GLdistYUEpetraLagMult : public Aristos::LagMult {
private:

  Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat_;

public:

  GLdistYUEpetraLagMult( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat );

  void getValue( const Aristos::Vector &x, const Aristos::Vector &g,
      Aristos::Vector &l, double &tol) const;

};

} // namespace GLdistApp

#endif
