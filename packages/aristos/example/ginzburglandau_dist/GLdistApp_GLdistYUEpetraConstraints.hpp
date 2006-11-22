#ifndef GLDISTAPP_YUEPETRACONSTRAINTS_H
#define GLDISTAPP_YUEPETRACONSTRAINTS_H

#include "Aristos_YUEpetraVector.hpp"
#include "Aristos_Constraints.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"

namespace GLdistApp {

/**
    \brief Implements the Aristos::Constraints interface module for the parallel
    Ginzburg-Landau (GLdist) application.
*/
class GLdistYUEpetraConstraints : public Aristos::Constraints {
private:

  Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat_;

public:

  GLdistYUEpetraConstraints( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool> dat);

  void getValue( const Aristos::Vector &x, Aristos::Vector &c ) const;

  void applyJacobian( bool Trans, const Aristos::Vector &x,
      const Aristos::Vector &v, Aristos::Vector &Jv) const;

  void applyNullSp( bool Trans, const Aristos::Vector &x,
      const Aristos::Vector &v, Aristos::Vector &Wv, double * tol) const;

};

} // namespace GLdistApp

#endif
