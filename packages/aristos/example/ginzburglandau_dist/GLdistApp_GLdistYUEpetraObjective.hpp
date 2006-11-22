#ifndef GLDISTAPP_YUEPETRAOBJECTIVE_H
#define GLDISTAPP_YUEPETRAOBJECTIVE_H

#include "Aristos_Objective.hpp"
#include "Aristos_YUEpetraVector.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"

namespace GLdistApp {

/**
    \brief Implements the Aristos::Objective interface module for the parallel
    Ginzburg-Landau (GLdist) application.
*/
class GLdistYUEpetraObjective : public Aristos::Objective {
private:

  Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat_;

public:

  GLdistYUEpetraObjective( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool>  dat );

  double getValue( const Aristos::Vector &x ) const;

  void getGradient( const Aristos::Vector &x, Aristos::Vector &g ) const;

};

} // namespace GLdistApp

#endif
