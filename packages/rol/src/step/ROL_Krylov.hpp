// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_KRYLOV_H
#define ROL_KRYLOV_H

/** \class ROL::Krylov
    \brief Provides definitions for Krylov solvers.
*/

#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class Krylov {

  EKrylov ekv_;

  Real tol1_;
  Real tol2_;
  int  maxit_;
  bool useInexact_;

public:
  Krylov( EKrylov &ekv = KRYLOV_CG, Real tol1 = 1.e-4, Real tol2 = 1.e-2, 
          int maxit = 100, bool useInexact = false ) 
    : ekv_(ekv), tol1_(tol1), tol2_(tol2), maxit_(maxit), useInexact_(useInexact) {}

  // Run Krylov Method
  void run( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
            Objective<Real> &obj, Teuchos::RCP<Secant<Real> > secant = Teuchos::null ) {
    switch ( this->ekv_ ) {
      case KRYLOV_CG: this->CG(s,iter,flag,g,x,obj,secant); break;
      case KRYLOV_CR: this->CG(s,iter,flag,g,x,obj,secant); break;
      case KRYLOV_LAST: break; // DO NOTHING
    }
  }

  // Use (inexact) CG to solve Newton system 
  void CG( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
           Objective<Real> &obj, Teuchos::RCP<Secant<Real> > secant = Teuchos::null ) {
    Real gnorm = g.norm(); 
    Real gtol = std::min(tol1_,tol2_*gnorm);
    Real itol = 0.0;

    s.zero(); 

    Teuchos::RCP<Vector<Real> > gnew = x.clone(); 
    gnew->set(g); 

    Teuchos::RCP<Vector<Real> > v = x.clone();  
    if ( secant != Teuchos::null ) {
      secant->applyB( *v, *gnew, x );
    }
    else {
      obj.precond( *v, *gnew, x );  
    }

    Teuchos::RCP<Vector<Real> > p = x.clone(); 
    p->set(*v); 

    Teuchos::RCP<Vector<Real> > Hp = x.clone();
    itol = 0.0;
    if ( this->useInexact_ ) {
      itol = gtol/(maxit_ * gnorm); 
      //itol = std::min( 0.5, gtol*p->norm()/(2.0 * this->maxit_ * gnorm * gnorm) ); 
    }
    obj.hessVec( *Hp, *p, x, itol );  

    iter = 0; 
    flag = 0;

    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v->dot(*gnew); 

    for (iter = 0; iter < this->maxit_; iter++) {
      kappa = p->dot(*Hp);
      if ( kappa <= 0.0 ) { 
        flag = 2;
        break;
      }
      alpha = gv/kappa;

      s.axpy(alpha,*p);

      gnew->axpy(-alpha,*Hp);
      gnorm = gnew->norm();
      if ( gnorm < gtol ) {
        break;
      }

      if ( secant != Teuchos::null ) {
        secant->applyB( *v, *gnew, x );
      }
      else {
        obj.precond( *v, *gnew, x );  
      }
      tmp  = gv;         
      gv   = v->dot(*gnew); 
      beta = gv/tmp;
 
      p->scale(beta);
      p->axpy(1.0,*v);

      itol = 0.0;
      if ( this->useInexact_ ) {
        itol = gtol/(this->maxit_ * gnorm); 
        //itol = std::min( 0.5, gtol*p->norm()/(2.0 * this->maxit_ * gnorm * gnorm) ); 
      }
      obj.hessVec( *Hp, *p, x, itol );
    }
    iter++;
    if ( iter == this->maxit_ ) {
      flag = 1;
    }    
  }

  // Use (inexact) CR to solve Newton system 
  void CR( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
           Objective<Real> &obj, Teuchos::RCP<Secant<Real> > secant = Teuchos::null ) {
    // Initialize
    Real gnorm = g.norm(); 
    Real gtol = std::min(tol1_,tol2_*gnorm);
    Real itol = 0.0;
    s.zero(); 

    // Apply preconditioner to residual
    Teuchos::RCP<Vector<Real> > v = x.clone();  
    if ( secant != Teuchos::null ) {
      secant->applyB( *v, g, x );
    }
    else {
      obj.precond( *v, g, x );  
    }

    // Initialize direction p
    Teuchos::RCP<Vector<Real> > p = x.clone(); 
    p->set(*v); 

    // Get Hessian tolerance
    itol = 0.0;
    if ( this->useInexact_ ) {
      itol = gtol/(maxit_ * gnorm); 
    }
    // Apply Hessian to residual
    Teuchos::RCP<Vector<Real> > Hv  = x.clone();
    obj.hessVec( *Hv, *v, x, itol );
    // Apply Hessian to direction p
    Teuchos::RCP<Vector<Real> > Hp  = x.clone();
    Teuchos::RCP<Vector<Real> > MHp = x.clone();  
    obj.hessVec( *Hp, *p, x, itol );  

    // Initialize scalar quantities
    iter = 0; 
    flag = 0;
    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real vHv   = Hv->dot(*v); 

    for (iter = 0; iter < this->maxit_; iter++) {
      if ( secant != Teuchos::null ) {
        secant->applyB( *MHp, *Hp, x );
      }
      else {
        obj.precond( *MHp, *Hp, x );  
      }
      kappa = Hp->dot(*MHp);
      if ( vHv <= 0.0 || kappa <= 0.0 ) { 
        flag = 2;
        break;
      }
      alpha = vHv/kappa;

      s.axpy(alpha,*p);

      v->axpy(-alpha,*MHp);
      gnorm = v->norm();
      if ( gnorm < gtol ) {
        break;
      }

      itol = 0.0;
      if ( this->useInexact_ ) {
        itol = gtol/(this->maxit_ * gnorm); 
      }
      obj.hessVec( *Hv, *v, x, itol );
      tmp  = vHv;
      vHv  = Hv->dot(*v);
      beta = vHv/tmp;

      p->scale(beta);
      p->axpy(1.0,*v);

      Hp->scale(beta);
      Hp->axpy(1.0,*Hv); 
    }
    iter++;
    if ( iter == this->maxit_ ) {
      flag = 1;
    }    
  }
};


}

#endif
