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

#ifndef ROL_BARZILAIBORWEIN_H
#define ROL_BARZILAIBORWEIN_H

/** \class ROL::BarzilaiBorwein
    \brief Provides definitions for Barzilai-Borwein operators.
*/

namespace ROL {

template<class Real>
class BarzilaiBorwein : public Secant<Real> {
private:

  int type_;

public:
  BarzilaiBorwein(int type = 1) : Secant<Real>(1), type_(type) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x ) {
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Hv.set(v);
    if ( state->iter != 0 && state->current != -1 ) {
      if ( type_ == 1 ) {
        Real yy = state->gradDiff[state->current]->dot(*(state->gradDiff[state->current]));
        Hv.scale(state->product[state->current]/yy);
      }
      else if ( type_ == 2 ) {
        Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
        Hv.scale(ss/state->product[state->current]);
      }
    }
  }

  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v, const Vector<Real> &x ) { 
    // Get Generic Secant State
    Teuchos::RCP<SecantState<Real> >& state = Secant<Real>::get_state();

    Bv.set(v);
    if ( state->iter != 0 && state->current != -1 ) {
      if ( type_ == 1 ) {
        Real yy = state->gradDiff[state->current]->dot(*(state->gradDiff[state->current]));
        Bv.scale(yy/state->product[state->current]);
      }
      else if ( type_ == 2 ) {
        Real ss = state->iterDiff[state->current]->dot(*(state->iterDiff[state->current]));
        Bv.scale(state->product[state->current]/ss);
      }
    }
  }

};

}

#endif
