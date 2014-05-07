//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef SIMPLEAPP_EPETRACONSTRAINTS_H
#define SIMPLEAPP_EPETRACONSTRAINTS_H

#include "Aristos_EpetraVector.hpp"
#include "Aristos_Constraints.hpp"
#include "SimpleApp_SimpleEpetraDataPool.hpp"

namespace SimpleApp {

class SimpleEpetraConstraints : public Aristos::Constraints {
private:

  Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat_;
  
public:

  SimpleEpetraConstraints( Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool>  dat );
  
  void getValue( const Aristos::Vector &x, Aristos::Vector &c ) const;

  void applyJacobian( bool Trans, const Aristos::Vector &x, const Aristos::Vector &v, Aristos::Vector &Jv) const;

  void applyNullSp( bool Trans, const Aristos::Vector &x, const Aristos::Vector &v, Aristos::Vector &Wv,
                    double * tol) const;

};

} // namespace SimpleApp

#endif
