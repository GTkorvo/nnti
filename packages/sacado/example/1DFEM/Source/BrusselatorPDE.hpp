// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef BRUSSELATORPDE_HPP
#define BRUSSELATORPDE_HPP

#include "Teuchos_RefCountPtr.hpp"

#include "AbstractPDE.hpp"
#include "AbstractSourceFunction.hpp"

template <typename ScalarT>
class BrusselatorPDE : public AbstractPDE<ScalarT> {
public:
  
  //! Constructor
  BrusselatorPDE(double alpha, double beta, double D1, double D2);

  //! Destructor
  virtual ~BrusselatorPDE();

  //! Number of discretized equations
  virtual unsigned int numEquations() const;

  //! Initialize PDE
  virtual void init(unsigned int numQuadPoints, unsigned int numNodes);

  //! Evaluate discretized PDE element-level residual
  virtual void
  evaluateElementResidual(const AbstractQuadrature& quadRule,
			  const AbstractElement& element,
			  const std::vector<ScalarT>& solution,
			  std::vector<ScalarT>& residual);

private:

  //! Private to prohibit copying
  BrusselatorPDE(const BrusselatorPDE&);

  //! Private to prohibit copying
  BrusselatorPDE& operator=(const BrusselatorPDE&);

protected:

  //! Number of quad points
  unsigned int num_qp;

  //! Number of nodes
  unsigned int num_nodes;

  //! Shape function values
  std::vector< std::vector<double> > phi;

  //! Shape function derivatives
  std::vector< std::vector<double> > dphi;

  //! Element transformation Jacobian
  std::vector<double> jac;

  //! Discretized solution
  std::vector<ScalarT> T;

  //! Discretized solution
  std::vector<ScalarT> C;

  //! Discretized solution derivative
  std::vector<ScalarT> dT;

  //! Discretized solution derivative
  std::vector<ScalarT> dC;

  //! Model parameters
  double alpha, beta, D1, D2;

};

// Include implementation
#include "BrusselatorPDEImpl.hpp"

#endif // BRUSSELATORPDE_HPP
