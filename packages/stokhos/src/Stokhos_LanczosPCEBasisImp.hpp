// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_TestForException.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_Lanczos.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(
   ordinal_type p,
   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
   const Stokhos::Quadrature<ordinal_type, value_type>& quad,
   bool normalize) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, normalize)
{
  // Evaluate PCE at quad points
  pce_weights = quad.getQuadWeights();
  ordinal_type nqp = pce_weights.size();
  pce_vals.resize(nqp);
  h0.resize(nqp);
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad.getQuadPoints();
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad.getBasisAtQuadPoints();
  for (ordinal_type i=0; i<nqp; i++) {
    pce_vals[i] = pce.evaluate(quad_points[i], basis_values[i]);
    h0[i] = value_type(1);
  }
  
  // Compute coefficients via Lanczos
  Teuchos::Array< Teuchos::Array<value_type> > h;
  Stokhos::Lanczos<ordinal_type,value_type>::computeDiag(
    nqp, p+1, pce_vals, pce_weights, h0, h, this->alpha, this->beta,
    this->norms);
  // lanczos(nqp, p+1, pce_weights, pce_vals, h0, this->alpha, this->beta,
  // 	  this->norms);
  for (ordinal_type i=0; i<=p; i++)
    this->delta[i] = value_type(1.0);
  
  // Setup rest of recurrence basis
  //this->setup();
  this->gamma[0] = value_type(1);
  for (ordinal_type k=1; k<=p; k++) {
    this->gamma[k] = value_type(1);
  } 

  //If you want normalized polynomials, set gamma and reset the norms to 1.
  if( normalize ) {
    for (ordinal_type k=0; k<=p; k++) {
      this->gamma[k] = value_type(1)/std::sqrt(this->norms[k]);
      this->norms[k] = value_type(1);
    }
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
~LanczosPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::LanczosPCEBasis -- compute Gauss points");
#endif

  // Call base class
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));

  // We can't reliably generate quadrature points of order > 2*p
  //std::cout << "quad_order = " << quad_order << ", 2*p = " << 2*this->p << std::endl;
  // if (quad_order > 2*this->p)
  //   quad_order = 2*this->p;
  Stokhos::RecurrenceBasis<ordinal_type,value_type>::getQuadPoints(quad_order, 
								   quad_points, 
								   quad_weights,
								   quad_values);

  // Fill in the rest of the points with zero weight
  if (quad_weights.size() < num_points) {
    ordinal_type old_size = quad_weights.size();
    quad_weights.resize(num_points);
    quad_points.resize(num_points);
    quad_values.resize(num_points);
    for (ordinal_type i=old_size; i<num_points; i++) {
      quad_weights[i] = value_type(0);
      quad_points[i] = quad_points[0];
      quad_values[i].resize(this->p+1);
      evaluateBases(quad_points[i], quad_values[i]);
    }
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::LanczosPCEBasis<ordinal_type,value_type>(
			 p,*this));
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta) const
{
  ordinal_type nqp = pce_weights.size();
  Teuchos::Array<value_type> nrm(n);
  Teuchos::Array< Teuchos::Array<value_type> > h;
  Stokhos::Lanczos<ordinal_type,value_type>::computeDiag(
    nqp, n, pce_vals, pce_weights, h0, h, alpha, beta, nrm);
  //lanczos(nqp, n, pce_weights, pce_vals, h0, alpha, beta, nrm);
  for (ordinal_type i=0; i<n; i++) {
    delta[i] = value_type(1.0);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
lanczos(ordinal_type n,
	ordinal_type nsteps,
	const Teuchos::Array<value_type>& w,
	const Teuchos::Array<value_type>& A,
	const Teuchos::Array<value_type>& h0,
	Teuchos::Array<value_type>& a,
	Teuchos::Array<value_type>& b,
	Teuchos::Array<value_type>& nrm_sqrd) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::LanczosPCEBasis -- Lanczos Procedure");
#endif

  Teuchos::Array<value_type> h(n), hm(n), y(n);
  value_type dp, nrm;
  for (ordinal_type j=0; j<n; j++) {
    h[j] = h0[j];
    hm[j] = value_type(0);
  }

  b[0] = 1.0;
  for (ordinal_type i=0; i<nsteps; i++) {

    // compute h^T*h
    nrm = value_type(0);
    for (ordinal_type j=0; j<n; j++)
      nrm += w[j]*h[j]*h[j];
    nrm_sqrd[i] = nrm;

    // compute y = A*h
    for (ordinal_type j=0; j<n; j++)
      y[j] = A[j]*h[j];

    // compute alpha = h^T*y / h^T*h
    dp = value_type(0);
    for (ordinal_type j=0; j<n; j++)
      dp += w[j]*y[j]*h[j];
    a[i] = dp/nrm_sqrd[i];

    // compute beta = h^T*h / hm^T*hm
    if (i > 0)
      b[i] = nrm_sqrd[i] / nrm_sqrd[i-1];

    // compute
    for (ordinal_type j=0; j<n; j++)
      y[j] = y[j] - a[i]*h[j] - b[i]*hm[j];

    // shift
    hm = h;
    h = y;

    std::cout << "i = " << i << " alpha = " << a[i] << " beta = " << b[i]
	      << " nrm = " << nrm_sqrd[i] << std::endl;
    TEST_FOR_EXCEPTION(nrm_sqrd[i] < 0, std::logic_error,
		       "Stokhos::LanczosPCEBasis::lanczos():  "
		       << " Polynomial " << i << " out of " << nsteps
		       << " has norm " << nrm_sqrd[i]
		       << "!  Try increasing number of quadrature points");
  }
}

template <typename ordinal_type, typename value_type> 
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(ordinal_type p, const LanczosPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, false),
  pce_weights(basis.pce_weights),
  pce_vals(basis.pce_vals),
  h0(basis.h0)
{
}
