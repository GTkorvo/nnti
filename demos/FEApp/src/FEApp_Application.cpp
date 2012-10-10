// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include "FEApp_Application.hpp"
#include "FEApp_ProblemFactory.hpp"
#include "FEApp_QuadratureFactory.hpp"
#include "FEApp_DiscretizationFactory.hpp"
#include "Epetra_LocalMap.h"
#include "FEApp_SGGaussQuadResidualGlobalFill.hpp"
#include "FEApp_SGGaussQuadJacobianGlobalFill.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Assert.hpp"

FEApp::Application::Application(
  const Teuchos::RCP<const Epetra_Comm>& comm,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  const Epetra_Vector* initial_soln) :
  params(params_)
{
  // Create parameter library
  paramLib = Teuchos::rcp(new ParamLib);

  // Create problem object
  Teuchos::RCP<Teuchos::ParameterList> problemParams = 
    Teuchos::rcp(&(params->sublist("Problem")),false);
  FEApp::ProblemFactory problemFactory(problemParams, paramLib);
  Teuchos::RCP<FEApp::AbstractProblem> problem = 
    problemFactory.create();

  // Get number of equations
  unsigned int num_equations = problem->numEquations();

  // Create quadrature object
  Teuchos::RCP<Teuchos::ParameterList> quadParams = 
    Teuchos::rcp(&(params->sublist("Quadrature")),false);
  FEApp::QuadratureFactory quadFactory(quadParams);
  quad = quadFactory.create();

  // Create mesh
  Teuchos::RCP<Teuchos::ParameterList> discParams = 
    Teuchos::rcp(&(params->sublist("Discretization")),false);
  int nelem = discParams->get("Number of Elements", 100);
  double h = 1.0/nelem;
  std::vector<double> coords(nelem+1);
  for (int i=0; i<=nelem; i++)
    coords[i] = h*i;

  // Create discretization object
  FEApp::DiscretizationFactory discFactory(discParams);
  disc = discFactory.create(coords, num_equations, comm);
  disc->createMesh();
  disc->createMaps();
  disc->createJacobianGraphs();

  // Create Epetra objects
  importer = Teuchos::rcp(new Epetra_Import(*(disc->getOverlapMap()), 
                                            *(disc->getMap())));
  exporter = Teuchos::rcp(new Epetra_Export(*(disc->getOverlapMap()), 
                                            *(disc->getMap())));
  overlapped_x = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_xdot = 
    Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_f = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_jac = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, 
                                      *(disc->getOverlapJacobianGraph())));

  // Initialize problem
  initial_x = Teuchos::rcp(new Epetra_Vector(*(disc->getMap())));
  problem->buildProblem(*(disc->getMap()), *(disc->getOverlapMap()), 
                        pdeTM, bc, responses, initial_x);
  if (initial_soln != NULL)
    initial_x->Scale(1.0, *initial_soln);
  typedef FEApp::AbstractPDE_TemplateManager<EvalTypes>::iterator iterator;
  int nqp = quad->numPoints();
  int nn = disc->getNumNodesPerElement();
  for (iterator it = pdeTM.begin(); it != pdeTM.end(); ++it)
    it->init(nqp, nn);

  // Create response map
  unsigned int total_num_responses = 0;
  for (unsigned int i=0; i<responses.size(); i++)
    total_num_responses += responses[i]->numResponses();
  if (total_num_responses > 0)
    response_map = Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(total_num_responses), 0,
                                                    *comm));
}

FEApp::Application::~Application()
{
}

Teuchos::RCP<const Epetra_Map>
FEApp::Application::getMap() const
{
  return disc->getMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::Application::getResponseMap() const
{
  return response_map;
}

Teuchos::RCP<const Epetra_CrsGraph>
FEApp::Application::getJacobianGraph() const
{
  return disc->getJacobianGraph();
}

Teuchos::RCP<const Epetra_Vector>
FEApp::Application::getInitialSolution() const
{
  return initial_x;
}

Teuchos::RCP<ParamLib> 
FEApp::Application::getParamLib()
{
  return paramLib;
}

void
FEApp::Application::init_sg(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
      const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad,
      const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& exp,
      const Teuchos::RCP<const EpetraExt::MultiComm>& multiComm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(basis == Teuchos::null, std::logic_error,
		     "Error!  FEApp::Application::init_sg():  " <<
		     "SG basis cannot be null!");
  TEUCHOS_TEST_FOR_EXCEPTION(exp == Teuchos::null, std::logic_error,
		     "Error!  FEApp::Application::init_sg():  " <<
		     "SG expansion cannot be null!");
  sg_basis = basis;
  sg_quad = quad;
  sg_expansion = exp;

  // Create Epetra orthogonal polynomial objects
  if (sg_overlapped_x == Teuchos::null) {
    sg_overlap_map =
      Teuchos::rcp(new Epetra_LocalMap(sg_basis->size(), 0, 
				       multiComm->TimeDomainComm()));

    sg_overlapped_x = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     sg_basis, sg_overlap_map, disc->getOverlapMap(),
		     multiComm));
    sg_overlapped_xdot = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     sg_basis, sg_overlap_map, disc->getOverlapMap(),
		     multiComm));
    sg_overlapped_f = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     sg_basis, sg_overlap_map, disc->getOverlapMap(),
		     multiComm));
    // Delay creation of sg_overlapped_jac until needed
  }

  for (unsigned int i=0; i<responses.size(); i++)
    responses[i]->init_sg(basis, quad, exp);
}

Teuchos::RCP<Epetra_Operator>
FEApp::Application::createW() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
					   *(disc->getJacobianGraph())));
}

Teuchos::RCP<Epetra_Operator>
FEApp::Application::createPrec() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
					   *(disc->getJacobianGraph())));
}

void
FEApp::Application::computeGlobalResidual(const Epetra_Vector* xdot,
					  const Epetra_Vector& x,
					  const Teuchos::Array<ParamVec>& p,
					  Epetra_Vector& f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalResidual");

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (xdot)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Zero out overlapped residual
  overlapped_f->PutScalar(0.0);
  f.PutScalar(0.0);

  // Create residual init/post op
  Teuchos::RCP<FEApp::ResidualOp> op = 
    Teuchos::rcp(new FEApp::ResidualOp(overlapped_xdot, overlapped_x, 
                                       overlapped_f));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::ResidualType> > pde = 
    pdeTM.getAsObject<FEApp::ResidualType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::ResidualType> globalFill(disc->getMesh(), quad, 
                                                    pde, bc, xdot!=NULL);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  f.Export(*overlapped_f, *exporter, Add);
}

void
FEApp::Application::computeGlobalJacobian(double alpha, double beta,
					  const Epetra_Vector* xdot,
					  const Epetra_Vector& x,
					  const Teuchos::Array<ParamVec>& p,
					  Epetra_Vector* f,
					  Epetra_Operator& jacOp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalJacobian");

  // Cast jacOp to an Epetra_CrsMatrix
  Epetra_CrsMatrix& jac = dynamic_cast<Epetra_CrsMatrix&>(jacOp);

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (xdot)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Zero out overlapped residual
  Teuchos::RCP<Epetra_Vector> overlapped_ff;
  if (f != NULL) {
    overlapped_ff = overlapped_f;
    overlapped_ff->PutScalar(0.0);
    f->PutScalar(0.0);
  }

  // Zero out Jacobian
  overlapped_jac->PutScalar(0.0);
  jac.PutScalar(0.0);

  // Create Jacobian init/post op
  Teuchos::RCP<FEApp::JacobianOp> op
    = Teuchos::rcp(new FEApp::JacobianOp(alpha, beta, overlapped_xdot, 
                                         overlapped_x, overlapped_ff, 
                                         overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > pde = 
    pdeTM.getAsObject<FEApp::JacobianType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::JacobianType> globalFill(disc->getMesh(), 
                                                    quad, pde, bc, 
                                                    xdot!=NULL);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (f != NULL)
    f->Export(*overlapped_f, *exporter, Add);

  // Assemble global Jacobian
  jac.Export(*overlapped_jac, *exporter, Add);

  jac.FillComplete(true);
}

void
FEApp::Application::computeGlobalPreconditioner(
  double alpha, double beta,
  const Epetra_Vector* xdot,
  const Epetra_Vector& x,
  const Teuchos::Array<ParamVec>& p,
  Epetra_Vector* f,
  Epetra_Operator& jacOp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalPreconditioner");

  // Cast jacOp to an Epetra_CrsMatrix
  Epetra_CrsMatrix& jac = dynamic_cast<Epetra_CrsMatrix&>(jacOp);

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (xdot)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Zero out overlapped residual
  Teuchos::RCP<Epetra_Vector> overlapped_ff;
  if (f != NULL) {
    overlapped_ff = overlapped_f;
    overlapped_ff->PutScalar(0.0);
    f->PutScalar(0.0);
  }

  // Zero out Jacobian
  overlapped_jac->PutScalar(0.0);
  jac.PutScalar(0.0);

  // Create Jacobian init/post op
  Teuchos::RCP<FEApp::JacobianOp> op = 
    Teuchos::rcp(new FEApp::JacobianOp(alpha, beta, overlapped_xdot, 
                                       overlapped_x, overlapped_ff, 
                                       overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > pde = 
    pdeTM.getAsObject<FEApp::JacobianType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::JacobianType> globalFill(disc->getMesh(), 
                                                    quad, pde, bc, 
                                                    xdot!=NULL);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (f != NULL)
    f->Export(*overlapped_f, *exporter, Add);

  // Assemble global Jacobian
  jac.Export(*overlapped_jac, *exporter, Add);

  jac.FillComplete(true);
}

void
FEApp::Application::computeGlobalTangent(
  double alpha, double beta,
  bool sum_derivs,
  const Epetra_Vector* xdot,
  const Epetra_Vector& x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::RCP<ParamVec>& deriv_p,
  const Epetra_MultiVector* Vx,
  const Teuchos::SerialDenseMatrix<int,double>* Vp,
  Epetra_Vector* f,
  Epetra_MultiVector* JVx,
  Epetra_MultiVector* fVp)
{
  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (xdot)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Zero out overlapped residual
  Teuchos::RCP<Epetra_Vector> overlapped_ff;
  if (f != NULL) {
    overlapped_ff = overlapped_f;
    overlapped_ff->PutScalar(0.0);
    f->PutScalar(0.0);
  }

  Teuchos::RCP<Epetra_MultiVector> overlapped_JVx;
  if (JVx != NULL) {
    overlapped_JVx = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          JVx->NumVectors()));
    overlapped_JVx->PutScalar(0.0);
    JVx->PutScalar(0.0);
  }
  
  Teuchos::RCP<Epetra_MultiVector> overlapped_fVp;
  if (fVp != NULL) {
    overlapped_fVp = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          fVp->NumVectors()));
    overlapped_fVp->PutScalar(0.0);
    fVp->PutScalar(0.0);
  }

  Teuchos::RCP<Epetra_MultiVector> overlapped_Vx;
  if (Vx != NULL) {
    overlapped_Vx = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          Vx->NumVectors()));
  }

  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> > vp =
    Teuchos::rcp(Vp, false);
  
  // Create Jacobian init/post op
  Teuchos::RCP<FEApp::TangentOp> op = 
    Teuchos::rcp(new FEApp::TangentOp(alpha, beta, sum_derivs,
                                      overlapped_xdot, 
                                      overlapped_x,
                                      deriv_p,
                                      overlapped_Vx,
                                      overlapped_Vx,
                                      vp,
                                      overlapped_ff, 
                                      overlapped_JVx,
                                      overlapped_fVp));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::TangentType> > pde = 
    pdeTM.getAsObject<FEApp::TangentType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::TangentType> globalFill(disc->getMesh(), 
                                                   quad, pde, bc, 
                                                   xdot!=NULL);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (f != NULL)
    f->Export(*overlapped_f, *exporter, Add);

  // Assemble derivatives
  if (JVx != NULL)
    JVx->Export(*overlapped_JVx, *exporter, Add);
  if (fVp != NULL)
    fVp->Export(*overlapped_fVp, *exporter, Add);
}

void
FEApp::Application::
evaluateResponses(const Epetra_Vector* xdot,
                  const Epetra_Vector& x,
                  const Teuchos::Array<ParamVec>& p,
                  Epetra_Vector& g)
{
  const Epetra_Comm& comm = x.Map().Comm();
  unsigned int offset = 0;
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(static_cast<int>(num_responses), 0, comm);

    // Create Epetra_Vector for response function
    Epetra_Vector local_g(local_response_map);

    // Evaluate response function
    responses[i]->evaluateResponses(xdot, x, p, local_g);

    // Copy result into combined result
    for (unsigned int j=0; j<num_responses; j++)
      g[offset+j] = local_g[j];

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::
evaluateResponseTangents(
  const Epetra_Vector* xdot,
  const Epetra_Vector& x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
  Epetra_Vector* g,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& gt)
{
  const Epetra_Comm& comm = x.Map().Comm();
  unsigned int offset = 0;
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > local_gt(gt.size());
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(static_cast<int>(num_responses), 0, comm);

    // Create Epetra_Vectors for response function
    Teuchos::RCP<Epetra_Vector> local_g;
    if (g != NULL)
      local_g = Teuchos::rcp(new Epetra_Vector(local_response_map));
    for (int j=0; j<gt.size(); j++)
      if (gt[j] != Teuchos::null)
	local_gt[j] = Teuchos::rcp(new Epetra_MultiVector(local_response_map, 
							  gt[j]->NumVectors()));

    // Evaluate response function
    responses[i]->evaluateTangents(xdot, x, p, deriv_p, dxdot_dp, dx_dp, 
				   local_g.get(), local_gt);

    // Copy results into combined result
    for (unsigned int j=0; j<num_responses; j++) {
      if (g != NULL)
        (*g)[offset+j] = (*local_g)[j];
      for (int l=0; l<gt.size(); l++)
	if (gt[l] != Teuchos::null)
	  for (int k=0; k<gt[l]->NumVectors(); k++)
	    (*gt[l])[k][offset+j] = (*local_gt[l])[k][j];
    }

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::
evaluateResponseGradients(
  const Epetra_Vector* xdot,
  const Epetra_Vector& x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  Epetra_Vector* g,
  Epetra_MultiVector* dg_dx,
  Epetra_MultiVector* dg_dxdot,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dg_dp)
{
  const Epetra_Comm& comm = x.Map().Comm();
  unsigned int offset = 0;
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > local_dgdp(dg_dp.size());
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(static_cast<int>(num_responses), 0, comm);

    // Create Epetra_Vectors for response function
    Teuchos::RCP<Epetra_Vector> local_g;
    if (g != NULL)
      local_g = Teuchos::rcp(new Epetra_Vector(local_response_map));
    Teuchos::RCP<Epetra_MultiVector> local_dgdx;
    if (dg_dx != NULL)
      local_dgdx = Teuchos::rcp(new Epetra_MultiVector(dg_dx->Map(), 
                                                       num_responses));
    Teuchos::RCP<Epetra_MultiVector> local_dgdxdot;
    if (dg_dxdot != NULL)
      local_dgdxdot = Teuchos::rcp(new Epetra_MultiVector(dg_dxdot->Map(), 
                                                          num_responses));

    for (int j=0; j<dg_dp.size(); j++)
      if (dg_dp[j] != Teuchos::null)
	local_dgdp[j] = Teuchos::rcp(new Epetra_MultiVector(local_response_map, 
							    dg_dp[j]->NumVectors()));

    // Evaluate response function
    responses[i]->evaluateGradients(xdot, x, p, deriv_p, local_g.get(), 
                                    local_dgdx.get(), local_dgdxdot.get(), 
                                    local_dgdp);

    // Copy results into combined result
    for (unsigned int j=0; j<num_responses; j++) {
      if (g != NULL)
        (*g)[offset+j] = (*local_g)[j];
      if (dg_dx != NULL)
        (*dg_dx)(offset+j)->Update(1.0, *((*local_dgdx)(j)), 0.0);
      if (dg_dxdot != NULL)
        (*dg_dxdot)(offset+j)->Update(1.0, *((*local_dgdxdot)(j)), 0.0);
      for (int l=0; l<dg_dp.size(); l++)
	if (dg_dp[l] != Teuchos::null)
	  for (int k=0; k<dg_dp[l]->NumVectors(); k++)
	    (*dg_dp[l])[k][offset+j] = (*local_dgdp[l])[k][j];
    }

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::computeGlobalSGResidual(
  const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
  const Stokhos::EpetraVectorOrthogPoly& sg_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& sg_p_index,
  const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
  Stokhos::EpetraVectorOrthogPoly& sg_f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGResidual");

  for (int i=0; i<sg_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*sg_overlapped_x)[i].Import(sg_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (sg_xdot)
      (*sg_overlapped_xdot)[i].Import((*sg_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    (*sg_overlapped_f)[i].PutScalar(0.0);
    sg_f[i].PutScalar(0.0);

  }

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Set SG parameters
  for (int i=0; i<sg_p_index.size(); i++) {
    int ii = sg_p_index[i];
    for (unsigned int j=0; j<p[ii].size(); j++)
	p[ii][j].family->setValue<FEApp::SGResidualType>(sg_p_vals[ii][j]);
    }

  // Create residual init/post op
  Teuchos::RCP<FEApp::SGResidualOp> sg_res_fill_op = 
    Teuchos::rcp(new FEApp::SGResidualOp(sg_expansion, 
					 sg_overlapped_xdot, 
					 sg_overlapped_x, 
					 sg_overlapped_f));
    
  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGResidualType> > pde = 
    pdeTM.getAsObject<FEApp::SGResidualType>();

  // Instantiate global fill
  if (sg_res_global_fill == Teuchos::null) {
    std::string method = params->get("SG Method", "AD");
    if (method == "AD") {
      sg_res_global_fill = 
        Teuchos::rcp(new FEApp::GlobalFill<FEApp::SGResidualType>(
	  disc->getMesh(), quad, pde, bc, sg_xdot!=NULL));
    }
    else if (method == "Gauss Quadrature") {
      Teuchos::RCP< FEApp::AbstractPDE<FEApp::ResidualType> > res_pde = 
	pdeTM.getAsObject<FEApp::ResidualType>();
      Teuchos::Array< Teuchos::RCP<const ParamVec> >sg_p(sg_p_index.size());
      for (int i=0; i<sg_p_index.size(); i++)
	sg_p[i] = Teuchos::rcp(&p[sg_p_index[i]], false);
      sg_res_global_fill = 
        Teuchos::rcp(new FEApp::SGGaussQuadResidualGlobalFill(disc->getMesh(), 
                                                              quad, pde, bc, 
                                                              sg_xdot!=NULL,
                                                              sg_basis,
							      sg_quad,
                                                              res_pde,
                                                              sg_p));
    }
  }

  // Do global fill
  sg_res_global_fill->computeGlobalFill(*sg_res_fill_op);

  // Assemble global residual
  for (int i=0; i<sg_f.size(); i++)
    sg_f[i].Export((*sg_overlapped_f)[i], *exporter, Add);
}

void
FEApp::Application::computeGlobalSGJacobian(
  double alpha, double beta,
  const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
  const Stokhos::EpetraVectorOrthogPoly& sg_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& sg_p_index,
  const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
  Stokhos::EpetraVectorOrthogPoly* sg_f,
  Stokhos::EpetraOperatorOrthogPoly& sg_jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGJacobian");

  for (int i=0; i<sg_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*sg_overlapped_x)[i].Import(sg_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (sg_xdot)
      (*sg_overlapped_xdot)[i].Import((*sg_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    if (sg_f != NULL) {
      (*sg_overlapped_f)[i].PutScalar(0.0);
      (*sg_f)[i].PutScalar(0.0);
    }

  }

  // Create, resize and initialize overlapped Jacobians
  if (sg_overlapped_jac == Teuchos::null || 
      sg_overlapped_jac->size() < sg_jac.size()) {
    sg_overlap_jac_map = 
      Teuchos::rcp(new Epetra_LocalMap(sg_basis->size(), 0, 
				       sg_overlap_map->Comm()));
    sg_overlapped_jac = 
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_CrsMatrix>(
		     sg_basis,  sg_overlap_jac_map, *overlapped_jac));
  }
  for (int i=0; i<sg_overlapped_jac->size(); i++)
    (*sg_overlapped_jac)[i].PutScalar(0.0);

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Set SG parameters
  for (int i=0; i<sg_p_index.size(); i++) {
    int ii = sg_p_index[i];
    for (unsigned int j=0; j<p[ii].size(); j++)
	p[ii][j].family->setValue<FEApp::SGJacobianType>(sg_p_vals[ii][j]);
    }

  // Create Jacobian init/post op
  Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > sg_overlapped_ff;
  if (sg_f != NULL)
    sg_overlapped_ff = sg_overlapped_f;
  Teuchos::RCP<FEApp::SGJacobianOp> sg_jac_fill_op = 
    Teuchos::rcp(new FEApp::SGJacobianOp(sg_expansion,
					 alpha, beta, 
					 sg_overlapped_xdot, 
					 sg_overlapped_x, 
					 sg_overlapped_ff, 
					 sg_overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGJacobianType> > pde = 
    pdeTM.getAsObject<FEApp::SGJacobianType>();

  // Instantiate global fill
  if (sg_jac_global_fill == Teuchos::null) {
    std::string method = params->get("SG Method", "AD");
    if (method == "AD") {
      sg_jac_global_fill = 
	Teuchos::rcp(new FEApp::GlobalFill<FEApp::SGJacobianType>(
	    disc->getMesh(), quad, pde, bc,sg_xdot!=NULL));
    }
    else if (method == "Gauss Quadrature") {
      Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > jac_pde = 
	pdeTM.getAsObject<FEApp::JacobianType>();
      Teuchos::Array< Teuchos::RCP<const ParamVec> >sg_p(sg_p_index.size());
      for (int i=0; i<sg_p_index.size(); i++)
	sg_p[i] = Teuchos::rcp(&p[sg_p_index[i]], false);
      sg_jac_global_fill = 
	Teuchos::rcp(new FEApp::SGGaussQuadJacobianGlobalFill(disc->getMesh(),
							      quad, pde, bc, 
							      sg_xdot!=NULL,
							      sg_basis,
							      sg_quad,
							      jac_pde,
							      sg_p, 
							      alpha, beta));
    }
  }

  // Do global fill
  sg_jac_global_fill->computeGlobalFill(*sg_jac_fill_op);
  
  // Assemble global residual
  if (sg_f != NULL)
    for (int i=0; i<sg_f->size(); i++)
      (*sg_f)[i].Export((*sg_overlapped_f)[i], *exporter, Add);
    
  // Assemble block Jacobians
  Teuchos::RCP<Epetra_CrsMatrix> jac;
  for (int i=0; i<sg_jac.size(); i++) {
    jac = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(sg_jac.getCoeffPtr(i), 
						      true);
    jac->PutScalar(0.0);
    jac->Export((*sg_overlapped_jac)[i], *exporter, Add);
    jac->FillComplete(true);
  }
}

void
FEApp::Application::computeGlobalSGTangent(
  double alpha, double beta, bool sum_derivs,
  const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
  const Stokhos::EpetraVectorOrthogPoly& sg_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& sg_p_index,
  const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
  const Teuchos::RCP<ParamVec>& deriv_p,
  const Epetra_MultiVector* Vx,
  const Teuchos::SerialDenseMatrix<int,double>* Vp,
  Stokhos::EpetraVectorOrthogPoly* sg_f,
  Stokhos::EpetraMultiVectorOrthogPoly* sg_JVx,
  Stokhos::EpetraMultiVectorOrthogPoly* sg_fVp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGTangent");

  for (int i=0; i<sg_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*sg_overlapped_x)[i].Import(sg_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (sg_xdot)
      (*sg_overlapped_xdot)[i].Import((*sg_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    if (sg_f != NULL) {
      (*sg_overlapped_f)[i].PutScalar(0.0);
      (*sg_f)[i].PutScalar(0.0);
    }

  }

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Set SG parameters
  for (int i=0; i<sg_p_index.size(); i++) {
    int ii = sg_p_index[i];
    for (unsigned int j=0; j<p[ii].size(); j++)
	p[ii][j].family->setValue<FEApp::SGTangentType>(sg_p_vals[ii][j]);
  }


  Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > sg_overlapped_JVx;
  if (sg_JVx != NULL) {
    sg_overlapped_JVx = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		     sg_basis, sg_overlap_map, disc->getOverlapMap(),
		     sg_x.productComm(),
		     (*sg_JVx)[0].NumVectors()));
    sg_JVx->init(0.0);
  }
  
  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly > sg_overlapped_fVp;
  if (sg_fVp != NULL) {
    sg_overlapped_fVp = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		     sg_basis, sg_overlap_map, disc->getOverlapMap(),
		     sg_x.productComm(), 
		     (*sg_fVp)[0].NumVectors()));
    sg_fVp->init(0.0);
  }

  Teuchos::RCP<Epetra_MultiVector> overlapped_Vx;
  if (Vx != NULL) {
    overlapped_Vx = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          Vx->NumVectors()));
  }

  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> > vp =
    Teuchos::rcp(Vp, false);

  // Create Jacobian init/post op
   Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > sg_overlapped_ff;
  if (sg_f != NULL)
    sg_overlapped_ff = sg_overlapped_f;
  Teuchos::RCP<FEApp::SGTangentOp> op = 
    Teuchos::rcp(new FEApp::SGTangentOp(sg_expansion,
					alpha, beta, sum_derivs,
					sg_overlapped_xdot, 
					sg_overlapped_x,
					deriv_p,
					overlapped_Vx,
					overlapped_Vx,
					vp,
					sg_overlapped_ff, 
					sg_overlapped_JVx,
					sg_overlapped_fVp));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGTangentType> > pde = 
    pdeTM.getAsObject<FEApp::SGTangentType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::SGTangentType> globalFill(disc->getMesh(), 
						     quad, pde, bc, 
						     sg_xdot!=NULL);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (sg_f != NULL)
    for (int i=0; i<sg_f->size(); i++)
      (*sg_f)[i].Export((*sg_overlapped_f)[i], *exporter, Add);

  // Assemble derivatives
  if (sg_JVx != NULL)
    for (int i=0; i<sg_JVx->size(); i++)
      (*sg_JVx)[i].Export((*sg_overlapped_JVx)[i], *exporter, Add);
  if (sg_fVp != NULL)
    for (int i=0; i<sg_fVp->size(); i++)
      (*sg_fVp)[i].Export((*sg_overlapped_fVp)[i], *exporter, Add);
}

void
FEApp::Application::
evaluateSGResponses(const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
		    const Stokhos::EpetraVectorOrthogPoly& sg_x,
		    const Teuchos::Array<ParamVec>& p,
		    const Teuchos::Array<int>& sg_p_index,
		    const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
		    Stokhos::EpetraVectorOrthogPoly& sg_g)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateSGResponses");

  const Epetra_Comm& comm = sg_x[0].Map().Comm();
  unsigned int offset = 0;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    sg_x.basis();
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Teuchos::RCP<Epetra_LocalMap> local_response_map = 
      Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(num_responses), 0, comm));

    // Create Epetra_Vector for response function
    Stokhos::EpetraVectorOrthogPoly local_sg_g(
      basis, sg_overlap_map, local_response_map,
      sg_x.productComm());

    // Evaluate response function
    responses[i]->evaluateSGResponses(sg_xdot, sg_x, p, sg_p_index, sg_p_vals, 
				      local_sg_g);

    // Copy result into combined result
    for (int k=0; k<sg_g.size(); k++)
      for (unsigned int j=0; j<num_responses; j++)
	sg_g[k][offset+j] = local_sg_g[k][j];

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::
evaluateSGResponseTangents(
      const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
      const Stokhos::EpetraVectorOrthogPoly& sg_x,
      const Teuchos::Array<ParamVec>& p,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
      const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
      const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
      Stokhos::EpetraVectorOrthogPoly* sg_g,
      const Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > >& sg_gt)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateSGResponseTangents");

  const Epetra_Comm& comm = sg_x[0].Map().Comm();
  unsigned int offset = 0;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    sg_x.basis();
  Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > > local_sg_gt(sg_gt.size());
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Teuchos::RCP<Epetra_LocalMap> local_response_map = 
      Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(num_responses), 0, comm));

    // Create Epetra_Vectors for response function
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly > local_sg_g;
    if (sg_g != NULL)
      local_sg_g = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, local_response_map,
		       sg_x.productComm()));
    for (int j=0; j<sg_gt.size(); j++)
      if (sg_gt[j] != Teuchos::null)
	local_sg_gt[j] = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			 basis, sg_overlap_map, local_response_map,
			 sg_x.productComm(), 
			 (*sg_gt[j])[0].NumVectors()));

    // Evaluate response function
    responses[i]->evaluateSGTangents(sg_xdot, sg_x, p, sg_p_index, sg_p_vals, 
				     deriv_p, dxdot_dp, dx_dp, 
				     local_sg_g.get(), local_sg_gt);

    // Copy results into combined result
    for (unsigned int j=0; j<num_responses; j++) {
      if (sg_g != NULL)
	for (int k=0; k<sg_g->size(); k++)
	  (*sg_g)[k][offset+j] = (*local_sg_g)[k][j];
      for (int l=0; l<sg_gt.size(); l++)
	if (sg_gt[l] != Teuchos::null)
	  for (int m=0; m<sg_gt[l]->size(); m++)
	    for (int k=0; k<(*sg_gt[l])[m].NumVectors(); k++)
	      (*sg_gt[l])[m][k][offset+j] = (*local_sg_gt[l])[m][k][j];
    }

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::
evaluateSGResponseGradients(
  const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
  const Stokhos::EpetraVectorOrthogPoly& sg_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& sg_p_index,
  const Teuchos::Array< Teuchos::Array<SGType> >& sg_p_vals,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  Stokhos::EpetraVectorOrthogPoly* sg_g,
  Stokhos::EpetraMultiVectorOrthogPoly* sg_dg_dx,
  Stokhos::EpetraMultiVectorOrthogPoly* sg_dg_dxdot,
  const Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > >& sg_dg_dp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateSGResponseGradients");

  const Epetra_Comm& comm = sg_x[0].Map().Comm();
  unsigned int offset = 0;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    sg_x.basis();
  Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > > local_sg_dgdp(sg_dg_dp.size());
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Teuchos::RCP<Epetra_LocalMap> local_response_map = 
      Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(num_responses), 0, comm));

    // Create Epetra_Vectors for response function
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > local_sg_g;
    if (sg_g != NULL)
      local_sg_g = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, local_response_map,
		       sg_x.productComm()));
    Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > local_sg_dgdx;
    if (sg_dg_dx != NULL)
      local_sg_dgdx = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		       basis, sg_overlap_map, disc->getMap(), 
		       sg_x.productComm(), num_responses));
    Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > local_sg_dgdxdot;
    if (sg_dg_dxdot != NULL)
      local_sg_dgdxdot = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		       basis, sg_overlap_map, disc->getMap(), 
		       sg_x.productComm(), num_responses));

    for (int j=0; j<sg_dg_dp.size(); j++)
      if (sg_dg_dp[j] != Teuchos::null)
	local_sg_dgdp[j] = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			 basis, sg_overlap_map, local_response_map,
			 sg_x.productComm(),
			 (*sg_dg_dp[j])[0].NumVectors()));

    // Evaluate response function
    responses[i]->evaluateSGGradients(sg_xdot, sg_x, p, sg_p_index, sg_p_vals, 
				      deriv_p,
				      local_sg_g.get(), 
				      local_sg_dgdx.get(), 
				      local_sg_dgdxdot.get(), 
				      local_sg_dgdp);

    // Copy results into combined result
    for (unsigned int j=0; j<num_responses; j++) {
      if (sg_g != NULL)
	for (int m=0; m<sg_g->size(); m++)
	  (*sg_g)[m][offset+j] = (*local_sg_g)[m][j];
      if (sg_dg_dx != NULL)
	for (int m=0; m<sg_dg_dx->size(); m++)
	  (*sg_dg_dx)[m](offset+j)->Update(1.0, *((*local_sg_dgdx)[m](j)), 0.0);
      if (sg_dg_dxdot != NULL)
	for (int m=0; m<sg_dg_dxdot->size(); m++)
	  (*sg_dg_dxdot)[m](offset+j)->Update(1.0, *((*local_sg_dgdxdot)[m](j)),
					      0.0);
      for (int l=0; l<sg_dg_dp.size(); l++)
	if (sg_dg_dp[l] != Teuchos::null)
	  for (int m=0; m<sg_dg_dp[l]->size(); m++)
	    for (int k=0; k<(*sg_dg_dp[l])[m].NumVectors(); k++)
	      (*sg_dg_dp[l])[m][k][offset+j] = (*local_sg_dgdp[l])[m][k][j];
    }

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::computeGlobalMPResidual(
  const Stokhos::ProductEpetraVector* mp_xdot,
  const Stokhos::ProductEpetraVector& mp_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& mp_p_index,
  const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
  Stokhos::ProductEpetraVector& mp_f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalMPResidual");

  // Create overlapped multi-point Epetra objects
  if (mp_overlapped_x == Teuchos::null || 
      mp_overlapped_x->size() != mp_x.size())
    mp_overlapped_x = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));
  if (mp_xdot && (mp_overlapped_xdot == Teuchos::null || 
		  mp_overlapped_xdot->size() != mp_x.size()))
    mp_overlapped_xdot = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));

  if (mp_overlapped_f == Teuchos::null || 
      mp_overlapped_f->size() != mp_x.size())
    mp_overlapped_f = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));

  for (int i=0; i<mp_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*mp_overlapped_x)[i].Import(mp_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (mp_xdot)
      (*mp_overlapped_xdot)[i].Import((*mp_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    (*mp_overlapped_f)[i].PutScalar(0.0);
    mp_f[i].PutScalar(0.0);

  }

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Set MP parameters
  for (int i=0; i<mp_p_index.size(); i++) {
    int ii = mp_p_index[i];
    for (unsigned int j=0; j<p[ii].size(); j++)
	p[ii][j].family->setValue<FEApp::MPResidualType>(mp_p_vals[ii][j]);
  }

  // Create residual init/post op
  Teuchos::RCP<FEApp::MPResidualOp> mp_res_fill_op = 
    Teuchos::rcp(new FEApp::MPResidualOp(mp_overlapped_xdot, 
					 mp_overlapped_x, 
					 mp_overlapped_f));
    
  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::MPResidualType> > pde = 
    pdeTM.getAsObject<FEApp::MPResidualType>();

  // Instantiate global fill
  if (mp_res_global_fill == Teuchos::null) {
    mp_res_global_fill = 
      Teuchos::rcp(new FEApp::GlobalFill<FEApp::MPResidualType>(
		     disc->getMesh(), quad, pde, bc, mp_xdot!=NULL));
  }

  // Do global fill
  mp_res_global_fill->computeGlobalFill(*mp_res_fill_op);

  // Assemble global residual
  for (int i=0; i<mp_f.size(); i++)
    mp_f[i].Export((*mp_overlapped_f)[i], *exporter, Add);
}

void
FEApp::Application::computeGlobalMPJacobian(
  double alpha, double beta,
  const Stokhos::ProductEpetraVector* mp_xdot,
  const Stokhos::ProductEpetraVector& mp_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& mp_p_index,
  const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
  Stokhos::ProductEpetraVector* mp_f,
  Stokhos::ProductEpetraOperator& mp_jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalMPJacobian");

  // Create overlapped multi-point Epetra objects
  if (mp_overlapped_x == Teuchos::null || 
      mp_overlapped_x->size() != mp_x.size())
    mp_overlapped_x = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));
  if (mp_xdot && (mp_overlapped_xdot == Teuchos::null || 
		  mp_overlapped_xdot->size() != mp_x.size()))
    mp_overlapped_xdot = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));

  if (mp_overlapped_f == Teuchos::null || 
      mp_overlapped_f->size() != mp_x.size())
    mp_overlapped_f = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));

  if (mp_overlapped_jac == Teuchos::null || 
      mp_overlapped_jac->size() != mp_x.size())
    mp_overlapped_jac = 
      Teuchos::rcp(new Stokhos::ProductContainer<Epetra_CrsMatrix>(
		     mp_x.map(), 
		     *overlapped_jac));

  for (int i=0; i<mp_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*mp_overlapped_x)[i].Import(mp_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (mp_xdot)
      (*mp_overlapped_xdot)[i].Import((*mp_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    if (mp_f != NULL) {
      (*mp_overlapped_f)[i].PutScalar(0.0);
      (*mp_f)[i].PutScalar(0.0);
    }

    // Zero out overlapped Jacobian
    (*mp_overlapped_jac)[i].PutScalar(0.0);

  } 

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Set MP parameters
  for (int i=0; i<mp_p_index.size(); i++) {
    int ii = mp_p_index[i];
    for (unsigned int j=0; j<p[ii].size(); j++)
	p[ii][j].family->setValue<FEApp::MPJacobianType>(mp_p_vals[ii][j]);
  }

  // Create Jacobian init/post op
  Teuchos::RCP< Stokhos::ProductEpetraVector > mp_overlapped_ff;
  if (mp_f != NULL)
    mp_overlapped_ff = mp_overlapped_f;
  Teuchos::RCP<FEApp::MPJacobianOp> mp_jac_fill_op = 
    Teuchos::rcp(new FEApp::MPJacobianOp(alpha, beta, 
					 mp_overlapped_xdot, 
					 mp_overlapped_x, 
					 mp_overlapped_ff, 
					 mp_overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::MPJacobianType> > pde = 
    pdeTM.getAsObject<FEApp::MPJacobianType>();

  // Instantiate global fill
  if (mp_jac_global_fill == Teuchos::null) {
    mp_jac_global_fill = 
      Teuchos::rcp(new FEApp::GlobalFill<FEApp::MPJacobianType>(
		     disc->getMesh(), quad, pde, bc,mp_xdot!=NULL));
  }

  // Do global fill
  mp_jac_global_fill->computeGlobalFill(*mp_jac_fill_op);
  
  // Assemble global residual
  if (mp_f != NULL)
    for (int i=0; i<mp_f->size(); i++)
      (*mp_f)[i].Export((*mp_overlapped_f)[i], *exporter, Add);
    
  // Assemble block Jacobians
  Teuchos::RCP<Epetra_CrsMatrix> jac;
  for (int i=0; i<mp_jac.size(); i++) {
    jac = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(mp_jac.getCoeffPtr(i), 
						      true);
    jac->PutScalar(0.0);
    jac->Export((*mp_overlapped_jac)[i], *exporter, Add);
    jac->FillComplete(true);
  }
}

void
FEApp::Application::computeGlobalMPTangent(
  double alpha, double beta, bool sum_derivs,
  const Stokhos::ProductEpetraVector* mp_xdot,
  const Stokhos::ProductEpetraVector& mp_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& mp_p_index,
  const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
  const Teuchos::RCP<ParamVec>& deriv_p,
  const Epetra_MultiVector* Vx,
  const Teuchos::SerialDenseMatrix<int,double>* Vp,
  Stokhos::ProductEpetraVector* mp_f,
  Stokhos::ProductEpetraMultiVector* mp_JVx,
  Stokhos::ProductEpetraMultiVector* mp_fVp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalMPTangent");

  // Create overlapped multi-point Epetra objects
  if (mp_overlapped_x == Teuchos::null || 
      mp_overlapped_x->size() != mp_x.size())
    mp_overlapped_x = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));
  if (mp_xdot && (mp_overlapped_xdot == Teuchos::null || 
		  mp_overlapped_xdot->size() != mp_x.size()))
    mp_overlapped_xdot = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));

  if (mp_overlapped_f == Teuchos::null || 
      mp_overlapped_f->size() != mp_x.size())
    mp_overlapped_f = 
      Teuchos::rcp(new Stokhos::ProductEpetraVector(mp_x.map(),
						    disc->getOverlapMap(),
						    mp_x.productComm()));

  for (int i=0; i<mp_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*mp_overlapped_x)[i].Import(mp_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (mp_xdot)
      (*mp_overlapped_xdot)[i].Import((*mp_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    if (mp_f != NULL) {
      (*mp_overlapped_f)[i].PutScalar(0.0);
      (*mp_f)[i].PutScalar(0.0);
    }

  }

  // Set parameters
  for (int i=0; i<p.size(); i++)
    for (unsigned int j=0; j<p[i].size(); j++)
      p[i][j].family->setRealValueForAllTypes(p[i][j].baseValue);

  // Set MP parameters
  for (int i=0; i<mp_p_index.size(); i++) {
    int ii = mp_p_index[i];
    for (unsigned int j=0; j<p[ii].size(); j++)
	p[ii][j].family->setValue<FEApp::MPTangentType>(mp_p_vals[ii][j]);
  }

  Teuchos::RCP< Stokhos::ProductEpetraMultiVector > mp_overlapped_JVx;
  if (mp_JVx != NULL) {
    mp_overlapped_JVx = 
      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
		     mp_x.map(), disc->getOverlapMap(), mp_x.productComm(), 
		     (*mp_JVx)[0].NumVectors()));
    mp_JVx->init(0.0);
  }
  
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector > mp_overlapped_fVp;
  if (mp_fVp != NULL) {
    mp_overlapped_fVp = 
      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
		     mp_x.map(), disc->getOverlapMap(), mp_x.productComm(),
		     (*mp_fVp)[0].NumVectors()));
    mp_fVp->init(0.0);
  }

  Teuchos::RCP<Epetra_MultiVector> overlapped_Vx;
  if (Vx != NULL) {
    overlapped_Vx = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          Vx->NumVectors()));
  }

  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> > vp =
    Teuchos::rcp(Vp, false);

  // Create Jacobian init/post op
   Teuchos::RCP< Stokhos::ProductEpetraVector > mp_overlapped_ff;
  if (mp_f != NULL)
    mp_overlapped_ff = mp_overlapped_f;
  Teuchos::RCP<FEApp::MPTangentOp> op = 
    Teuchos::rcp(new FEApp::MPTangentOp(alpha, beta, sum_derivs,
					mp_overlapped_xdot, 
					mp_overlapped_x,
					deriv_p,
					overlapped_Vx,
					overlapped_Vx,
					vp,
					mp_overlapped_ff, 
					mp_overlapped_JVx,
					mp_overlapped_fVp));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::MPTangentType> > pde = 
    pdeTM.getAsObject<FEApp::MPTangentType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::MPTangentType> globalFill(disc->getMesh(), 
						     quad, pde, bc, 
						     mp_xdot!=NULL);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (mp_f != NULL)
    for (int i=0; i<mp_f->size(); i++)
      (*mp_f)[i].Export((*mp_overlapped_f)[i], *exporter, Add);

  // Assemble derivatives
  if (mp_JVx != NULL)
    for (int i=0; i<mp_JVx->size(); i++)
      (*mp_JVx)[i].Export((*mp_overlapped_JVx)[i], *exporter, Add);
  if (mp_fVp != NULL)
    for (int i=0; i<mp_fVp->size(); i++)
      (*mp_fVp)[i].Export((*mp_overlapped_fVp)[i], *exporter, Add);
}

void
FEApp::Application::
evaluateMPResponses(
  const Stokhos::ProductEpetraVector* mp_xdot,
  const Stokhos::ProductEpetraVector& mp_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& mp_p_index,
  const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
  Stokhos::ProductEpetraVector& mp_g)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateMPResponses");
  const Epetra_Vector* xdot = NULL;

  // Create a copy of p (since we shouldn't modify p below)
  Teuchos::Array<ParamVec> p2(p);

  for (int i=0; i<mp_x.size(); i++) {
    if (mp_xdot != NULL)
      xdot = mp_xdot->getCoeffPtr(i).get();

    // Set the base value for each MP parameter
    for (int j=0; j<mp_p_index.size(); j++) {
      int jj = mp_p_index[j];
      for (unsigned int l=0; l<p2[jj].size(); l++)
	p2[jj][l].baseValue = mp_p_vals[jj][l].coeff(i);
    }

    evaluateResponses(xdot, mp_x[i], p2, mp_g[i]);
  }
}

void
FEApp::Application::
evaluateMPResponseTangents(
  const Stokhos::ProductEpetraVector* mp_xdot,
  const Stokhos::ProductEpetraVector& mp_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& mp_p_index,
  const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
  Stokhos::ProductEpetraVector* mp_g,
  const Teuchos::Array< Teuchos::RCP< Stokhos::ProductEpetraMultiVector > >& mp_gt)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateMPResponseTangents");
  const Epetra_Vector* xdot = NULL;
  Epetra_Vector* g = NULL;
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > gt(mp_gt.size());

  // Create a copy of p (since we shouldn't modify p below)
  Teuchos::Array<ParamVec> p2(p);

  // Create a copy of deriv_p (since we shouldn't modify p below)
  Teuchos::Array< Teuchos::RCP<ParamVec> > deriv_p2(deriv_p.size());
  for (int i=0; i<deriv_p.size(); i++)
    if (deriv_p[i] != Teuchos::null)
      deriv_p2[i] = Teuchos::rcp(new ParamVec(*(deriv_p[i])));

  for (int i=0; i<mp_x.size(); i++) {
    if (mp_xdot != NULL)
      xdot = mp_xdot->getCoeffPtr(i).get();
    if (mp_g != NULL)
      g = mp_g->getCoeffPtr(i).get();
    for (int j=0; j<mp_gt.size(); j++)
      if (mp_gt[j] != Teuchos::null)
	gt[j] = mp_gt[j]->getCoeffPtr(i);

    // Set the base value for each MP parameter
    for (int j=0; j<mp_p_index.size(); j++) {
      int jj = mp_p_index[j];
      for (unsigned int l=0; l<p2[jj].size(); l++) {
	p2[jj][l].baseValue = mp_p_vals[jj][l].coeff(i);
	if (deriv_p[jj] != Teuchos::null) {
	  for (unsigned int k=0; k<deriv_p2[jj]->size(); k++) {
	    if ((*(deriv_p2[jj]))[k].family->getName() ==
		p2[jj][l].family->getName())
	      (*(deriv_p2[jj]))[k].baseValue = mp_p_vals[jj][l].coeff(i);
	  }
	}
      }
    }
    
    evaluateResponseTangents(xdot, mp_x[i], p2, deriv_p2, dxdot_dp, dx_dp, g, 
			     gt);
  }
}

void
FEApp::Application::
evaluateMPResponseGradients(
  const Stokhos::ProductEpetraVector* mp_xdot,
  const Stokhos::ProductEpetraVector& mp_x,
  const Teuchos::Array<ParamVec>& p,
  const Teuchos::Array<int>& mp_p_index,
  const Teuchos::Array< Teuchos::Array<MPType> >& mp_p_vals,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  Stokhos::ProductEpetraVector* mp_g,
  Stokhos::ProductEpetraMultiVector* mp_dg_dx,
  Stokhos::ProductEpetraMultiVector* mp_dg_dxdot,
  const Teuchos::Array< Teuchos::RCP< Stokhos::ProductEpetraMultiVector > >& mp_dg_dp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateMPResponseGradients");

  const Epetra_Vector* xdot = NULL;
  Epetra_Vector* g = NULL;
  Epetra_MultiVector* dg_dx = NULL;
  Epetra_MultiVector* dg_dxdot = NULL;
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > dg_dp(mp_dg_dp.size());

  // Create a copy of p (since we shouldn't modify p below)
  Teuchos::Array<ParamVec> p2(p);

  // Create a copy of deriv_p (since we shouldn't modify p below)
  Teuchos::Array< Teuchos::RCP<ParamVec> > deriv_p2(deriv_p.size());
  for (int i=0; i<deriv_p.size(); i++)
    if (deriv_p[i] != Teuchos::null)
      deriv_p2[i] = Teuchos::rcp(new ParamVec(*(deriv_p[i])));

  for (int i=0; i<mp_x.size(); i++) {
    if (mp_xdot != NULL)
      xdot = mp_xdot->getCoeffPtr(i).get();
    if (mp_g != NULL)
      g = mp_g->getCoeffPtr(i).get();
    if (mp_dg_dx != NULL)
      dg_dx = mp_dg_dx->getCoeffPtr(i).get();
    if (mp_dg_dxdot != NULL)
      dg_dxdot = mp_dg_dxdot->getCoeffPtr(i).get();
    for (int j=0; j<mp_dg_dp.size(); j++)
      if (mp_dg_dp[j] != Teuchos::null)
	dg_dp[j] = mp_dg_dp[j]->getCoeffPtr(i);

    // Set the base value for each MP parameter
    for (int j=0; j<mp_p_index.size(); j++) {
      int jj = mp_p_index[j];
      for (unsigned int l=0; l<p2[jj].size(); l++) {
	p2[jj][l].baseValue = mp_p_vals[jj][l].coeff(i);
	if (deriv_p[jj] != Teuchos::null) {
	  for (unsigned int k=0; k<deriv_p2[jj]->size(); k++) {
	    if ((*(deriv_p2[jj]))[k].family->getName() ==
		p2[jj][l].family->getName())
	      (*(deriv_p2[jj]))[k].baseValue = mp_p_vals[jj][l].coeff(i);
	  }
	}
      }
    }

    evaluateResponseGradients(xdot, mp_x[i], p2, deriv_p2, g, dg_dx, dg_dxdot,
			      dg_dp);
  }
}
