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

#include "FEApp_ModelEvaluator.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Assert.hpp"

FEApp::ModelEvaluator::ModelEvaluator(
  const Teuchos::RCP<FEApp::Application>& app_,
  const Teuchos::RCP<Teuchos::ParameterList>& appParams
) 
  : app(app_),
    eval_W_with_f(false)
{
  Teuchos::ParameterList& problemParams = appParams->sublist("Problem");

  // Parameters (e.g., for sensitivities, SG expansions, ...)
  Teuchos::ParameterList& parameterParams = 
    problemParams.sublist("Parameters");
  int num_param_vecs = 
    parameterParams.get("Number of Parameter Vectors", 0);
  param_names.resize(num_param_vecs);
  for (int i=0; i<num_param_vecs; i++) {
    std::stringstream ss;
    ss << "Parameter Vector " << i;
    Teuchos::ParameterList& pList = parameterParams.sublist(ss.str());
    int numParameters = pList.get<int>("Number");
    TEUCHOS_TEST_FOR_EXCEPTION(
      numParameters == 0, 
      Teuchos::Exceptions::InvalidParameter,
      std::endl << "Error!  FEApp::ModelEvaluator::ModelEvaluator():  " <<
      "Parameter vector " << i << " has zero parameters!" << std::endl);
    param_names[i] = 
      Teuchos::rcp(new Teuchos::Array<std::string>(numParameters));
    for (int j=0; j<numParameters; j++) {
      std::stringstream ss2;
      ss2 << "Parameter " << j;
      (*param_names[i])[j] = pList.get<std::string>(ss2.str());
    }
  }

  // Setup sacado and epetra storage for parameters  
  sacado_param_vec.resize(num_param_vecs);
  epetra_param_map.resize(num_param_vecs);
  epetra_param_vec.resize(num_param_vecs);
  p_sg_vals.resize(num_param_vecs);
  p_mp_vals.resize(num_param_vecs);
  const Epetra_Comm& comm = app->getMap()->Comm();
  for (int i=0; i<num_param_vecs; i++) {

    // Initialize Sacado parameter vector
    app->getParamLib()->fillVector<FEApp::ResidualType>(*(param_names[i]), 
							sacado_param_vec[i]);

    // Create Epetra map for parameter vector
    epetra_param_map[i] = 
      Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(sacado_param_vec[i].size()), 0, comm));

    // Create Epetra vector for parameters
    epetra_param_vec[i] = 
      Teuchos::rcp(new Epetra_Vector(*(epetra_param_map[i])));
    for (unsigned int j=0; j<sacado_param_vec[i].size(); j++)
      (*(epetra_param_vec[i]))[j] = sacado_param_vec[i][j].baseValue;

    p_sg_vals[i].resize(sacado_param_vec[i].size());
    p_mp_vals[i].resize(sacado_param_vec[i].size());
  }

  supports_g = (app->getResponseMap() != Teuchos::null);
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_x_map() const
{
  return app->getMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_f_map() const
{
  return app->getMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= static_cast<int>(epetra_param_map.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return epetra_param_map[l];
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_g_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(supports_g == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_g_map():  " <<
                     "No response functions have been supplied.  " <<
                     "Supplied index l = " << l << std::endl);
  TEUCHOS_TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_g_map() only " <<
                     " supports 1 response vector.  Supplied index l = " << 
                     l << std::endl);

  return app->getResponseMap();
}

Teuchos::RCP<const Teuchos::Array<std::string> >
FEApp::ModelEvaluator::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= static_cast<int>(param_names.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return param_names[l];
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_x_init() const
{
  return app->getInitialSolution();
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= static_cast<int>(param_names.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return epetra_param_vec[l];
}

Teuchos::RCP<Epetra_Operator>
FEApp::ModelEvaluator::create_W() const
{
  my_W = app->createW();
  return my_W;
}

EpetraExt::ModelEvaluator::InArgs
FEApp::ModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_t,true);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);
  inArgs.set_Np(param_names.size());
  
  inArgs.setSupports(IN_ARG_x_sg,true);
  inArgs.setSupports(IN_ARG_x_dot_sg,true);
  for (int i=0; i<param_names.size(); i++)
    inArgs.setSupports(IN_ARG_p_sg, i, true);
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);
  
  inArgs.setSupports(IN_ARG_x_mp,true);
  inArgs.setSupports(IN_ARG_x_dot_mp,true);
  for (int i=0; i<param_names.size(); i++)
    inArgs.setSupports(IN_ARG_p_mp, i, true); 
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
FEApp::ModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  int n_g = 0;
  if (supports_g)
    n_g = 1;

  // Deterministic
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(DERIV_LINEARITY_UNKNOWN ,DERIV_RANK_FULL ,true));
  outArgs.set_Np_Ng(param_names.size(), n_g);

  for (int i=0; i<param_names.size(); i++)
    outArgs.setSupports(OUT_ARG_DfDp, i, DerivativeSupport(DERIV_MV_BY_COL));
  if (supports_g) {
    outArgs.setSupports(OUT_ARG_DgDx, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    outArgs.setSupports(OUT_ARG_DgDx_dot, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    for (int i=0; i<param_names.size(); i++)
      outArgs.setSupports(OUT_ARG_DgDp, 0, i, 
			  DerivativeSupport(DERIV_MV_BY_COL));
  }
  
  // Stochastic
  outArgs.setSupports(OUT_ARG_f_sg,true);
  outArgs.setSupports(OUT_ARG_W_sg,true);
  for (int i=0; i<param_names.size(); i++)
    outArgs.setSupports(OUT_ARG_DfDp_sg, i, DerivativeSupport(DERIV_MV_BY_COL));
  if (supports_g) {
    outArgs.setSupports(OUT_ARG_g_sg, 0, true);
    outArgs.setSupports(OUT_ARG_DgDx_sg, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    outArgs.setSupports(OUT_ARG_DgDx_dot_sg, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    for (int i=0; i<param_names.size(); i++)
      outArgs.setSupports(OUT_ARG_DgDp_sg, 0, i, 
			  DerivativeSupport(DERIV_MV_BY_COL));
  }
      
  // Multi-point
  outArgs.setSupports(OUT_ARG_f_mp,true);
  outArgs.setSupports(OUT_ARG_W_mp,true);
  for (int i=0; i<param_names.size(); i++)
    outArgs.setSupports(OUT_ARG_DfDp_mp, i, DerivativeSupport(DERIV_MV_BY_COL));
  if (supports_g) {
    outArgs.setSupports(OUT_ARG_g_mp, 0, true);
    outArgs.setSupports(OUT_ARG_DgDx_mp, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    outArgs.setSupports(OUT_ARG_DgDx_dot_mp, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    for (int i=0; i<param_names.size(); i++)
      outArgs.setSupports(OUT_ARG_DgDp_mp, 0, i, 
			  DerivativeSupport(DERIV_MV_BY_COL));
  }

  return outArgs;
}

void 
FEApp::ModelEvaluator::evalModel(const InArgs& inArgs, 
				 const OutArgs& outArgs) const
{
  //
  // Get the input arguments
  //
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  Teuchos::RCP<const Epetra_Vector> x_dot = inArgs.get_x_dot();
  double alpha = 0.0;
  double beta = 1.0;
  if (x_dot != Teuchos::null) {
    alpha = inArgs.get_alpha();
    beta = inArgs.get_beta();
  }
  for (int i=0; i<inArgs.Np(); i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i);
    if (p != Teuchos::null) {
      for (unsigned int j=0; j<sacado_param_vec[i].size(); j++)
	sacado_param_vec[i][j].baseValue = (*p)[j];
    }
  }

  //
  // Get the output arguments
  //
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  
  //
  // Compute the functions
  //
  bool f_computed = false;

  // W matrix
  if (f_out != Teuchos::null && eval_W_with_f) {
    app->computeGlobalJacobian(alpha, beta, x_dot.get(), *x, sacado_param_vec,
			       f_out.get(), *my_W);
    f_computed = true;
  }
  else if (W_out != Teuchos::null && !eval_W_with_f) {
    if (f_out.getType() == EVAL_TYPE_EXACT ||
        f_out.getType() == EVAL_TYPE_APPROX_DERIV)
      app->computeGlobalJacobian(alpha, beta, x_dot.get(), *x, sacado_param_vec,
                                 f_out.get(), *W_out);
    else
      app->computeGlobalPreconditioner(alpha, beta, x_dot.get(), *x, 
                                       sacado_param_vec, f_out.get(), *W_out);
    f_computed = true;
  }
  
  // df/dp
  for (int i=0; i<outArgs.Np(); i++) {
    Teuchos::RCP<Epetra_MultiVector> dfdp_out = 
      outArgs.get_DfDp(i).getMultiVector();
    if (dfdp_out != Teuchos::null) {
      Teuchos::Array<int> p_indexes = 
	outArgs.get_DfDp(i).getDerivativeMultiVector().getParamIndexes();
      Teuchos::RCP<ParamVec> p_vec;
      if (p_indexes.size() == 0)
	p_vec = Teuchos::rcp(&sacado_param_vec[i],false);
      else {
	p_vec = Teuchos::rcp(new ParamVec);
	for (int j=0; j<p_indexes.size(); j++)
	  p_vec->addParam(sacado_param_vec[i][p_indexes[j]].family, 
			  sacado_param_vec[i][p_indexes[j]].baseValue);
      }
      
      app->computeGlobalTangent(0.0, 0.0, false, x_dot.get(), *x, 
				sacado_param_vec, p_vec,
				NULL, NULL, f_out.get(), NULL, 
				dfdp_out.get());

      f_computed = true;
    }
  }

  // f
  if (f_out != Teuchos::null && !f_computed) {
    app->computeGlobalResidual(x_dot.get(), *x, sacado_param_vec, *f_out);
  }

  // Response functions
  if (outArgs.Ng() > 0) {
    Teuchos::RCP<Epetra_Vector> g_out = outArgs.get_g(0);
    Teuchos::RCP<Epetra_MultiVector> dgdx_out = 
      outArgs.get_DgDx(0).getMultiVector();
    Teuchos::RCP<Epetra_MultiVector> dgdxdot_out = 
      outArgs.get_DgDx_dot(0).getMultiVector();
    
    Teuchos::Array< Teuchos::RCP<ParamVec> > p_vec(outArgs.Np());
    Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > dgdp_out(outArgs.Np());
    bool have_dgdp = false;
    for (int i=0; i<outArgs.Np(); i++) {
      dgdp_out[i] = outArgs.get_DgDp(0,i).getMultiVector();
      if (dgdp_out[i] != Teuchos::null) {
	have_dgdp = true;
	Teuchos::Array<int> p_indexes = 
	  outArgs.get_DgDp(0,i).getDerivativeMultiVector().getParamIndexes();
	if (p_indexes.size() == 0)
	  p_vec[i] = Teuchos::rcp(&sacado_param_vec[i],false);
	else {
	  p_vec[i] = Teuchos::rcp(new ParamVec);
	  for (int j=0; j<p_indexes.size(); j++)
	    p_vec[i]->addParam(sacado_param_vec[i][p_indexes[j]].family, 
			       sacado_param_vec[i][p_indexes[j]].baseValue);
	}
      }
    }

    if (have_dgdp ||dgdx_out != Teuchos::null || dgdxdot_out != Teuchos::null) {
      app->evaluateResponseGradients(x_dot.get(), *x, sacado_param_vec, 
				     p_vec, g_out.get(), dgdx_out.get(), 
                                     dgdxdot_out.get(), dgdp_out);
    }
    else if (g_out != Teuchos::null)
      app->evaluateResponses(x_dot.get(), *x, sacado_param_vec, *g_out);
  }

  // SG calculation
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();
  if (x_sg != Teuchos::null) {
    app->init_sg(inArgs.get_sg_basis(), inArgs.get_sg_quadrature(), 
		 inArgs.get_sg_expansion(), x_sg->productComm());
    InArgs::sg_const_vector_t x_dot_sg = inArgs.get_x_dot_sg();
    Teuchos::Array<int> p_sg_index;
    for (int i=0; i<inArgs.Np(); i++) {
      InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(i);
      if (p_sg != Teuchos::null) {
	p_sg_index.push_back(i);
	for (int j=0; j<p_sg_vals[i].size(); j++) {
	  int num_sg_blocks = p_sg->size();
	  p_sg_vals[i][j].reset(app->getStochasticExpansion(), num_sg_blocks);
	  p_sg_vals[i][j].copyForWrite();
	  for (int l=0; l<num_sg_blocks; l++) {
	    p_sg_vals[i][j].fastAccessCoeff(l) = (*p_sg)[l][j];
	  }
	}
      }
    }

    OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
    OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
    bool f_sg_computed = false;
      
    // W_sg
    if (W_sg != Teuchos::null) {
      app->computeGlobalSGJacobian(alpha, beta, x_dot_sg.get(), *x_sg, 
				   sacado_param_vec, p_sg_index, p_sg_vals, 
				   f_sg.get(), *W_sg);
      f_sg_computed = true;
    }

    // df/dp_sg
    for (int i=0; i<outArgs.Np(); i++) {
      Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dfdp_sg 
	= outArgs.get_DfDp_sg(i).getMultiVector();
      if (dfdp_sg != Teuchos::null) {
	Teuchos::Array<int> p_indexes = 
	  outArgs.get_DfDp_sg(i).getDerivativeMultiVector().getParamIndexes();
	Teuchos::RCP<ParamVec> p_vec;
	if (p_indexes.size() == 0)
	  p_vec = Teuchos::rcp(&sacado_param_vec[i],false);
	else {
	  p_vec = Teuchos::rcp(new ParamVec);
	  for (int j=0; j<p_indexes.size(); j++)
	    p_vec->addParam(sacado_param_vec[i][p_indexes[j]].family, 
			    sacado_param_vec[i][p_indexes[j]].baseValue);
	}
	
	app->computeGlobalSGTangent(0.0, 0.0, false, x_dot_sg.get(), *x_sg, 
				    sacado_param_vec, p_sg_index, p_sg_vals, 
				    p_vec, NULL, NULL, 
				    f_sg.get(), NULL, dfdp_sg.get());
	
	f_sg_computed = true;
      }
    }

    if (f_sg != Teuchos::null && !f_sg_computed)
      app->computeGlobalSGResidual(x_dot_sg.get(), *x_sg, sacado_param_vec, 
				   p_sg_index, p_sg_vals, *f_sg);

    // Response functions
    if (outArgs.Ng() > 0) {
      Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > g_sg 
	= outArgs.get_g_sg(0);
      Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdx_sg 
	= outArgs.get_DgDx_sg(0).getMultiVector();
      Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdxdot_sg = 
	outArgs.get_DgDx_dot_sg(0).getMultiVector();
    
      Teuchos::Array< Teuchos::RCP<ParamVec> > p_vec(outArgs.Np());
      Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > > dgdp_sg(outArgs.Np());
      bool have_dgdp = false;
      for (int i=0; i<outArgs.Np(); i++) {
	dgdp_sg[i] = outArgs.get_DgDp_sg(0,i).getMultiVector();
	if (dgdp_sg[i] != Teuchos::null) {
	  have_dgdp = true;
	  Teuchos::Array<int> p_indexes = 
	    outArgs.get_DgDp_sg(0,i).getDerivativeMultiVector().getParamIndexes();
	  if (p_indexes.size() == 0)
	    p_vec[i] = Teuchos::rcp(&sacado_param_vec[i],false);
	  else {
	    p_vec[i] = Teuchos::rcp(new ParamVec);
	    for (int j=0; j<p_indexes.size(); j++)
	      p_vec[i]->addParam(sacado_param_vec[i][p_indexes[j]].family, 
				 sacado_param_vec[i][p_indexes[j]].baseValue);
	  }
	}
      }

      if (have_dgdp ||dgdx_sg != Teuchos::null || 
	  dgdxdot_sg != Teuchos::null) {
	app->evaluateSGResponseGradients(
	  x_dot_sg.get(), *x_sg, sacado_param_vec, p_sg_index, p_sg_vals, 
	  p_vec, g_sg.get(), dgdx_sg.get(), dgdxdot_sg.get(), dgdp_sg);
      }
      else if (g_sg != Teuchos::null)
	app->evaluateSGResponses(
	  x_dot_sg.get(), *x_sg, sacado_param_vec, p_sg_index, p_sg_vals, 
	  *g_sg);
    }
  }
    
  //
  // Multi-point evaluation
  //
  mp_const_vector_t x_mp = inArgs.get_x_mp();
  if (x_mp != Teuchos::null) {
    mp_const_vector_t x_dot_mp = inArgs.get_x_dot_mp();
    Teuchos::Array<int> p_mp_index;
    for (int i=0; i<inArgs.Np(); i++) {
      mp_const_vector_t p_mp = inArgs.get_p_mp(i);
      if (p_mp != Teuchos::null) {
	p_mp_index.push_back(i);
	for (int j=0; j<p_mp_vals[i].size(); j++) {
	  int num_mp_blocks = p_mp->size();
	  p_mp_vals[i][j].reset(num_mp_blocks);
	  p_mp_vals[i][j].copyForWrite();
	  for (int l=0; l<num_mp_blocks; l++) {
	    p_mp_vals[i][j].fastAccessCoeff(l) = (*p_mp)[l][j];
	  }
	}
      }
    }
    
    mp_vector_t f_mp = outArgs.get_f_mp();
    mp_operator_t W_mp = outArgs.get_W_mp();
    bool f_mp_computed = false;
      
    // W_mp
    if (W_mp != Teuchos::null) {
      app->computeGlobalMPJacobian(alpha, beta, x_dot_mp.get(), *x_mp, 
				   sacado_param_vec, p_mp_index, p_mp_vals, 
				   f_mp.get(), *W_mp);
      f_mp_computed = true;
    }
    
    // df/dp_mp
    for (int i=0; i<outArgs.Np(); i++) {
      Teuchos::RCP< Stokhos::ProductEpetraMultiVector > dfdp_mp 
	= outArgs.get_DfDp_mp(i).getMultiVector();
      if (dfdp_mp != Teuchos::null) {
	Teuchos::Array<int> p_indexes = 
	  outArgs.get_DfDp_mp(i).getDerivativeMultiVector().getParamIndexes();
	Teuchos::RCP<ParamVec> p_vec;
	if (p_indexes.size() == 0)
	  p_vec = Teuchos::rcp(&sacado_param_vec[i],false);
	else {
	  p_vec = Teuchos::rcp(new ParamVec);
	  for (int j=0; j<p_indexes.size(); j++)
	    p_vec->addParam(sacado_param_vec[i][p_indexes[j]].family, 
			    sacado_param_vec[i][p_indexes[j]].baseValue);
	}
	    
	app->computeGlobalMPTangent(0.0, 0.0, false, x_dot_mp.get(), *x_mp, 
				    sacado_param_vec, p_mp_index, p_mp_vals, 
				    p_vec, NULL, NULL, 
				    f_mp.get(), NULL, dfdp_mp.get());
	
	f_mp_computed = true;
      }
    }

    if (f_mp != Teuchos::null && !f_mp_computed)
      app->computeGlobalMPResidual(x_dot_mp.get(), *x_mp, sacado_param_vec, 
				   p_mp_index, p_mp_vals, *f_mp);

    // Response functions
    if (outArgs.Ng() > 0) {
      Teuchos::RCP< Stokhos::ProductEpetraVector > g_mp 
	= outArgs.get_g_mp(0);
      Teuchos::RCP< Stokhos::ProductEpetraMultiVector > dgdx_mp 
	= outArgs.get_DgDx_mp(0).getMultiVector();
      Teuchos::RCP< Stokhos::ProductEpetraMultiVector > dgdxdot_mp = 
	outArgs.get_DgDx_dot_mp(0).getMultiVector();
      
      Teuchos::Array< Teuchos::RCP<ParamVec> > p_vec(outArgs.Np());
      Teuchos::Array< Teuchos::RCP< Stokhos::ProductEpetraMultiVector > > dgdp_mp(outArgs.Np());
      bool have_dgdp = false;
      for (int i=0; i<outArgs.Np(); i++) {
	dgdp_mp[i] = outArgs.get_DgDp_mp(0,i).getMultiVector();
	if (dgdp_mp[i] != Teuchos::null) {
	  have_dgdp = true;
	  Teuchos::Array<int> p_indexes = 
	    outArgs.get_DgDp_mp(0,i).getDerivativeMultiVector().getParamIndexes();
	  if (p_indexes.size() == 0)
	    p_vec[i] = Teuchos::rcp(&sacado_param_vec[i],false);
	  else {
	    p_vec[i] = Teuchos::rcp(new ParamVec);
	    for (int j=0; j<p_indexes.size(); j++)
	      p_vec[i]->addParam(sacado_param_vec[i][p_indexes[j]].family, 
				 sacado_param_vec[i][p_indexes[j]].baseValue);
	  }
	}
      }

      if (have_dgdp ||dgdx_mp != Teuchos::null || 
	  dgdxdot_mp != Teuchos::null) {
	app->evaluateMPResponseGradients(
	  x_dot_mp.get(), *x_mp, sacado_param_vec, p_mp_index, p_mp_vals, 
	  p_vec, g_mp.get(), dgdx_mp.get(), dgdxdot_mp.get(), dgdp_mp);
      }
      else if (g_mp != Teuchos::null)
	app->evaluateMPResponses(
	  x_dot_mp.get(), *x_mp, sacado_param_vec, p_mp_index, p_mp_vals, 
	  *g_mp);
    }
  }
}
