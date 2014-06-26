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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_Assert.hpp"
#include "FEApp_ProblemFactory.hpp"
#include "FEApp_BrusselatorProblem.hpp"
#include "FEApp_HeatNonlinearSourceProblem.hpp"
#include "FEApp_LinearConvDiffProblem.hpp"

FEApp::ProblemFactory::ProblemFactory(
       const Teuchos::RCP<Teuchos::ParameterList>& problemParams_,
       const Teuchos::RCP<ParamLib>& paramLib_) :
  problemParams(problemParams_),
  paramLib(paramLib_)
{
}

Teuchos::RCP<FEApp::AbstractProblem>
FEApp::ProblemFactory::create()
{
  Teuchos::RCP<FEApp::AbstractProblem> strategy;

  std::string& method = problemParams->get("Name", "Brusselator");
  if (method == "Brusselator") {
    strategy = Teuchos::rcp(new FEApp::BrusselatorProblem(problemParams, 
                                                          paramLib));
  }
  else if (method == "Heat Nonlinear Source") {
    strategy = 
      Teuchos::rcp(new FEApp::HeatNonlinearSourceProblem(problemParams,
                                                         paramLib));
  }
  else if (method == "Convection Diffusion") {
    strategy = 
      Teuchos::rcp(new FEApp::LinearConvDiffProblem(problemParams,
                                                    paramLib));
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                       std::endl << 
                       "Error!  Unknown problem " << method << 
                       "!" << std::endl << "Supplied parameter list is " << 
                       std::endl << *problemParams);
  }

  return strategy;
}
