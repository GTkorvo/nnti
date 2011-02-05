/* @HEADER@ */
/* ***********************************************************************
// 
//           Playa: Trilinos Solver Framework Extended
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#include "PlayaExceptions.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "PlayaAmesosSolver.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaBelosSolver.hpp"
#include "PlayaBICGSTABSolverDecl.hpp"
#include "PlayaBlockTriangularSolverDecl.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaBICGSTABSolverImpl.hpp"
#include "PlayaBlockTriangularSolverImpl.hpp"
#endif

using namespace Playa;
using namespace PlayaOps;
using namespace Teuchos;


LinearSolver<double> LinearSolverBuilder::createSolver(const std::string& filename)
{
  ParameterXMLFileReader reader(filename);
  ParameterList solverParams = reader.getParameters();
  return createSolver(solverParams);
}


LinearSolver<double> LinearSolverBuilder::createSolver(const ParameterList& params, int verb)
{
  TEST_FOR_EXCEPTION(!params.isSublist("Linear Solver"), std::runtime_error,
                     "did not find Linear Solver sublist in " << params);

  ParameterList solverSublist = params.sublist("Linear Solver");

  const std::string& solverType = getParameter<string>(solverSublist, "Type");

  Tabs tab;
  SUNDANCE_MSG1(verb, tab << "Solver builder creating a solver of type="
    << solverType);
  Tabs tab2;
  SUNDANCE_MSG2(verb, tab2 << "params = " << solverSublist);

  if (solverType=="Aztec")
    {
      return new AztecSolver(solverSublist);
    }
  else if (solverType=="Playa")
    {
      const std::string& solverMethod = getParameter<string>(solverSublist, "Method");
      if (solverMethod=="BICGSTAB") 
        {
          return new BICGSTABSolver<double>(solverSublist);
        }
      else if (solverMethod=="GMRES")
        {
          TEST_FOR_EXCEPTION(true, std::runtime_error, "Playa GMRES solver not implemented");
        }
    }
  else if (solverType=="Amesos")
    {
      return new AmesosSolver(solverSublist);
    }
  else if (solverType=="Belos")
    {
      return new BelosSolver(solverSublist);
    }
  else if (solverType=="Block Triangular")
    {
      ParameterList subSolverParams = solverSublist.sublist("Sub Solver");
      LinearSolver<double> subSolver = createSolver(subSolverParams);
      return new BlockTriangularSolver<double>(subSolver);
    }

  TEST_FOR_EXCEPTION(true, std::runtime_error, 
                     "Could not create a solver from parameter list " 
                     << params);
  return LinearSolver<double>();
    
}

