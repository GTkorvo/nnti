//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include <iostream>
#include <sstream>

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
#include "FEApp_ModelEvaluator.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// TriKota
#include "TriKota_Driver.hpp"
#include "TriKota_DirectApplicInterface.hpp"
#include "TriKota_MPDirectApplicInterface.hpp"

// Piro solvers
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_StokhosMPSolver.hpp"

// Utilities
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;

  int MyPID;
  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs the FEApp Dakota solver.\n");

    std::string xml_filename = "input.xml";
    CLP.setOption("input", &xml_filename, "Input XML filename");
    bool use_multi_point = false;
    CLP.setOption("mp", "no-mp", &use_multi_point, "Use Multi-point method");
    CLP.parse( argc, argv );

    // Parse parameters
    RCP<ParameterList> appParams = 
      Teuchos::getParametersFromXmlFile(xml_filename);
    RCP<ParameterList> piroParams = 
      rcp(&(appParams->sublist("Piro")),false);
    ParameterList& dakotaParams = piroParams->sublist("Dakota");
    std::string dakota_input_file = 
      dakotaParams.get("Input File", "dakota.in");
    std::string dakota_output_file = 
      dakotaParams.get("Output File", "dakota.out");
    std::string dakota_restart_file = 
      dakotaParams.get("Restart File", "dakota.res");
    std::string dakota_error_file = 
      dakotaParams.get("Error File", "dakota.err");
    int p_index = dakotaParams.get("Parameter Vector Index", 0);
    int g_index = dakotaParams.get("Response Vector Index", 0);

    // Construct driver
    TriKota::Driver dakota(dakota_input_file.c_str(),
			   dakota_output_file.c_str(),
			   dakota_restart_file.c_str(),
			   dakota_error_file.c_str());

    // Create a communicator for Epetra objects
    RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = rcp(new Epetra_MpiComm(dakota.getAnalysisComm()));
#else
    globalComm = rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Construct a concrete Dakota interface with an EpetraExt::ModelEvaluator
    RCP<Dakota::DirectApplicInterface> trikota_interface;
    if (use_multi_point) {

      // Create MP solver
      RCP<ParameterList> mpParams = 
	rcp(&(dakotaParams.sublist("Multi-Point")),false);
      RCP<Piro::Epetra::StokhosMPSolver> mp_solver =
	rcp(new Piro::Epetra::StokhosMPSolver(
	      piroParams, mpParams, globalComm,
	      mpParams->get("Block Size", 10),
	      mpParams->get("Number of Spatial Processors", -1)));

      // Create application
      RCP<FEApp::Application> app = 
	rcp(new FEApp::Application(mp_solver->getSpatialComm(), appParams));
      
      // Create application model evaluator
      RCP<EpetraExt::ModelEvaluator> model = 
	rcp(new FEApp::ModelEvaluator(app, appParams)); 

      // Setup rest of solver
      mp_solver->setup(model);
      
      trikota_interface = 
	rcp(new TriKota::MPDirectApplicInterface(dakota.getProblemDescDB(), 
						 mp_solver, p_index, g_index), 
	    false);
    }
    else {
      // Create application
      RCP<FEApp::Application> app = 
	rcp(new FEApp::Application(globalComm, appParams));
      
      // Create application model evaluator
      RCP<EpetraExt::ModelEvaluator> model = 
	rcp(new FEApp::ModelEvaluator(app, appParams)); 

      // Create solver to map p -> g
      RCP<Piro::Epetra::NOXSolver> solver =
	rcp(new Piro::Epetra::NOXSolver(piroParams, model));

      trikota_interface = 
	rcp(new TriKota::DirectApplicInterface(dakota.getProblemDescDB(), 
					       solver), false);
    }

    // Run the requested Dakota strategy using this interface
    {
      TEUCHOS_FUNC_TIME_MONITOR("Total Dakota Analysis Time");
      dakota.run(trikota_interface.get());
    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

}
		  
		  
