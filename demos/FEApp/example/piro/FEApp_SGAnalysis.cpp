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

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Piro solvers
#include "Piro_Epetra_StokhosSolverFactory.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_PerformAnalysis.hpp"
#include "Thyra_VectorBase.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

int testResponses(const Epetra_Vector& g, 
		  const Teuchos::Array<double> testValues,
		  double absTol, double relTol,
		  const std::string& tag,
		  bool print_proc)
{
  int failures = 0;
  TEUCHOS_TEST_FOR_EXCEPTION(g.MyLength() != testValues.size(),
		     std::logic_error,
		     tag << " Test Values array has size " << 
		     testValues.size() << "but expected size " <<
		     g.MyLength());
  for (int i=0; i<testValues.size(); i++) {
    if (std::abs(g[i]-testValues[i]) > 
	relTol*std::abs(testValues[i])+absTol) {
      ++failures;
      if (print_proc)
	std::cout << tag << " test " << i << " failed:  Expected:  " 
		  << testValues[i] << ", got:  " << g[i]
		  << ", abs tol = " << absTol << ", rel tol = " << relTol
		  << "." << std::endl;
    }
  }

  return failures;
}

int main(int argc, char *argv[]) {
  int MyPID;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs the FEApp stochastic Galerkin solver.\n");

    std::string xml_filename = "input.xml";
    CLP.setOption("input", &xml_filename, "Input XML filename");

    CLP.parse( argc, argv );

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total FEAppSG Run Time");

    // Create a communicator for Epetra objects
    Teuchos::RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = Teuchos::rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Set up application parameters
    Teuchos::RCP<Teuchos::ParameterList> appParams = 
      Teuchos::getParametersFromXmlFile(xml_filename);
    
    // Create stochastic Galerkin solver factory
    Teuchos::RCP< Teuchos::ParameterList> piroParams = 
      Teuchos::rcp(&(appParams->sublist("Piro")),false);
    Piro::Epetra::StokhosSolverFactory sg_solver_factory(piroParams, 
							 globalComm);

    // Get comm for spatial problem
    Teuchos::RCP<const Epetra_Comm> app_comm = 
      sg_solver_factory.getSpatialComm();

    // Create application
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(app_comm, appParams));
    
    // Create application model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, appParams)); 

    // Setup rest of solver
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model = 
      sg_solver_factory.createSGModel(model);
    Teuchos::RCP<EpetraExt::ModelEvaluator> sg_solver =
      sg_solver_factory.createSGSolver(sg_model);
    Teuchos::RCP<EpetraExt::ModelEvaluator> rs_model =
      sg_solver_factory.createRSModel(sg_solver);

    Thyra::EpetraModelEvaluator rs_model_thyra;
    rs_model_thyra.initialize(rs_model, Teuchos::null);

    Teuchos::RCP< ::Thyra::VectorBase<double> > p;
    Teuchos::ParameterList& analysisParams = piroParams->sublist("Analysis");
    int status = 
      Piro::PerformAnalysis(rs_model_thyra, analysisParams, p); 

    std::cout << *p;

    int failures = 0;
    /*
    // Regression tests
    int failures = 0;
    Teuchos::ParameterList& testParams = appParams->sublist("Regression Tests");
    double relTol = testParams.get("Relative Tolerance", 1.0e-3);
    double absTol = testParams.get("Absolute Tolerance", 1.0e-8);

    for (int i=0; i<ng-1; i++) {
      // Don't loop over last g which is x, since it is a long vector
      // to print out.
      if (sg_outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_g_sg, i)) {

	// Print mean and standard deviation      
	Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> g_sg = 
	  sg_outArgs.get_g_sg(i);
	Epetra_Vector g_mean(*(sg_solver->get_g_map(i)));
	Epetra_Vector g_std_dev(*(sg_solver->get_g_map(i)));
	g_sg->computeMean(g_mean);
	g_sg->computeStandardDeviation(g_std_dev);
	std::cout.precision(12);
	std::cout << "Response " << i << " Mean =      " << std::endl 
		  << g_mean << std::endl;
	std::cout << "Response " << i << " Std. Dev. = " << std::endl 
		  << g_std_dev << std::endl;

	// Test mean
	std::stringstream ss1;
	ss1 << "Response " << i << " Mean Test Values";
	bool testMean = 
	  testParams.isType< Teuchos::Array<double> >(ss1.str());
	if (testMean) { 
	  Teuchos::Array<double> testValues =
	    testParams.get<Teuchos::Array<double> >(ss1.str());
	  failures += testResponses(g_mean, testValues, absTol, relTol, "Mean", 
				    MyPID==0);
	}

	// Test std. dev.
	std::stringstream ss2;
	ss2 << "Response " << i << " Standard Deviation Test Values";
	bool testSD = 
	  testParams.isType< Teuchos::Array<double> >(ss2.str());
	if (testSD) { 
	  Teuchos::Array<double> testValues =
	    testParams.get<Teuchos::Array<double> >(ss2.str());
	  failures += testResponses(g_std_dev, testValues, absTol, relTol, 
				    "Standard Deviation", MyPID==0);
	}

      }
    }
    */

    if (MyPID == 0) {
      if (failures == 0)
	std::cout << "Test Passed!" << std::endl;
      else
	std::cout << "Test Failed!" << std::endl;
    }

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (std::string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
		  
		  
