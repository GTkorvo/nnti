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

// Piro solver
#include "Piro_Epetra_Factory.hpp"

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
      "This example runs the FEApp solver.\n");

    std::string xml_filename = "input.xml";
    CLP.setOption("input", &xml_filename, "Input XML filename");

    CLP.parse( argc, argv );

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total FEApp Run Time");

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

    // Create application
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(globalComm, appParams));

    // Create application model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, appParams)); 
    
    // Create Piro solver
    Teuchos::RCP< Teuchos::ParameterList> piroParams = 
      Teuchos::rcp(&(appParams->sublist("Piro")),false);
    Teuchos::RCP<EpetraExt::ModelEvaluator> solver =
      Piro::Epetra::Factory::createSolver(piroParams, model);

    // Evaluate SG responses at SG parameters
    EpetraExt::ModelEvaluator::InArgs inArgs = solver->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = solver->createOutArgs();
    for (int i=0; i<inArgs.Np(); i++) {
      Teuchos::RCP<const Epetra_Vector> p = solver->get_p_init(i);
      inArgs.set_p(i, p);
    }
    for (int i=0; i<outArgs.Ng(); i++) {
      Teuchos::RCP<Epetra_Vector> g = 
	Teuchos::rcp(new Epetra_Vector(*(solver->get_g_map(i))));
      outArgs.set_g(i, g);
    }
    solver->evalModel(inArgs, outArgs);

    // Print responses (not last one since that is x)
    std::cout << std::endl;
    std::cout.precision(8);
    for (int i=0; i<outArgs.Ng()-1; i++) {
      if (outArgs.get_g(i) != Teuchos::null)
	std::cout << "Response " << i << " = " << std::endl 
		  << *(outArgs.get_g(i)) << std::endl;
    }

    // Regression tests
    Teuchos::ParameterList& testParams = appParams->sublist("Regression Tests");
    int failures = 0;
    double relTol = testParams.get("Relative Tolerance", 1.0e-3);
    double absTol = testParams.get("Absolute Tolerance", 1.0e-8);
    
    // Test responses
    for (int i=0; i<outArgs.Ng(); i++) {
      std::stringstream ss1, ss2;
      ss1 << "Response Test Values " << i;
      ss2 << "Response " << i;
      bool testResp = 
	testParams.isType< Teuchos::Array<double> >(ss1.str());
      if (testResp) { 
	Teuchos::Array<double> testValues =
	  testParams.get<Teuchos::Array<double> >(ss1.str());
	failures += testResponses(*(outArgs.get_g(i)), testValues, 
				  absTol, relTol, ss2.str(), MyPID==0);
      }
    }

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
		  
		  
