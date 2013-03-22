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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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

// NOX
#include "NOX.H"
#include "NOX_Epetra.H"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  int nelem = 100;
  int num_KL = 3;
  int p = 5;
  bool full_expansion = false;
  bool matrix_free = true;
  bool write_linear_system = false;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  int MyPID;

  try {

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

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
      Teuchos::rcp(new Teuchos::ParameterList);

    // Problem
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    problemParams.set("Name", "Heat Nonlinear Source");

    // Boundary conditions
    problemParams.set("Left BC", 0.0);
    problemParams.set("Right BC", 0.0);

    // Source function
    Teuchos::ParameterList& sourceParams = 
      problemParams.sublist("Source Function");
    sourceParams.set("Name", "Constant");
    sourceParams.set("Constant Value", 1.0);

    // Material
    Teuchos::ParameterList& matParams = 
      problemParams.sublist("Material Function");
    matParams.set("Name", "KL Exponential Random Field");
    matParams.set("Mean", 1.0);
    matParams.set("Standard Deviation", 0.5);
    matParams.set("Number of KL Terms", num_KL);
    Teuchos::Array<double> a(1), b(1), L(1);
    a[0] = 0.0; b[0] = 1.0; L[0] = 1.0;
    matParams.set("Domain Lower Bounds", a);
    matParams.set("Domain Upper Bounds", b);
    matParams.set("Correlation Lengths", L);

    // Response functions
    Teuchos::ParameterList& responseParams =
      problemParams.sublist("Response Functions");
    responseParams.set("Number", 1);
    responseParams.set("Response 0", "Solution Average");

    // Stochastic parameters
    Teuchos::ParameterList& parameterParams = 
      problemParams.sublist("Parameters");
    parameterParams.set("Number of Parameter Vectors", 1);
    Teuchos::ParameterList& pParams = 
      parameterParams.sublist("Parameter Vector 0");
    pParams.set("Number", num_KL);
    for (int i=0; i<num_KL; i++) {
      std::stringstream ss1, ss2;
      ss1 << "Parameter " << i;
      ss2 << "KL Exponential Function Random Variable " << i;
      pParams.set(ss1.str(), ss2.str());
    }

    // Mesh
    Teuchos::ParameterList& discParams = appParams->sublist("Discretization");
    discParams.set("Number of Elements", nelem);
    
    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (full_expansion)
      Cijk = basis->computeTripleProductTensor();
    else
      Cijk = basis->computeLinearTripleProductTensor();
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis,
									 Cijk));
    
    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    int num_spatial_procs = -1;
    if (argc > 1)
      num_spatial_procs = std::atoi(argv[1]);
    Teuchos::ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data =
      Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, globalComm,
					     parallelParams));
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = 
      sg_parallel_data->getMultiComm();
    Teuchos::RCP<const Epetra_Comm> app_comm = 
      sg_parallel_data->getSpatialComm();

    // Create application
    appParams->set("SG Method", "AD");
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(app_comm, appParams));
    
    // Create application model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, appParams));
    
    // Setup stochastic Galerkin algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(&(appParams->sublist("SG Parameters")),false);
    if (!full_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    Teuchos::ParameterList& sgOpParams = 
      sgParams->sublist("SG Operator");
    if (matrix_free)
      sgOpParams.set("Operator Method", "Matrix Free");
    else
      sgOpParams.set("Operator Method", "Fully Assembled");
    Teuchos::ParameterList& sgPrecParams = 
      sgParams->sublist("SG Preconditioner");
    sgPrecParams.set("Preconditioner Method", "Mean-based");
    sgPrecParams.set("Mean Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      sgPrecParams.sublist("Mean Preconditioner Parameters");
    precParams.set("default values", "SA");

    // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
						 expansion, sg_parallel_data, 
						 sgParams));

    // Set up stochastic parameters
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p_init =
      sg_model->create_p_sg(0);
    for (int i=0; i<num_KL; i++) {
      sg_p_init->term(i,0)[i] = 0.0;
      sg_p_init->term(i,1)[i] = 1.0;
    }
    sg_model->set_p_sg_init(0, *sg_p_init);

    // Setup stochastic initial guess
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init = 
      sg_model->create_x_sg();
    sg_x_init->init(0.0);
    sg_model->set_x_sg_init(*sg_x_init);

    // Set up NOX parameters
    Teuchos::RCP<Teuchos::ParameterList> noxParams =
      Teuchos::rcp(&(appParams->sublist("NOX")),false);

    // Set the nonlinear solver method
    noxParams->set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", 
                    NOX::Utils::OuterIteration + 
                    NOX::Utils::OuterIterationStatusTest + 
                    NOX::Utils::InnerIteration +
                    //NOX::Utils::Parameters + 
                    //NOX::Utils::Details + 
                    NOX::Utils::LinearSolverDetails +
                    NOX::Utils::Warning + 
                    NOX::Utils::Error);

    // Create printing utilities
    NOX::Utils utils(printParams);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = noxParams->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = noxParams->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 100);
    lsParams.set("Size of Krylov Subspace", 100);
    lsParams.set("Tolerance", 1e-12); 
    lsParams.set("Output Frequency", 10);
    lsParams.set("Preconditioner", "None");
    lsParams.set("Write Linear System", write_linear_system);

    // Sublist for convergence tests
    Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
    statusParams.set("Test Type", "Combo");
    statusParams.set("Number of Tests", 2);
    statusParams.set("Combo Type", "OR");
    Teuchos::ParameterList& normF = statusParams.sublist("Test 0");
    normF.set("Test Type", "NormF");
    normF.set("Tolerance", 1e-10);
    normF.set("Scale Type", "Scaled");
    Teuchos::ParameterList& maxIters = statusParams.sublist("Test 1");
    maxIters.set("Test Type", "MaxIters");
    maxIters.set("Maximum Iterations", 1);

    // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface = 
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));

    // Create NOX linear system object
    Teuchos::RCP<const Epetra_Vector> u = sg_model->get_x_init();
    Teuchos::RCP<Epetra_Operator> A = sg_model->create_W();
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
    if (matrix_free) {
      Teuchos::RCP<Epetra_Operator> M = sg_model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = nox_interface;
      lsParams.set("Preconditioner", "User Defined");
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							  iJac, A, iPrec, M,
							  *u));
    }
    else {
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							  iReq, iJac, A, 
							  *u));
    }

    // Build NOX group
    Teuchos::RCP<NOX::Epetra::Group> grp = 
      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));

    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
      NOX::StatusTest::buildStatusTests(statusParams, utils);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(grp, statusTests, noxParams);

    // Solve the system
    NOX::StatusTest::StatusType status = solver->solve();

    // Get final solution
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
      
    // Evaluate SG responses at SG parameters
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
      sg_model->createOutArgs();
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(0))));
    sg_inArgs.set_p(1, sg_p);
    sg_inArgs.set_x(Teuchos::rcp(&finalSolution,false));
    sg_outArgs.set_g(0, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    // Print mean and standard deviation
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_g_poly =
      sg_model->create_g_sg(0, View, sg_g.get());
    Epetra_Vector g_mean(*(model->get_g_map(0)));
    Epetra_Vector g_std_dev(*(model->get_g_map(0)));
    sg_g_poly->computeMean(g_mean);
    sg_g_poly->computeStandardDeviation(g_std_dev);
    std::cout << "\nResponse Expansion = " << std::endl;
    std::cout.precision(12);
    sg_g_poly->print(std::cout);
    std::cout << "\nResponse Mean =      " << std::endl << g_mean << std::endl;
    std::cout << "Response Std. Dev. = " << std::endl << g_std_dev << std::endl;
      
    if (status == NOX::StatusTest::Converged && MyPID == 0) 
      utils.out() << "Test Passed!" << std::endl;

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
