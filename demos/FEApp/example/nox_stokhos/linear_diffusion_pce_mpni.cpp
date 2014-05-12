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
//@HEADER

#include <iostream>
#include <sstream>

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
#include "FEApp_ModelEvaluator.hpp"
#include "Piro_Epetra_NOXSolver.hpp"
#include "BelosTypes.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// NOX
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "NOX_Epetra_LinearSystem_MPBD.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  int nelem = 100;
  int num_KL = 3;
  int p = 5;
  bool use_solver = false;
  std::string solver_type = "GMRES";

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

    // Alternative linear solver list for Stratimikos
    Teuchos::ParameterList& stratLinSolParams =
      newtonParams.sublist("Stratimikos Linear Solver");
    Teuchos::ParameterList& stratParams = 
      stratLinSolParams.sublist("Stratimikos");

    // Sublist for linear solver for the Newton method
    stratParams.set("Linear Solver Type", "Belos");
    Teuchos::ParameterList& belosParams = 
      stratParams.sublist("Linear Solver Types").sublist("Belos");
    Teuchos::ParameterList* belosSolverParams = NULL;
    if (solver_type == "GMRES") {
      belosParams.set("Solver Type","Block GMRES");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("Block GMRES"));
    }
    else if (solver_type == "CG") {
      belosParams.set("Solver Type","Block CG");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("Block CG"));
    }
    else if (solver_type == "RGMRES") {
      belosParams.set("Solver Type","GCRODR");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("GCRODR"));
      belosSolverParams->set("Num Recycled Blocks", 20);
    }
    else if (solver_type == "RCG") {
      belosParams.set("Solver Type","RCG");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("RCG"));
      Teuchos::RCP<const Teuchos::ParameterList> ortho_params = 
	  Teuchos::rcp(new Teuchos::ParameterList);
      belosSolverParams->set("Num Recycled Blocks", 10);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			 "Unknown solver type " << solver_type);
    belosSolverParams->set("Convergence Tolerance", 1e-12);
    belosSolverParams->set("Maximum Iterations", 1000);
    if (solver_type != "CG")
      belosSolverParams->set("Num Blocks", 100);
    belosSolverParams->set("Block Size", 1);
    belosSolverParams->set("Output Frequency",1);
    belosSolverParams->set("Output Style",1);
    //belosSolverParams->set("Verbosity",33);
    belosSolverParams->set("Verbosity", 
			   Belos::Errors + 
			   Belos::Warnings +
			   Belos::StatusTestDetails);
    //stratLinSolParams.set("Preconditioner", "User Defined");
    Teuchos::ParameterList& verboseParams = 
      belosParams.sublist("VerboseObject");
    verboseParams.set("Verbosity Level", "medium");

    stratParams.set("Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      stratParams.sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");
    precParams.set("default values", "SA");
    precParams.set("ML output", 0);
    precParams.set("max levels",5);
    precParams.set("increasing or decreasing","increasing");
    precParams.set("aggregation: type", "Uncoupled");
    precParams.set("smoother: type","ML symmetric Gauss-Seidel");
    precParams.set("smoother: sweeps",2);
    precParams.set("smoother: pre or post", "both");
    precParams.set("coarse: max size", 200);
#ifdef HAVE_ML_AMESOS
    precParams.set("coarse: type","Amesos-KLU");
#else
    precParams.set("coarse: type","Jacobi");
#endif

    //stratLinSolParams.set("Write Linear System", true);

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
    
    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    int sz = basis->size();
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    int num_mp = quad->size();
    
    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    int num_spatial_procs = -1;
    if (argc > 1)
      num_spatial_procs = std::atoi(argv[1]);
    Teuchos::RCP<const EpetraExt::MultiComm> multi_comm = 
      Stokhos::buildMultiComm(*globalComm, num_mp, num_spatial_procs);
    Teuchos::RCP<const Epetra_Comm> mp_comm = 
      Stokhos::getStochasticComm(multi_comm);
    Teuchos::RCP<const Epetra_Comm> app_comm = 
      Stokhos::getSpatialComm(multi_comm);
    Teuchos::RCP<const Epetra_Map> mp_block_map = 
      Teuchos::rcp(new Epetra_Map(num_mp, 0, *mp_comm));

    // Create application
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(app_comm, appParams));
    
    // Create application model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, appParams));
    
    // Wrap application model evaluator with an MP adapter
    Teuchos::RCP<EpetraExt::ModelEvaluator> mp_model = model;
      //Teuchos::rcp(new Stokhos::MPModelEvaluatorAdapter(model, num_mp));

    // Turn mp_model into an MP-nonlinear problem
    Teuchos::RCP<Teuchos::ParameterList> mpParams =
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& mpPrecParams = 
      mpParams->sublist("MP Preconditioner");
    //mpPrecParams.set("Preconditioner Method", "Block Diagonal");
    mpPrecParams.set("Preconditioner Method", "Mean-based");
    mpPrecParams.set("MP Preconditioner Type", "ML");
    Teuchos::ParameterList& pointPrecParams = 
      mpPrecParams.sublist("MP Preconditioner Parameters");
    pointPrecParams = precParams;
    Teuchos::RCP<Stokhos::MPModelEvaluator> mp_nonlinear_model =
      Teuchos::rcp(new Stokhos::MPModelEvaluator(mp_model, multi_comm,
						 mp_block_map, mpParams));

    // Create custom MPBD linear solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface;
    if (use_solver) {
      nox_interface = 
	Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(mp_nonlinear_model));
      Teuchos::RCP<Epetra_Operator> A = 
	mp_nonlinear_model->create_W();
      Teuchos::RCP<Epetra_Operator> M = 
	mp_nonlinear_model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = 
	nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = 
	nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = 
	nox_interface;
      Teuchos::ParameterList& outerSolParams = 
	newtonParams.sublist("Linear Solver");
      outerSolParams.sublist("Deterministic Solver Parameters") = 
	stratLinSolParams;
      outerSolParams.set("Preconditioner Strategy", "Mean");
      //outerSolParams.set("Preconditioner Strategy", "On the fly");
      Teuchos::RCP<Epetra_Operator> inner_A = model->create_W();
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> inner_nox_interface = 
	Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));
      Teuchos::RCP<NOX::Epetra::Interface::Required> inner_iReq = 
	inner_nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> inner_iJac = 
	inner_nox_interface;
      Teuchos::RCP<const Epetra_Vector> inner_u = model->get_x_init();
      Teuchos::RCP<NOX::Epetra::LinearSystem> inner_linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		       printParams, 
		       stratLinSolParams,
		       inner_iJac, inner_A, *inner_u));
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemMPBD(printParams, 
						       outerSolParams,
						       inner_linsys,
						       iReq, iJac, A,
						       model->get_x_map()));
    }

    // Create solver to map p -> g
    Teuchos::RCP<Piro::Epetra::NOXSolver> mp_solver =
      Teuchos::rcp(new Piro::Epetra::NOXSolver(appParams, mp_nonlinear_model,
					       Teuchos::null, nox_interface,
					       linsys));

    // Create MP inverse model evaluator to map p_mp -> g_mp
    Teuchos::Array<int> mp_p_index_map = 
      mp_nonlinear_model->get_p_mp_map_indices();
    Teuchos::Array<int> mp_g_index_map = 
      mp_nonlinear_model->get_g_mp_map_indices();
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_g_maps = 
      mp_nonlinear_model->get_g_mp_base_maps();
    mp_g_index_map.push_back(base_g_maps.size());
    base_g_maps.push_back(model->get_x_map());
    Teuchos::RCP<EpetraExt::ModelEvaluator> mp_inverse_solver =
      Teuchos::rcp(new Stokhos::MPInverseModelEvaluator(mp_solver,
							mp_p_index_map,
							mp_g_index_map,
							base_g_maps));

    // Create MP-based SG Quadrature model evaluator to calculate g_sg
    Teuchos::RCP<EpetraExt::ModelEvaluator> sg_solver =
      Teuchos::rcp(new Stokhos::SGQuadMPModelEvaluator(mp_inverse_solver, 
						       multi_comm, 
						       mp_block_map));

    // Set up stochastic parameters
    Teuchos::ParameterList parallelParams;
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data =
      Teuchos::rcp(new Stokhos::ParallelData(basis, Teuchos::null, globalComm,
					     parallelParams));
    Teuchos::RCP<Epetra_LocalMap> stoch_map =
      Teuchos::rcp(new Epetra_LocalMap(sz, 0, *globalComm));
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p_init = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, stoch_map, sg_solver->get_p_map(0), 
		     sg_parallel_data->getMultiComm()));
    for (int i=0; i<num_KL; i++) {
      sg_p_init->term(i,0)[i] = 0.0;
      sg_p_init->term(i,1)[i] = 1.0;
    }

    // Evaluate SG responses at SG parameters
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = 
      sg_solver->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
      sg_solver->createOutArgs();
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_g = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, stoch_map, sg_solver->get_g_map(0), 
		     sg_parallel_data->getMultiComm()));
    sg_inArgs.set_sg_basis(basis);
    sg_inArgs.set_sg_quadrature(quad);
    sg_inArgs.set_p_sg(0, sg_p_init);
    sg_outArgs.set_g_sg(0, sg_g);
    sg_solver->evalModel(sg_inArgs, sg_outArgs);

    // Print mean and standard deviation
    std::cout.precision(12);
    std::cout << "SG expansion of response:" << std::endl << *sg_g;
    Epetra_Vector mean(*(sg_solver->get_g_map(0)));
    Epetra_Vector std_dev(*(sg_solver->get_g_map(0)));
    sg_g->computeMean(mean);
    sg_g->computeStandardDeviation(std_dev);
    std::cout << "Mean =      " << mean[0] << std::endl;
    std::cout << "Std. Dev. = " << std_dev[0] << std::endl;
   
    if (!sg_outArgs.isFailed() && MyPID == 0) 
      std::cout << "Test Passed!" << std::endl;

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
