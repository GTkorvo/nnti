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

// AztecOO solver
#include "AztecOO.h"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  int nelem = 100;
  int num_KL = 3;
  int p = 5;
  bool full_expansion = false;

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
    sgOpParams.set("Operator Method", "Matrix Free");
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

    // Create vectors and operators
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> sg_x = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    *sg_x = *(sg_model->get_x_init());
    Teuchos::RCP<Epetra_Vector> sg_f = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_f_map())));
    Teuchos::RCP<Epetra_Vector> sg_dx = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    Teuchos::RCP<Epetra_Operator> sg_J = sg_model->create_W();
    Teuchos::RCP<Epetra_Operator> sg_M = sg_model->create_WPrec()->PrecOp;

    // Setup InArgs and OutArgs
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = sg_model->createOutArgs();
    sg_inArgs.set_p(1, sg_p);
    sg_inArgs.set_x(sg_x);
    sg_outArgs.set_f(sg_f);
    sg_outArgs.set_W(sg_J);
    sg_outArgs.set_WPrec(sg_M);

    // Evaluate model
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    // Print initial residual norm
    double norm_f;
    sg_f->Norm2(&norm_f);
    if (MyPID == 0)
      std::cout << "\nInitial residual norm = " << norm_f << std::endl;

    // Setup AztecOO solver
    AztecOO aztec;
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
    aztec.SetAztecOption(AZ_precond, AZ_none);
    aztec.SetAztecOption(AZ_kspace, 20);
    aztec.SetAztecOption(AZ_conv, AZ_r0);
    aztec.SetAztecOption(AZ_output, 10);
    aztec.SetUserOperator(sg_J.get());
    aztec.SetPrecOperator(sg_M.get());
    aztec.SetLHS(sg_dx.get());
    aztec.SetRHS(sg_f.get());

    // Solve linear system
    aztec.Iterate(100, 1e-12);

    // Update x
    sg_x->Update(-1.0, *sg_dx, 1.0);

    // Compute new residual & response function
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(0))));
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs2 = sg_model->createOutArgs();
    sg_outArgs2.set_f(sg_f);
    sg_outArgs2.set_g(0, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs2);

    // Print initial residual norm
    sg_f->Norm2(&norm_f);
    if (MyPID == 0)
      std::cout << "\nFinal residual norm = " << norm_f << std::endl;

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

    if (norm_f < 1.0e-10 && MyPID == 0)
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
