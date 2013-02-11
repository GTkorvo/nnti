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
//@HEADER

#include <iostream>
#include <sstream>

#include "NOX.H"
#include "NOX_Epetra.H"

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
#include "FEApp_ModelEvaluator.hpp"
#include "Piro_Epetra_NOXSolver.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Stokhos_Epetra.hpp"
#include "Stokhos_IfpackPreconditionerFactory.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "AztecOO.h"

#include "Stokhos_PCEAnasaziKL.hpp"
#include "Stokhos_StieltjesGramSchmidtBuilder.hpp"
#include "EpetraExt_MultiVectorOut.h"
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "EpetraExt_BlockUtility.h"

int main(int argc, char *argv[]) {
  int nelem = 100;
  double alpha = 2.0;
  double leftBC = 0.0;
  double rightBC = 0.1;
  int numalpha = 3;
  unsigned int p = 7;
  unsigned int d = numalpha;

  int MyPID;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
    Teuchos::RCP<Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = Teuchos::rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Create SG basis to setup parallel correctly
    typedef Stokhos::LegendreBasis<int,double> basis_type;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (unsigned int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Create multi-level comm and spatial comm
    int num_spatial_procs = -1;
    if (argc > 1)
      num_spatial_procs = std::atoi(argv[1]);
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm =
      Stokhos::buildMultiComm(*globalComm, basis->size(), num_spatial_procs);
    Teuchos::RCP<const Epetra_Comm> app_comm = Stokhos::getSpatialComm(sg_comm);

    // Set up application parameters
    Teuchos::RCP<Teuchos::ParameterList> appParams = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Problem
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    problemParams.set("Name", "Heat Nonlinear Source");

    // Boundary conditions
    problemParams.set("Left BC", leftBC);
    problemParams.set("Right BC", rightBC);

    // Source function
    Teuchos::ParameterList& sourceParams = 
      problemParams.sublist("Source Function");
    sourceParams.set("Name", "Multi-Variate Exponential");
    sourceParams.set<unsigned int>("Nonlinear Factor Dimensions", numalpha);
    for (int i=0; i<numalpha; i++) {
      std::stringstream ss;
      ss << "Nonlinear Factor " << i;
      sourceParams.set(ss.str(), alpha/numalpha);
    }

    // Material
    Teuchos::ParameterList& matParams = 
      problemParams.sublist("Material Function");
    matParams.set("Name", "Constant");
    matParams.set("Constant Value", 1.0);

    // Response functions
    Teuchos::ParameterList& responseParams =
      problemParams.sublist("Response Functions");
    responseParams.set("Number", 1);
    responseParams.set("Response 0", "Solution Average");

    // Free parameters (determinisic, e.g., for sensitivities)
    Teuchos::ParameterList& parameterParams = 
      problemParams.sublist("Parameters");
    parameterParams.set("Number of Parameter Vectors", 1);
    Teuchos::ParameterList& pParams = 
      parameterParams.sublist("Parameter Vector 0");
    pParams.set("Number", 1);
    pParams.set("Parameter 0", "Constant Function Value");

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
    lsParams.set("Max Iterations", 20);
    lsParams.set("Size of Krylov Subspace", 20);
    lsParams.set("Tolerance", 1e-4); 
    lsParams.set("Output Frequency", 50);
    lsParams.set("Preconditioner", "Ifpack");
    lsParams.set("RCM Reordering", "Enabled");

    // Sublist for convergence tests
    Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
    statusParams.set("Test Type", "Combo");
    statusParams.set("Number of Tests", 2);
    statusParams.set("Combo Type", "OR");
    Teuchos::ParameterList& comboParams = statusParams.sublist("Test 0");
    comboParams.set("Test Type", "Combo");
    comboParams.set("Number of Tests", 2);
    comboParams.set("Combo Type", "AND");
    Teuchos::ParameterList& normF = comboParams.sublist("Test 0");
    normF.set("Test Type", "NormF");
    normF.set("Tolerance", 1e-10);
    normF.set("Scale Type", "Scaled");
    Teuchos::ParameterList& normWRMS = comboParams.sublist("Test 1");
    normWRMS.set("Test Type", "NormWRMS");
    normWRMS.set("Relative Tolerance", 1e-6);
    normWRMS.set("Absolute Tolerance", 1e-8);
    Teuchos::ParameterList& maxIters = statusParams.sublist("Test 1");
    maxIters.set("Test Type", "MaxIters");
    maxIters.set("Maximum Iterations", 10);

    // Create application
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(app_comm, appParams));

    // Create model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, appParams));

    // Create NOX solver
    Piro::Epetra::NOXSolver solver(appParams, model);

    // Evaluate responses at parameters
    EpetraExt::ModelEvaluator::InArgs inArgs = solver.createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = solver.createOutArgs();
    Teuchos::RCP<const Epetra_Vector> p_init = solver.get_p_init(0);
    Teuchos::RCP<Epetra_Vector> g = 
      Teuchos::rcp(new Epetra_Vector(*(solver.get_g_map(0))));
    Teuchos::RCP<Epetra_Vector> finalSolution = 
      Teuchos::rcp(new Epetra_Vector(*(solver.get_g_map(1))));
    Teuchos::RCP<Epetra_MultiVector> dgdp = 
      Teuchos::rcp(new Epetra_MultiVector(*(solver.get_g_map(0)),
					  p_init->MyLength()));
    inArgs.set_p(0, p_init);
    outArgs.set_g(0, g);
    outArgs.set_g(1, finalSolution);
    outArgs.set_DgDp(0, 0, dgdp);
    solver.evalModel(inArgs, outArgs);

    g->Print(utils.out());
    dgdp->Print(utils.out());

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

    //TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create SG quadrature and expansion
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis, p));
    unsigned int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, 
								    Cijk,
								    quad));
    // Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
    // 	Teuchos::rcp(new Stokhos::ForUQTKOrthogPolyExpansion<int,double>(basis, 
    // 								      Stokhos::ForUQTKOrthogPolyExpansion<int,double>::INTEGRATION, 1e-6));
    
    utils.out() << "sz = " << sz << std::endl;
    appParams->set("SG Method", "AD");

    // Create stochastic parallel distribution
    Teuchos::ParameterList parallelParams;
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data =
      Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, sg_comm,
					     parallelParams));
    Teuchos::RCP<const EpetraExt::MultiComm> multi_comm = 
      sg_parallel_data->getMultiComm();
    Teuchos::RCP<const Epetra_Comm> stoch_comm =
      sg_parallel_data->getStochasticComm();
    Teuchos::RCP<const Epetra_BlockMap> stoch_overlap_map = 
      Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(sz), 0, *stoch_comm));

    // Stochastic parameters
     parameterParams.set("Number of Parameter Vectors", 2);
     Teuchos::ParameterList& pParams2 = 
       parameterParams.sublist("Parameter Vector 1");
    pParams2.set("Number", numalpha);
    for (int i=0; i<numalpha; i++) {
      std::stringstream ss1, ss2;
      ss1 << "Parameter " << i;
      ss2 << "Exponential Source Function Nonlinear Factor " << i;
      pParams2.set(ss1.str(), ss2.str());
    }
      
    // Create new app for Stochastic Galerkin solve
    app = Teuchos::rcp(new FEApp::Application(app_comm, appParams,
					      finalSolution.get()));

    model = Teuchos::rcp(new FEApp::ModelEvaluator(app, appParams));

    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(&(appParams->sublist("SG Parameters")),false);
     Teuchos::ParameterList& sgOpParams = 
      sgParams->sublist("SG Operator");
    Teuchos::ParameterList& sgPrecParams = 
      sgParams->sublist("SG Preconditioner");
    sgOpParams.set("Operator Method", "Matrix Free");
    sgPrecParams.set("Preconditioner Method", "Mean-based");
     sgPrecParams.set("Mean Preconditioner Type", "Ifpack");
    Teuchos::ParameterList& precParams = 
      sgPrecParams.sublist("Mean Preconditioner Parameters");
    precParams.set("Ifpack Preconditioner", "ILU");
    sgParams->set("Evaluate W with F", false);
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, quad,
						 expansion, sg_parallel_data, 
						 sgParams));

    // Set up stochastic parameters
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p = 
      sg_model->create_p_sg(1);
    for (unsigned int i=0; i<d; i++) {
      sg_p->term(i,0)[i] = 2.0;
      sg_p->term(i,1)[i] = 1.0;
    }
    sg_model->set_p_sg_init(1, *sg_p);

    // Setup stochastic initial guess
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init =
      sg_model->create_x_sg();
    sg_x_init->init(0.0);
    if (sg_x_init->myGID(0))
      (*sg_x_init)[0] = *finalSolution;
    sg_model->set_x_sg_init(*sg_x_init);
  
    // Evaluate SG responses at SG parameters
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = sg_model->createOutArgs();
    Teuchos::RCP<const Epetra_Vector> sg_p_init = sg_model->get_p_init(2);
    Teuchos::RCP<Epetra_Vector> sg_x = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    *sg_x = *(sg_model->get_x_init());
    Teuchos::RCP<Epetra_Vector> sg_f = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_f_map())));
    Teuchos::RCP<Epetra_Vector> sg_dx = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    Teuchos::RCP<Epetra_Operator> sg_J = sg_model->create_W();
    Teuchos::RCP<Epetra_Operator> sg_M = sg_model->create_WPrec()->PrecOp;
    sg_inArgs.set_p(2, sg_p_init);
    sg_inArgs.set_x(sg_x);
    sg_outArgs.set_f(sg_f);
    sg_outArgs.set_W(sg_J);
    sg_outArgs.set_WPrec(sg_M);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    AztecOO aztec;
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
    aztec.SetAztecOption(AZ_precond, AZ_none);
    aztec.SetAztecOption(AZ_kspace, 20);
    aztec.SetAztecOption(AZ_conv, AZ_r0);
    aztec.SetAztecOption(AZ_output, 1);
    aztec.SetUserOperator(sg_J.get());
    aztec.SetPrecOperator(sg_M.get());
    aztec.SetLHS(sg_dx.get());
    aztec.SetRHS(sg_f.get());
    aztec.Iterate(20, 1e-4);

    sg_x->Update(-1.0, *sg_dx, 1.0);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    /*
    sg_dx->PutScalar(0.0);
    aztec.SetUserOperator(sg_J.get());
    aztec.SetPrecOperator(sg_M.get());
    aztec.SetLHS(sg_dx.get());
    aztec.SetRHS(sg_f.get());
    aztec.Iterate(20, 1e-4);
    */

    double normf;
    sg_f->Norm2(&normf);
    int nit = 0;
    utils.out() << "nit = " << nit << " normf = " << normf << std::endl;
    double norm_dx = 1.0;
    while (norm_dx > 1e-5 && nit < 5) {

    // Get Jacobian blocks
    Teuchos::RCP<Stokhos::SGOperator> sg_J_op = 
      Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(sg_J);
    Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_J_poly =
      sg_J_op->getSGPolynomial();
    Teuchos::RCP<const Epetra_CrsMatrix> J_mean = 
      Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(
	sg_J_poly->getCoeffPtr(0));

    // Build a vector polynomial out of matrix nonzeros
    Epetra_Map J_vec_map(-1, J_mean->NumMyNonzeros(), 0, 
			 sg_comm->SubDomainComm());
    Stokhos::VectorOrthogPoly<Epetra_Vector> sg_J_vec_poly(
      basis, stoch_overlap_map, Stokhos::EpetraVectorCloner(J_vec_map));
    for (unsigned int coeff=0; coeff<sz; coeff++) {
      Teuchos::RCP<const Epetra_CrsMatrix> J_coeff = 
	Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>
	(sg_J_poly->getCoeffPtr(coeff));
      int row = 0;
      for (int i=0; i<J_mean->NumMyRows(); i++) {
	int num_col;
	J_mean->NumMyRowEntries(i, num_col);
	for (int j=0; j<num_col; j++)
	  sg_J_vec_poly[coeff][row++] = (*J_coeff)[i][j];
      }
    }

    // Compute KL expansion of solution sg_J_vec_poly
    int num_KL = 3;
    Stokhos::PCEAnasaziKL pceKL(sg_J_vec_poly, num_KL);
    Teuchos::ParameterList anasazi_params = pceKL.getDefaultParams();
    bool result = pceKL.computeKL(anasazi_params);
    if (!result)
      utils.out() << "KL Eigensolver did not converge!" << std::endl;
    Teuchos::Array<double> evals = pceKL.getEigenvalues();
    Teuchos::RCP<Epetra_MultiVector> evecs = pceKL.getEigenvectors();

    //EpetraExt::MultiVectorToMatrixMarketFile("eigenvectors.mm", *evecs);

    // Compute PCE expansions of new random variables
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > rv_pce(num_KL);
    for (int rv=0; rv < num_KL; rv++) {
      rv_pce[rv].reset(basis);
      rv_pce[rv][0] = 0.0;
      for (unsigned int coeff=1; coeff < sz; coeff++) {
	double dot;
	sg_J_vec_poly[coeff].Dot(*((*evecs)(rv)), &dot);
	rv_pce[rv][coeff] = dot/std::sqrt(evals[rv]);
      }
      //utils.out() << "rv[" << rv << "] = " << rv_pce[rv] << std::endl;
    }

    // Compute Stieltjes-Gram-Schmidt basis for KL expansion
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad2 = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad2 = 
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,p));
    //int p2 = p;
    unsigned int p2 = 2;
    if (p2 > p)
      p2 = p;
    Stokhos::StieltjesGramSchmidtBuilder<int,double> gs_builder(quad2, rv_pce, 
								p2, true, 
								false);
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > kl_basis =
      gs_builder.getReducedBasis();
    Teuchos::RCP<Stokhos::Quadrature<int,double> > kl_quad =
      gs_builder.getReducedQuadrature();
    unsigned int sz2 = kl_basis->size();
    utils.out() << "sz2 = " << sz2 << std::endl;

    // Compute new block vectors
    Teuchos::RCP<const EpetraExt::MultiComm> kl_comm =
      Stokhos::buildMultiComm(*globalComm, sz2, num_spatial_procs);
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > kl_Cijk =
      kl_basis->computeLinearTripleProductTensor();
    Teuchos::RCP<Stokhos::EpetraSparse3Tensor> kl_epetraCijk =
      Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(kl_basis, kl_Cijk, 
						    kl_comm));
    kl_epetraCijk->transformToLocal();
    Teuchos::RCP<const Epetra_BlockMap> kl_stoch_row_map =
      kl_epetraCijk->getStochasticRowMap();
    Teuchos::RCP<const Epetra_BlockMap> kl_ov_stoch_row_map = 
      Teuchos::rcp(new Epetra_LocalMap(static_cast<int>(sz2), 0, kl_comm->TimeDomainComm()));
    
    Teuchos::RCP<const Epetra_Map> base_x_map = model->get_x_map();
    Teuchos::RCP<const Epetra_Map> base_f_map = model->get_f_map();
    Teuchos::RCP<Epetra_Map> sg_x_kl_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*base_x_map,
    							     *kl_stoch_row_map,
    							     *kl_comm));
    Teuchos::RCP<Epetra_Map> sg_ov_x_kl_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*base_x_map,
    							     *kl_ov_stoch_row_map,
    							     *kl_comm));
    Teuchos::RCP<Epetra_Map> sg_f_kl_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*base_f_map,
    							     *kl_stoch_row_map,
    							     *kl_comm));
    Teuchos::RCP<Epetra_Map> sg_ov_f_kl_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*base_f_map,
    							     *kl_ov_stoch_row_map,
    							     *kl_comm));
    // Teuchos::RCP<const Epetra_Map> sg_x_kl_map = sg_model->get_x_map();
    // Teuchos::RCP<const Epetra_Map> sg_f_kl_map = sg_model->get_f_map();

    Epetra_Import sg_x_kl_importer(*sg_ov_x_kl_map, *sg_x_kl_map);
    Epetra_Import sg_f_kl_importer(*sg_ov_f_kl_map, *sg_f_kl_map);

    // Map Jacobian to KL basis
    utils.out() << "Mapping Jacobian to KL basis..." << std::endl;
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > rv_kl_pce(num_KL);
    gs_builder.computeReducedPCEs(rv_pce, rv_kl_pce);
    // for (unsigned int rv=0; rv < num_KL; rv++)
    //   utils.out() << "rv_kl[" << rv << "] = " << rv_kl_pce[rv] << std::endl;
    Stokhos::VectorOrthogPoly<Epetra_Vector> sg_J_kl_vec_poly(
      kl_basis, kl_ov_stoch_row_map, Stokhos::EpetraVectorCloner(J_vec_map));
    for (int rv=0; rv < num_KL; rv++) {
      for (int i=0; i<sg_J_kl_vec_poly[0].MyLength(); i++) {
	sg_J_kl_vec_poly[0][i] = sg_J_vec_poly[0][i] + 
	  std::sqrt(evals[rv])*(*evecs)[rv][i]*rv_kl_pce[rv][0];
	for (int j=1; j<num_KL; j++)
	  sg_J_kl_vec_poly[j][i] +=  
	    std::sqrt(evals[rv])*(*evecs)[rv][i]*rv_kl_pce[rv][j];
      }
    }
    Stokhos::EpetraOperatorOrthogPoly sg_J_kl_poly(
      kl_basis, kl_ov_stoch_row_map, base_x_map, base_f_map, 
      sg_f_kl_map, kl_comm);
    for (int coeff = 0; coeff < kl_basis->size(); coeff++) {
      Teuchos::RCP<Epetra_CrsMatrix> mat = 
	Teuchos::rcp(new Epetra_CrsMatrix(*J_mean));
      int row = 0;
      for (int i=0; i<J_mean->NumMyRows(); i++) {
	int num_col;
	J_mean->NumMyRowEntries(i, num_col);
	for (int j=0; j<num_col; j++)
	  (*mat)[i][j] = sg_J_kl_vec_poly[coeff][row++];
      }
      mat->FillComplete();
      sg_J_kl_poly.setCoeffPtr(coeff, mat);
    }

    // Map RHS to KL basis
    utils.out() << "Mapping RHS to KL basis..." << std::endl;
    const Teuchos::Array<double>& weights = quad2->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
      quad2->getBasisAtQuadPoints();
    const Teuchos::Array< Teuchos::Array<double> >& kl_basis_vals = 
      kl_quad->getBasisAtQuadPoints();
    int nqp = weights.size();
    const Teuchos::Array<double>& kl_norms = kl_basis->norm_squared();
    Teuchos::RCP<EpetraExt::BlockVector> sg_f_ov = 
      sg_model->import_residual(*sg_f);
    Stokhos::EpetraVectorOrthogPoly sg_f_poly(
      basis, stoch_overlap_map, base_f_map, 
      Teuchos::rcp(&(sg_f_ov->Map()), false),
      multi_comm, View, *sg_f_ov);
    Stokhos::EpetraVectorOrthogPoly sg_f_kl_poly(
      kl_basis, kl_ov_stoch_row_map, base_f_map, sg_ov_f_kl_map, kl_comm);
    Epetra_Vector f_val(*base_f_map);
    for (int qp=0; qp < nqp; qp++) {
      sg_f_poly.evaluate(basis_vals[qp], f_val);
      sg_f_kl_poly.sumIntoAllTerms(weights[qp], kl_basis_vals[qp], kl_norms, 
				   f_val);
    }
    Teuchos::RCP<EpetraExt::BlockVector> sg_ov_f_kl_block = 
      sg_f_kl_poly.getBlockVector();
    EpetraExt::BlockVector sg_f_kl_block(*base_f_map, *sg_f_kl_map);
    sg_f_kl_block.Export(*sg_ov_f_kl_block, sg_f_kl_importer, Insert);

    // Create linear system operator
    utils.out() << "Setting up operators..." << std::endl;
    Teuchos::RCP<Teuchos::ParameterList> pl = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Stokhos::MatrixFreeOperator kl_jac_op(
      kl_comm, kl_basis, kl_epetraCijk, base_x_map, base_f_map, 
      sg_x_kl_map, sg_f_kl_map, pl);
    kl_jac_op.setupOperator(Teuchos::rcp(&sg_J_kl_poly,false));
    Teuchos::RCP<Stokhos::IfpackPreconditionerFactory> ifpack_factory = 
      Teuchos::rcp(new Stokhos::IfpackPreconditionerFactory(Teuchos::rcp(&precParams,false)));
    Stokhos::MeanBasedPreconditioner kl_prec_op(
      kl_comm, kl_basis, kl_epetraCijk, base_x_map, sg_x_kl_map, ifpack_factory,
      pl);
    kl_prec_op.setupPreconditioner(Teuchos::rcp(&kl_jac_op,false), 
				   *sg_x);
				   

    // Solve KL linear system
    utils.out() << "Solving linear system..." << std::endl;
    EpetraExt::BlockVector sg_dx_kl_block(*base_x_map, *sg_x_kl_map);
    AztecOO kl_aztec;
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
    aztec.SetAztecOption(AZ_precond, AZ_none);
    aztec.SetAztecOption(AZ_kspace, 20);
    aztec.SetAztecOption(AZ_conv, AZ_r0);
    aztec.SetAztecOption(AZ_output, 1);
    aztec.SetUserOperator(&kl_jac_op);
    aztec.SetPrecOperator(&kl_prec_op);
    aztec.SetLHS(&sg_dx_kl_block);
    aztec.SetRHS(&sg_f_kl_block);
    aztec.Iterate(70, 1e-4);

    // Map solution back to original basis
    utils.out() << "Mapping solution back to original basis..." << std::endl;
    const Teuchos::Array<double>& norms = basis->norm_squared();
    EpetraExt::BlockVector sg_ov_dx_kl_block(*base_x_map, *sg_ov_x_kl_map);
    sg_ov_dx_kl_block.Import(sg_dx_kl_block, sg_x_kl_importer, Insert);
    Teuchos::RCP<EpetraExt::BlockVector> sg_ov_dx_block = 
      sg_model->import_solution(*sg_x);
    sg_ov_dx_block->PutScalar(0.0);
    Stokhos::EpetraVectorOrthogPoly sg_dx_poly2(
      basis, stoch_overlap_map, base_x_map, 
      Teuchos::rcp(&(sg_ov_dx_block->Map()), false), multi_comm, 
      View, *sg_ov_dx_block);
    Stokhos::EpetraVectorOrthogPoly sg_dx_kl_poly(
      kl_basis, kl_ov_stoch_row_map, base_x_map, sg_ov_x_kl_map, 
      kl_comm, View, sg_ov_dx_kl_block);
    Epetra_Vector dx_val(*base_x_map);
    for (int qp=0; qp < nqp; qp++) {
      sg_dx_kl_poly.evaluate(kl_basis_vals[qp], dx_val);
      sg_dx_poly2.sumIntoAllTerms(weights[qp], basis_vals[qp], norms, dx_val);
    }
    Teuchos::RCP<EpetraExt::BlockVector> sg_dx_block2 =
      sg_model->export_solution(*sg_ov_dx_block);

    // Epetra_Vector sg_dx_err(*sg_dx);
    // sg_dx_err.Update(1.0, *sg_dx, -1.0, sg_dx_block2, 0.0);
    // sg_dx_err.Print(utils.out());

    sg_x->Update(-1.0, *sg_dx_block2, 1.0);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    sg_f->Norm2(&normf);
    sg_dx_block2->Norm2(&norm_dx);
    nit++;
    utils.out() << "nit = " << nit << " norm_f = " << normf << " norm_dx = " 
	      << norm_dx << std::endl;

    }

    EpetraExt::ModelEvaluator::OutArgs sg_outArgs2 = 
	sg_model->createOutArgs();
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(0))));
    sg_outArgs2.set_g(0, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs2);

    // Print mean and standard deviation
    utils.out().precision(12);
    sg_g->Print(utils.out());
    double mean = (*sg_g)[0];
    double std_dev = 0.0;
    const Teuchos::Array<double>& nrm2 = basis->norm_squared();
    for (int i=1; i<basis->size(); i++)
      std_dev += (*sg_g)[i]*(*sg_g)[i]*nrm2[i];
    std_dev = std::sqrt(std_dev);
    
    utils.out() << "Mean =      " << mean << std::endl;
    utils.out() << "Std. Dev. = " << std_dev << std::endl;

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif

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

}
