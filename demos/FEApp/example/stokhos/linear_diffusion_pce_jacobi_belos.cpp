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

// AztecOO solver
#include "AztecOO.h"

// Belos solver
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
//#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
//#include "BelosGCRODRSolMgr.hpp"

//Ifpack preconditioner
#include "Ifpack.h"

// Stokhos Stochastic Galerkin
#include "Stokhos.hpp"
#include "EpetraExt_BlockVector.h"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) {
  int nelem = 100;
  double h = 1.0/nelem;
  int num_KL = 3;
  int p = 5;
  bool full_expansion = true;


  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;


// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  int MyPID;

  try {

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create a communicator for Epetra objects
    Teuchos::RCP<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    MyPID = Comm->MyPID();
    
    // Create mesh
    std::vector<double> x(nelem+1);
    for (int i=0; i<=nelem; i++)
      x[i] = h*i;

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

    // Free parameters (determinisic, e.g., for sensitivities)
    Teuchos::RefCountPtr< Teuchos::Array<std::string> > free_param_names =
	Teuchos::rcp(new Teuchos::Array<std::string>);
    free_param_names->push_back("Constant Source Function Value");
    
    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (full_expansion)
      Cijk = basis->computeTripleProductTensor(sz);
    else
      Cijk = basis->computeTripleProductTensor(num_KL+1);
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis,
									 Cijk));
    std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create application
    appParams->set("SG Method", "AD");
    Teuchos::RCP<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(x, Comm, appParams, false));
    
    // Set up stochastic parameters
    Epetra_LocalMap p_sg_map(num_KL, 0, *Comm);
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p_init = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis, p_sg_map));
    for (int i=0; i<num_KL; i++) {
      sg_p_init->term(i,0)[i] = 0.0;
      sg_p_init->term(i,1)[i] = 1.0;
    }
    Teuchos::RefCountPtr< Teuchos::Array<std::string> > sg_param_names =
      Teuchos::rcp(new Teuchos::Array<std::string>);
    for (int i=0; i<num_KL; i++) {
      std::stringstream ss;
      ss << "KL Exponential Function Random Variable " << i;
      sg_param_names->push_back(ss.str());
    }

    // Setup stochastic initial guess
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis, 
						       *(app->getMap())));
    sg_x_init->init(0.0);
    
    // Create application model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names,
					     sg_param_names, sg_x_init, 
					     sg_p_init));
    
    // Setup stochastic Galerkin algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(&(appParams->sublist("SG Parameters")),false);
    if (!full_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    sgParams->set("Jacobian Method", "Matrix Free");
    sgParams->set("Mean Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      sgParams->sublist("Preconditioner Parameters");
    precParams.set("default values", "DD");

    // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
						 expansion, Cijk, sgParams,
						 Comm));

    // Create vectors and operators
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(2);
    Teuchos::RCP<Epetra_Vector> sg_x = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    Teuchos::RCP<Epetra_Vector> sg_y = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    *sg_x = *(sg_model->get_x_init());
    Teuchos::RCP<Epetra_Vector> sg_f = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_f_map())));
    Teuchos::RCP<Epetra_Vector> sg_df = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_f_map())));
    Teuchos::RCP<Epetra_Vector> sg_dx = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_x_map())));
    Teuchos::RCP<Epetra_Operator> sg_J = sg_model->create_W();
    Teuchos::RCP<Epetra_Operator> sg_M = sg_model->create_WPrec()->PrecOp;
    
    // std::cout << "parameter PC expansion:" << std::endl 
//            << *sg_p_init << std::endl;
  //  std::cout << "parameter block vector:" << std::endl 
//            << *sg_p << std::endl;
//
  //  std::cout << "parameter PC expansion:" << std::endl 
//            << *sg_x_init << std::endl;
  //  std::cout << "parameter block vector:" << std::endl 
//            << *sg_x << std::endl;

  //  EpetraExt::BlockVector sg_p_block(View, p_sg_map, *sg_p);
  //  Teuchos::RCP<const Epetra_Vector> sg_p_vec = sg_p_block.GetBlock(2);

//    EpetraExt::BlockVector sg_x_block(View, *(app->getMap()), *sg_x);
  //  Teuchos::RCP<const Epetra_Vector> sg_x_vec = sg_x_block.GetBlock(2);


    // Setup InArgs and OutArgs
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = sg_model->createOutArgs();
    sg_inArgs.set_p(2, sg_p);
    sg_inArgs.set_x(sg_x);
    sg_outArgs.set_f(sg_f);
    sg_outArgs.set_W(sg_J);
    sg_outArgs.set_WPrec(sg_M);

    // Evaluate model
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    Teuchos::RCP<Stokhos::MatrixFreeEpetraOp> stokhos_op =
      Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(sg_J, true);
    Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_J_poly
      = stokhos_op->getOperatorBlocks();
    // to get k-th matrix:  (*sg_J_poly)[k]

    Teuchos::RCP<Stokhos::MeanEpetraOp> mean_op = 
      Teuchos::rcp_dynamic_cast<Stokhos::MeanEpetraOp>(sg_M, true);

    Teuchos::RCP<Epetra_Operator> mean_prec
      = mean_op->getMeanOperator();

    //const Stokhos::VectorOrthogPoly<Epetra_Operator>& sg_J_poly_ref = *sg_J_poly;

    std::vector< Teuchos::RCP<const Epetra_Vector> > sg_p_vec_all ;
    std::vector< Teuchos::RCP< Epetra_MultiVector> > sg_x_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_dx_vec_all ; 
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_f_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_df_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_kx_vec_all ;

    // Extract blocks
    EpetraExt::BlockVector sg_p_block(View, p_sg_map, *sg_p);
    EpetraExt::BlockVector sg_x_block(View, *(app->getMap()), *sg_x);
    EpetraExt::BlockVector sg_dx_block(View, *(app->getMap()), *sg_dx);
    EpetraExt::BlockVector sg_f_block(View, *(app->getMap()), *sg_f);

    // sg_p_vec_all.push_back(sg_p_block.GetBlock(0));
    for (int i=0; i<sz; i++) {
      sg_p_vec_all.push_back(sg_p_block.GetBlock(i));
      sg_x_vec_all.push_back(sg_x_block.GetBlock(i));
      sg_dx_vec_all.push_back(sg_dx_block.GetBlock(i));
      sg_f_vec_all.push_back(sg_f_block.GetBlock(i));
      sg_df_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(*(app->getMap()))));
      sg_kx_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(*(app->getMap()))));
    }

//    Teuchos::RCP<Epetra_MultiVector> sg_x_mv =
//      Teuchos::rcp_dynamic_cast<Epetra_MultiVector>(sg_x_block, true);    

    //create epetra multivector
    Teuchos::RCP<Epetra_MultiVector> sg_x_mv = 
	Teuchos::rcp(new Epetra_MultiVector(*(app->getMap()),sz));
    Teuchos::RCP<Epetra_MultiVector> sg_f_mv = 
	Teuchos::rcp(new Epetra_MultiVector(*(app->getMap()),sz));
    Teuchos::RCP<Epetra_MultiVector> sg_dx_mv = 
	Teuchos::rcp(new Epetra_MultiVector(*(app->getMap()),sz));

    std::vector<int> vecind(1);
    for (int i=0; i<sz; i++) {
      vecind[0] = i;
      MVT::SetBlock(*(sg_x_vec_all[i]), vecind, *sg_x_mv);
      MVT::SetBlock(*(sg_f_vec_all[i]), vecind, *sg_f_mv);
      MVT::SetBlock(*(sg_dx_vec_all[i]), vecind, *sg_dx_mv);
    }

    Teuchos::RCP<Epetra_CrsMatrix> sg_J_poly_Crs =
        Teuchos::rcp_dynamic_cast< Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(0),true);

  bool proc_verbose = true;  
  bool leftprec = true;      // left preconditioning or right.
  int numrhs = sz;
  //
  // ************Construct preconditioner*************
  //
  ParameterList ifpackList;

  // allocates an IFPACK factory. No data is associated
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  std::string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.
 
  Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType,&*sg_J_poly_Crs, OverlapLevel) );
  assert(Prec != Teuchos::null);
  
  // specify parameters for ILU
  ifpackList.set("fact: drop tolerance", 1e-10);
  ifpackList.set("fact: level-of-fill", 1);
  // the combine mode is on the following:
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  // Their meaning is as defined in file Epetra_CombineMode.h
  ifpackList.set("schwarz: combine mode", "Add");
  // sets the parameters
  IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));
  
  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec->Initialize());
  
  // Builds the preconditioners, by looking for the values of
  // the matrix.
  IFPACK_CHK_ERR(Prec->Compute());
  
  // Create the Belos preconditioned operator from the Ifpack preconditioner.
  // NOTE:  This is necessary because Belos expects an operator to apply the
  //        preconditioner with Apply() NOT ApplyInverse().
  Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( Prec ) );

/*    Teuchos::RCP< Belos::LinearProblem<double,MV,OP> > myProblem 
	= Teuchos::rcp( new Belos::LinearProblem<double,MV,OP>(sg_J_poly_Crs, sg_dx_mv, sg_f_mv) );

//    myProblem->setHermitian();
    myProblem->setProblem();
*/
    int verbosity = Belos::Warnings + Belos::Errors + Belos::FinalSummary + Belos::TimingDetails;

    // Create the parameter list for the eigensolver
    Teuchos::RCP<Teuchos::ParameterList> myPL = Teuchos::rcp( new Teuchos::ParameterList() );
    myPL->set( "Verbosity", verbosity );
    myPL->set( "Block Size", 1 );
    myPL->set( "Convergence Tolerance", 1.0e-10 );

  myPL->set( "Num Blocks", 101 );            // Maximum number of blocks in Krylov factorization
  myPL->set( "Maximum Iterations", 1000 );       // Maximum number of iterations allowed

//  myPL->set( "Maximum Restarts", 50 );      // Maximum number of restarts allowed
//  myPL->set( "Convergence Tolerance", 1.0e-10 );         // Relative convergence tolerance requested
//  myPL->set( "Num Recycled Blocks", 50 );       // Number of vectors in recycle space
//  myPL->set( "Orthogonalization", "IMGS" );           // Orthogonalization type
  if (numrhs > 1) {
    myPL->set( "Show Maximum Residual Norm Only", true );  // Show only the maximum residual norm 
  }


/*    // Create the Block CG solver
    // This takes as inputs the linear problem and the solver parameters
    Belos::BlockGmresSolMgr<double,MV,OP> mySolver(myProblem, myPL);

    // Solve the linear problem, and save the return code

    // Solve the linear problem, and save the return code
    Belos::ReturnType solverRet = mySolver.solve(); 

*/
    
//    std::cout << "residula solution:" << std::endl
  //            << *sg_dx_mv << std::endl;
//    std::cout << "parameter block vector:" << std::endl
  //            << ((*sg_J_poly).getCoefficients())[0] << std::endl;
//    std::cout << "Block Matrices:" << std::endl
  //            << (*sg_J_poly)[0] << std::endl;

   // std::cout << "force block vector:" << std::endl
     //         << *sg_f_vec_all[0] << std::endl;

    // Print initial residual norm
    double norm_f,norm_df;
    sg_f->Norm2(&norm_f);
    sg_J->Apply(*(sg_dx),*(sg_y));
    sg_df->Update(1.0,*sg_y,-1.0,*(sg_f),0.0);
    sg_df->Norm2(&norm_df);

    std::cout << "\nInitial residual norm = " << norm_df/norm_f << std::endl;
   
    const Teuchos::Array<double>& norms = sg_J_poly->basis()->norm_squared(); 
//  Teuchos::RCP< Epetra_Vector> kx ;
//    Teuchos::RCP<Epetra_Vector> kx =
//      Teuchos::rcp(new Epetra_Vector(*(app->getMap())));
    Teuchos::RCP<Epetra_MultiVector> dx =
      Teuchos::rcp(new Epetra_MultiVector(*(app->getMap()),1));
    Teuchos::RCP<Epetra_Vector> df =
      Teuchos::rcp(new Epetra_Vector(*(app->getMap())));

   // (*sg_J_poly)[0].Apply(*(sg_f_vec_all[0]),*kx);
   // std::cout << "f(0):" << std::endl
   //           << *(sg_f_vec_all[0]) << std::endl;
   // std::cout << "kx:" << std::endl
  //            << *kx << std::endl;

sg_dx->PutScalar(0.0);

int iter = 0;
//for (int iter=0;iter<60;iter++){
while ((norm_df/norm_f)>1e-12) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++; 
     // Extract blocks
   // EpetraExt::BlockVector sg_f_block(View, *(app->getMap()), *sg_f);

    // sg_p_vec_all.push_back(sg_p_block.GetBlock(0));
//    for (int i=0; i<sz; i++) {
  //   vecind[0] = i;
    // sg_f_vec_all.push_back(sg_f_block.GetBlock(i));
    //  MVT::SetBlock(*(sg_f_vec_all[i]), vecind, *sg_f_mv);
    // }

    for (int i=0; i<sz; i++) {
      (sg_df_vec_all[i])->Update(1.0, *sg_f_vec_all[i], 0.0);
    }
//      sg_df->Update(1.0, *sg_f, 0.0);

//    double c0kk;
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
  //  ordinal_type Cijk_size = Cijk.size();
    for (int k=1; k<num_KL+1; k++) {
    //  df->Update(1.0, *sg_f_vec_all[k], 0.0);
      int nj = Cijk->num_j(k);
      const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
      //std::cout << "j_indices for k =" << k << j_indices << std::endl;
      for (int jj=0; jj<nj; jj++) {
        int j = j_indices[jj];
//        vecind[0] = j;
  //      sg_dx_vec_all[j] = MVT::CloneCopy(*sg_dx_mv, vecind);
        (*sg_J_poly)[k].Apply(*(sg_dx_vec_all[j]),*(sg_kx_vec_all[j]));
      }
      for (int jj=0; jj<nj; jj++) {
        int j = j_indices[jj];
        const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
        const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
        int ni = i_indices.size();
        for (int ii=0; ii<ni; ii++) {
          int i = i_indices[ii];
          double c = cijk_values[ii];  // C(i,j,k)
          sg_df_vec_all[i]->Update(-1.0*c,*(sg_kx_vec_all[j]),1.0);
    /*      if (i !=0 && i<num_KL+1) {
           (*sg_J_poly)[i].Apply(*dx,*kx);
           df->Update(-1.0*c,*kx,1.0);          
          }
          if (i==0 && j==k) {
            c0kk = c;
	    // std::cout << "c0kk in [" << i << j << k << "]th iteration = " << c0kk << std::endl;
          }*/
        }
      }    


//      aztec.SetRHS(df.get());
      // Solve linear system
//      aztec.Iterate(100, 1e-12);

      // Update x
//      sg_dx_vec_all[k]->Update(1.0, *dx, 0.0);

    } //End of k loop

    for(int i=0; i<sz; i++) {    
      sg_df_vec_all[i]->Scale(1/norms[i]);
      vecind[0] = i;
      MVT::SetBlock(*(sg_df_vec_all[i]), vecind, *sg_f_mv);  
    }

    Teuchos::RCP< Belos::LinearProblem<double,MV,OP> > myProblem 
        = Teuchos::rcp( new Belos::LinearProblem<double,MV,OP>(sg_J_poly_Crs, sg_dx_mv, sg_f_mv) );

//    myProblem->setHermitian();
/*  if (leftprec) {
    myProblem->setLeftPrec( belosPrec );
  }
  else {
    myProblem->setRightPrec( belosPrec );
  }*/
  myProblem->setLeftPrec( belosPrec );
  myProblem->setProblem();
/*  bool set = myProblem->setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
*/

    //int verbosity = Belos::Warnings + Belos::Errors + Belos::FinalSummary + Belos::TimingDetails;

/*    // Create the parameter list for the eigensolver
    Teuchos::RCP<Teuchos::ParameterList> myPL = Teuchos::rcp( new Teuchos::ParameterList() );
    myPL->set( "Verbosity", verbosity );
    myPL->set( "Block Size", 1 );
    myPL->set( "Maximum Iterations", 1000 );
    myPL->set( "Convergence Tolerance", 1.0e-10 );
  */  
    // Create the Block Gmres solver
    // This takes as inputs the linear problem and the solver parameters
    Belos::BlockGmresSolMgr<double,MV,OP> mySolver(myProblem, myPL);
//    Belos::GCRODRSolMgr<double,MV,OP> mySolver(myProblem, myPL);
    

    // Solve the linear problem, and save the return code
//    Belos::ReturnType solverRet = mySolver.solve();
   {
    TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");
    mySolver.solve();
   }

    for (int k=0; k<sz; k++) {
      vecind[0] = k;
      dx = MVT::CloneCopy(*sg_dx_mv, vecind);
      (sg_dx_vec_all[k])->Update(1.0, *dx, 0.0);
    }

    sg_J->Apply(*(sg_dx),*(sg_y));
    sg_df->Update(1.0,*sg_y,-1.0,*sg_f,0.0);
    sg_df->Norm2(&norm_df);

    std::cout << "rel residual norm at iteration "<< iter <<" is " << norm_df/norm_f << std::endl;

  } //End of iter loop 

  sg_x_mv->Update(-1.0, *sg_dx_mv, 1.0);

  for (int k=0; k<sz; k++) {
//    vecind[0] = k;
//    sg_dx_vec_all[k] = MVT::CloneCopy(*sg_dx_mv, vecind);
    sg_x_vec_all[k]->Update(-1.0, *sg_dx_vec_all[k], 1.0);
  }

  //  std::cout << "f(0):" << std::endl
    //          << *df << std::endl;
//    std::cout << "x(0):" << std::endl
//              << *(sg_x_vec_all[0]) << std::endl;

    //  int LoadBlockValues(const Epetra_Vector & BaseVec, int BlockRow)
    //    sg_x_block->LoadBlockValues(*(sg_x_vec_all[0]),0);

    std::cout << "sg_x_mv" << *sg_x << std::endl;

//    sg_inArgs.set_x(sg_x);

    // Compute new residual & response function
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(1))));
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs2 = sg_model->createOutArgs();
    sg_outArgs2.set_f(sg_f);
    sg_outArgs2.set_g(1, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs2);

    // Print Final residual norm
    sg_f->Norm2(&norm_f);
    std::cout << "\nFinal residual norm = " << norm_f << std::endl;

    // Print mean and standard deviation
    Stokhos::EpetraVectorOrthogPoly sg_g_poly(basis, View, 
					      *(model->get_g_map(0)), *sg_g);
    Epetra_Vector mean(*(model->get_g_map(0)));
    Epetra_Vector std_dev(*(model->get_g_map(0)));
    sg_g_poly.computeMean(mean);
    sg_g_poly.computeStandardDeviation(std_dev);
    std::cout << "\nResponse Expansion = " << std::endl;
    std::cout.precision(12);
    sg_g_poly.print(std::cout);
    std::cout << "\nResponse Mean =      " << std::endl << mean << std::endl;
    std::cout << "Response Std. Dev. = " << std::endl << std_dev << std::endl;

    if (norm_f < 1.0e-10)
      std::cout << "Test Passed!" << std::endl;

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
