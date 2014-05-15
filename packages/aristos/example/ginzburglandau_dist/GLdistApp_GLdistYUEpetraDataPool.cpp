//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Amesos.h"
#include "Ifpack.h"
#include "Amesos_Klu.h"
#include "Amesos_Umfpack.h"
// #include "Amesos_Superludist.h"
// #include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_Reindex_LinearProblem.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "GLdistApp_SchurOp.hpp"
#include "GLdistApp_GLdistYUEpetraDataPool.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <fstream>

namespace GLdistApp {

bool GLdistYUEpetraDataPool::useStratimikos = false;
std::string GLdistYUEpetraDataPool::stratimikosXmlFile = "";

GLdistYUEpetraDataPool::GLdistYUEpetraDataPool( Epetra_Comm * commptr,
  double beta, const std::string &myfile )
  : commptr_(commptr),
    beta_(beta)
{
  ipcoords_ = Teuchos::rcp( new Epetra_SerialDenseMatrix() );
  ipindx_ = Teuchos::rcp( new Epetra_IntSerialDenseVector() );
  pcoords_ = Teuchos::rcp( new Epetra_SerialDenseMatrix() );
  pindx_ = Teuchos::rcp( new Epetra_IntSerialDenseVector() );
  t_ = Teuchos::rcp( new Epetra_IntSerialDenseMatrix() );
  e_ = Teuchos::rcp( new Epetra_IntSerialDenseMatrix() );

  strcpy(geomfile_, myfile.c_str());  

  // Read subdomain info.
  meshreader(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_,
    myfile.c_str());

  // Assemble volume and boundary mass and stiffness matrices, and the right-hand side of the PDE.
  assemble(*commptr, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, A_, H_, b_);
  assemble_bdry(*commptr, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, *e_, B_, R_);

  // Set desired state q.
  Epetra_Map standardmap(A_->DomainMap());
  q_ = Teuchos::rcp(new Epetra_FEVector(standardmap));
  int * qintvalues = new int[standardmap.NumMyElements()];
  double * qdvalues = new double[standardmap.NumMyElements()];
  standardmap.MyGlobalElements(qintvalues);
  for (int i = 0; i < standardmap.NumMyElements(); i++)
      qdvalues[i]=cos( M_PI* ((*ipcoords_)(i,0)) ) * cos( M_PI* ((*ipcoords_)(i,1)) );
  q_->ReplaceGlobalValues(standardmap.NumMyElements(), qintvalues, qdvalues);
  q_->GlobalAssemble();
}

void GLdistYUEpetraDataPool::computeAll( const Aristos::Vector &x )
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> ey =
        (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getYVector();

  computeNy(ey);

  computeNpy(ey);

  computeAugmat();

  computePrec();
  
}


int GLdistYUEpetraDataPool::solveAugsys( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                                       const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                                       const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                                       const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                                       const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                                       const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                                       double * tol )
{
  int systemChoice = 1;   // 1 for full KKT system solve, 2 for Schur complement solve

  // These options are for the full KKT system solve.
  // 11 for AztecOO with built-in Schwarz DD preconditioning and ILU on subdomains
  // 12 for AztecOO with IFPACK Schwarz DD preconditioning and Umfpack on subdomains
  // 13 for a direct sparse solver (Umfpack, KLU)
  int solverChoice = 12;

  std::cout << "\nsolverChoice = " << solverChoice << "\n";
  
  if (systemChoice == 1) {
    // We're using the full KKT system formulation to solve the augmented system.
   
    Epetra_Map standardmap(A_->DomainMap());
    Epetra_Map bdryctrlmap(B_->DomainMap());
    Epetra_Vector rhs( (Epetra_BlockMap&)Augmat_->RangeMap() );
    Epetra_Vector soln( (Epetra_BlockMap&)Augmat_->RangeMap() );
    soln.PutScalar(1.0);  

    Teuchos::Array<double> values(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
    Teuchos::Array<int> indices(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
    ((Epetra_BlockMap&)Augmat_->RangeMap()).MyGlobalElements(indices.getRawPtr());

    for (int i=0; i<rhsy->MyLength(); i++) {
      values[i] = (*((*rhsy)(0)))[i];
    }
    for (int i=0; i<rhsu->MyLength(); i++) {
      values[i+rhsy->MyLength()] = (*((*rhsu)(0)))[i];
    }
    for (int i=0; i<rhsp->MyLength(); i++) {
      values[i+rhsy->MyLength()+rhsu->MyLength()] = (*((*rhsp)(0)))[i];
    }

    rhs.ReplaceGlobalValues(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength(), values.getRawPtr(), indices.getRawPtr());

    if (solverChoice == 11) {
      int Overlap = 3;
      int ival = 4;

      AztecOO::AztecOO kktsolver(&(*Augmat_), &soln, &rhs);
      kktsolver.SetAztecOption( AZ_solver, AZ_gmres );
      kktsolver.SetAztecOption( AZ_precond, AZ_dom_decomp );
      kktsolver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
      //kktsolver.SetAztecOption( AZ_kspace, 2*numstates+numcontrols );
      kktsolver.SetAztecOption( AZ_kspace, 1000 );
      kktsolver.SetAztecOption(AZ_overlap,Overlap);
      kktsolver.SetAztecOption(AZ_graph_fill,ival);
      //kktsolver.SetAztecOption(AZ_poly_ord, ival);
      //kktsolver.SetAztecParam(AZ_drop, 1e-9);
      kktsolver.SetAztecParam(AZ_athresh, 1e-5);
      //kktsolver.SetAztecParam(AZ_rthresh, 0.0);
      kktsolver.SetAztecOption( AZ_reorder, 0 );
      //kktsolver.SetAztecParam44( AZ_ilut_fill, 1.5 );
      kktsolver.SetAztecOption( AZ_output, AZ_last );
      //kktsolver.Iterate(2*numstates+numcontrols,1e-12);
      kktsolver.Iterate(1000,1e-11);
      //cout << soln;
    }
    else if (solverChoice == 12) {
      // =============================================================== //
      // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
      // =============================================================== //

      Teuchos::ParameterList List;

      // allocates an IFPACK factory. No data is associated
      // to this object (only method Create()).
      Ifpack Factory;

      // create the preconditioner. For valid PrecType values,
      // please check the documentation
      string PrecType = "Amesos";
      int OverlapLevel = 2; // must be >= 0. If Comm.NumProc() == 1,
                            // it is ignored.
  
      Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &(*Augmat_), OverlapLevel);
      assert(Prec != 0);

      // specify the Amesos solver to be used.
      // If the selected solver is not available,
      // IFPACK will try to use Amesos' KLU (which is usually always
      // compiled). Amesos' serial solvers are:
      // "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
      // List.set("amesos: solver type", "Amesos_Umfpack");
      List.set("amesos: solver type", "Amesos_Klu");

      // sets the parameters
      IFPACK_CHK_ERR(Prec->SetParameters(List));

      // initialize the preconditioner. At this point the matrix must
      // have been FillComplete()'d, but actual values are ignored.
      // At this call, Amesos will perform the symbolic factorization.
      IFPACK_CHK_ERR(Prec->Initialize());

      // Builds the preconditioners, by looking for the values of
      // the matrix. At this call, Amesos will perform the
      // numeric factorization.
      IFPACK_CHK_ERR(Prec->Compute());

      // =================================================== //
      // E N D   O F   I F P A C K   C O N S T R U C T I O N //
      // =================================================== //

      // need an Epetra_LinearProblem to define AztecOO solver
      Epetra_LinearProblem Problem;
      Problem.SetOperator(&(*Augmat_));
      Problem.SetLHS(&soln);
      Problem.SetRHS(&rhs);

      // now we can allocate the AztecOO solver
      AztecOO kktsolver(Problem);

      // specify solver
      kktsolver.SetAztecOption(AZ_solver,AZ_gmres);
      kktsolver.SetAztecOption(AZ_kspace, 300 );
      kktsolver.SetAztecOption(AZ_output,AZ_none);

      // HERE WE SET THE IFPACK PRECONDITIONER
      kktsolver.SetPrecOperator(Prec);

      // .. and here we solve
      kktsolver.Iterate(300,1e-12);

      // delete the preconditioner
      delete Prec;
    }
    else if (solverChoice == 13) {
      Epetra_LinearProblem Problem;
      Problem.SetOperator(&(*Augmat_));
      Problem.SetLHS(&soln);
      Problem.SetRHS(&rhs);
      
      EpetraExt::LinearProblem_Reindex reindex(NULL);
      Epetra_LinearProblem newProblem = reindex(Problem);
      
      //Amesos_Umfpack kktsolver(newProblem);
      Amesos_Klu kktsolver(newProblem);
   
      AMESOS_CHK_ERR(kktsolver.SymbolicFactorization());
      AMESOS_CHK_ERR(kktsolver.NumericFactorization());
      AMESOS_CHK_ERR(kktsolver.Solve());
      kktsolver.PrintTiming();
    }
    
    
    for (int i=0; i<rhsy->MyLength(); i++) {
      (*((*y)(0)))[i] = soln[i];
    }
    for (int i=0; i<rhsu->MyLength(); i++) {
      (*((*u)(0)))[i] = soln[i+rhsy->MyLength()];
    }
    for (int i=0; i<rhsp->MyLength(); i++) {
      (*((*p)(0)))[i] = soln[i+rhsy->MyLength()+rhsu->MyLength()];
    }
    
  }
  else if (systemChoice == 2) {
    // We're using the Schur complement formulation to solve the augmented system.
  
    // Form linear operator.
    GLdistApp::SchurOp schurop(A_, B_, Npy_);
  
    // Form Schur complement right-hand side.
    Epetra_MultiVector ny( (Epetra_BlockMap&)Npy_->RangeMap(), 1);
    Epetra_MultiVector ay( (Epetra_BlockMap&)A_->RangeMap(), 1);
    Epetra_MultiVector schurrhs( (Epetra_BlockMap&)A_->RangeMap(), 1);
    Epetra_MultiVector bu( (Epetra_BlockMap&)B_->RangeMap(), 1);
    A_->Multiply(false, *rhsy, ay);
    Npy_->Multiply(false, *rhsy, ny);
    B_->Multiply(false, *rhsu, bu);
    schurrhs.Update(1.0, ny, 1.0, ay, 0.0);
    schurrhs.Update(-1.0, *rhsp, 1.0, bu, 1.0);
  
    p->PutScalar(0.0);
    Epetra_LinearProblem linprob(&schurop, &(*p), &schurrhs);
    AztecOO::AztecOO Solver(linprob);
    Solver.SetAztecOption( AZ_solver, AZ_cg );
    Solver.SetAztecOption( AZ_precond, AZ_none );
    Solver.SetAztecOption( AZ_output, AZ_none );
    Solver.Iterate(8000,tol[0]);
  
    Epetra_MultiVector bp( (Epetra_BlockMap&)B_->DomainMap(), 1);
    B_->Multiply(true, *p, bp);
    u->Update(1.0, *rhsu, -1.0, bp, 0.0);

    Epetra_MultiVector ap( (Epetra_BlockMap&)A_->DomainMap(), 1);
    Epetra_MultiVector np( (Epetra_BlockMap&)A_->DomainMap(), 1);
    A_->Multiply(true, *p, ap);
    Npy_->Multiply(true, *p, np);
    y->Update(1.0, *rhsy,0.0);
    y->Update(-1.0, ap, -1.0, np, 1.0);
  }
  
  return 0;
}


int GLdistYUEpetraDataPool::solveAugsysDyn( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                                         const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                                         const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                                         const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                                         const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                                         const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                                         double * tol )
{

  const RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //*out << "\nsolveAugsysDyn(...) ...\n";

  Teuchos::OSTab tab(out);

  double minthreshold = 1e-20;  // tolerance threshold for "exact" solves, anything below is computed w/ minprecision
  double minprecision = 1e-12;  // relative tolerance for "exact" solves
  int    maxit        = 200;    // maximum number of Krylov solver iterations
  int    innerit      = 50;     // number of Krylov solver iterations, per restart, for the first dynamic solve
  int    numcycles    = 4;      // number of restarts for first dynamic solve
  bool   wantstats    = false;  // choose true if output of solver info in file stats.txt is desired
  int mypid = y->Comm().MyPID();

  ofstream outfile("stats.txt", ios_base::app);

  double tolerance = tol[0];

  // will need to set up state-control-adjoint splitting and corresponding maps

  Epetra_Map standardmap(A_->DomainMap());
  Epetra_Map bdryctrlmap(B_->DomainMap());
  Epetra_Vector rhs( (Epetra_BlockMap&)Augmat_->RangeMap() );
  Epetra_Vector soln( (Epetra_BlockMap&)Augmat_->RangeMap() );
  // Set initial iterate.
  soln.PutScalar(0.0);  

  Teuchos::Array<double> values(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
  Teuchos::Array<int> indices(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength());
  ((Epetra_BlockMap&)Augmat_->RangeMap()).MyGlobalElements(indices.getRawPtr());

  for (int i=0; i<rhsy->MyLength(); i++) {
    values[i] = (*((*rhsy)(0)))[i];
  }
  for (int i=0; i<rhsu->MyLength(); i++) {
    values[i+rhsy->MyLength()] = (*((*rhsu)(0)))[i];
  }
  for (int i=0; i<rhsp->MyLength(); i++) {
    values[i+rhsy->MyLength()+rhsu->MyLength()] = (*((*rhsp)(0)))[i];
  }

  rhs.ReplaceGlobalValues(rhsy->MyLength() + rhsu->MyLength() + rhsp->MyLength(), values.getRawPtr(), indices.getRawPtr());


/*
  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // allocates an IFPACK factory. No data is associated
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  string PrecType = "Amesos";
  int OverlapLevel = 4; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.
  
  Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &(*Augmat_), OverlapLevel);
  assert(Prec != 0);

  // specify the Amesos solver to be used.
  // If the selected solver is not available,
  // IFPACK will try to use Amesos' KLU (which is usually always
  // compiled). Amesos' serial solvers are:
  // "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
  List.set("amesos: solver type", "Amesos_Umfpack");

  // sets the parameters
  IFPACK_CHK_ERR(Prec->SetParameters(List));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  // At this call, Amesos will perform the symbolic factorization.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Builds the preconditioners, by looking for the values of
  // the matrix. At this call, Amesos will perform the
  // numeric factorization.
  IFPACK_CHK_ERR(Prec->Compute());

  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //
*/



/*
  // =============================== //
  // Create the ML + IFPACK smoother //
  // =============================== //

  Teuchos::ParameterList MLList;
                                                                                
  ML_Epetra::SetDefaults("DD",MLList);
  MLList.set("smoother: pre or post", "post");
  MLList.set("max levels", 2);
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("aggregation: nodes per aggregate", 9);
  MLList.set("PDE equations", 1);

  // fix the smoother to be IFPACK; can be set using (level X) syntax
  MLList.set("smoother: type","IFPACK");

  // now we have to specify which IFPACK preconditioner should be
  // built. Any value that is valid for the IFPACK factory. We also need
  // to define the overlap (>= 0).
                                                                                
  MLList.set("smoother: ifpack type", "Amesos");
  MLList.set("smoother: ifpack overlap", 4);

  // Then, all parameters can will control the definition of the IFPACK
  // smoother are inserted in IFPACKList. In this case, we specify the fill-in
  // factor. For a list of supported parameters, please consult the IFPACK
  // documentation. For example, IFPACK preconditioner "Amesos" or
  // "Amesos stand-alone" can be used to solve with an LU
  // factorization on each domain.
  MLList.sublist("smoother: ifpack list").set("amesos: solver type", "Amesos_Umfpack");

  // we can now build the preconditioner...
                                                                                
  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner( (*Augmat_), MLList, true);

  // =============================== //
  //   End of ML + IFPACK smoother   //
  // =============================== //
*/

  // Solve the system

  if (useStratimikos) {
    
    *out << "\nUsing Stratimikos!\n";

    using Teuchos::RCP; using Teuchos::FancyOStream; using Teuchos::describe;
    using Teuchos::rcpFromRef;

    // Wrap the Epetra objects as Thyra objects

    RCP<Thyra::VectorBase<double> > thyra_soln =
      Thyra::create_Vector( rcpFromRef(soln), thyra_Augmat_->domain() );
    RCP<const Thyra::VectorBase<double> > thyra_rhs =
      Thyra::create_Vector( rcpFromRef(rhs), thyra_Augmat_->range() );

    // Create the SolveCriteria and solve

    Thyra::SolveCriteria<double> solveCriteria;
    solveCriteria.solveMeasureType.numerator = Thyra::SOLVE_MEASURE_NORM_RESIDUAL;
    //solveCriteria.numeratorReductionFunc = createMockNormReductionFunctional<Scalar>();
    solveCriteria.solveMeasureType.denominator = Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL;
    //solveCriteria.denominatorReductionFunc = createMockMaxNormInfEpsReductionFunctional<Scalar>();
    solveCriteria.requestedTol = 1e-6;

    Thyra::SolveStatus<double> status = Augmat_lows_->solve(
      Thyra::NOTRANS, *thyra_rhs, thyra_soln.ptr(),
      optInArg(solveCriteria));
    *out << "\nSolve status:\n" << status;

  }
  else {

    *out << "\nUsing AztecOO!\n";

    // need an Epetra_LinearProblem to define AztecOO solver
    Epetra_LinearProblem Problem;
    Problem.SetOperator(&(*Augmat_));
    Problem.SetLHS(&soln);
    Problem.SetRHS(&rhs);

    // now we can allocate the AztecOO solver
    AztecOO kktsolver(Problem);

    // HERE WE SET THE IFPACK PRECONDITIONER
    //kktsolver.SetPrecOperator(Prec);
    //kktsolver.SetPrecOperator(MLPrec);
    kktsolver.SetPrecOperator(&(*Prec_));

    // specify solver
    kktsolver.SetAztecOption(AZ_solver,AZ_gmres);
    //kktsolver.SetAztecOption(AZ_output,AZ_none);
    kktsolver.SetAztecOption(AZ_output,AZ_all);

    // Set up tolerances.
    if ((mypid==0) && wantstats)
      outfile << "TOL: " << tolerance << endl;

    if (fabs(tolerance) < minthreshold) {

      // Almost zero tolerance.
      kktsolver.SetAztecOption(AZ_conv,AZ_r0);
      kktsolver.SetAztecOption(AZ_kspace, maxit);
      // .. and here we solve
      kktsolver.Iterate(maxit, minprecision);
      if ((mypid==0) && wantstats)
        outfile << kktsolver.NumIters() << endl;

    }

    else if (tolerance < 0) {

      if ((mypid==0) && wantstats) {
        outfile << "\nFIRST CG ITER\n";
        cout << "\nFIRST CG ITER\n";
      }
    
      double normg  = tol[1];
      double delta  = tol[2];
      double smallc = tol[3];
      // First iterative solve with dynamic tolerance.
      // do one solve and then adjust if needed
      kktsolver.SetAztecOption(AZ_conv,AZ_noscaled);
      kktsolver.SetAztecOption(AZ_kspace, innerit);
      // .. and here we solve
      kktsolver.Iterate(innerit, minprecision);
      for (int i=1; i <= numcycles; i++) {
        double normsoln;
        soln.Norm2(&normsoln);
        // dynamic stopping criterion
        double dynmin = 0.0;
        dynmin = min(normsoln/normg, delta/normg);
        dynmin = min(dynmin, smallc);
        dynmin = dynmin*min(normsoln, 1.0);
        if ((mypid==0) && wantstats) {
          outfile << "\nTOL: " << dynmin << endl;
          cout << "\nTOL: " << dynmin << endl;
        }
        if (kktsolver.TrueResidual() > dynmin) {
          if ((mypid==0) && wantstats) {
            outfile << "\nCYCLE\n";
            cout << "\nCYCLE\n";
          }
          kktsolver.Iterate(innerit, minprecision);
        }
      }
      if ((mypid==0) && wantstats)
        outfile << kktsolver.NumIters() << endl;
      
    }
    else {
      
      // Solve to desired preset tolerance, dynamic or fixed.
      if (tol[1] < 0.1)  // if absolute tolerance chosen, for dynamic solves
        kktsolver.SetAztecOption(AZ_conv,AZ_noscaled);
      else               // if relative tolerance chosen, for fixed solves
        kktsolver.SetAztecOption(AZ_conv,AZ_r0);
      kktsolver.SetAztecOption(AZ_kspace, maxit);
      kktsolver.Iterate(maxit, tolerance);
      if ((mypid==0) && wantstats)
        outfile << kktsolver.NumIters() << endl;

    }

  }

  // Pull the solution out into three vectors y, u, and p
    
  for (int i=0; i<rhsy->MyLength(); i++) {
    (*((*y)(0)))[i] = soln[i];
  }
  for (int i=0; i<rhsu->MyLength(); i++) {
    (*((*u)(0)))[i] = soln[i+rhsy->MyLength()];
  }
  for (int i=0; i<rhsp->MyLength(); i++) {
    (*((*p)(0)))[i] = soln[i+rhsy->MyLength()+rhsu->MyLength()];
  }

  outfile.close();
  return 0;

}

Epetra_Comm * GLdistYUEpetraDataPool::getCommPtr()   { return commptr_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLdistYUEpetraDataPool::getA()  { return A_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLdistYUEpetraDataPool::getB()  { return B_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLdistYUEpetraDataPool::getH()  { return H_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLdistYUEpetraDataPool::getR()  { return R_; }

Teuchos::RefCountPtr<Epetra_CrsMatrix> GLdistYUEpetraDataPool::getAugmat()  { return Augmat_; }

Teuchos::RefCountPtr<Epetra_FECrsMatrix> GLdistYUEpetraDataPool::getNpy()  { return Npy_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLdistYUEpetraDataPool::getb()  { return b_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLdistYUEpetraDataPool::getq()  { return q_; }

Teuchos::RefCountPtr<Epetra_FEVector> GLdistYUEpetraDataPool::getNy()  { return Ny_; }

double GLdistYUEpetraDataPool::getbeta()  { return beta_; }

Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> GLdistYUEpetraDataPool::getipcoords()  { return ipcoords_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> GLdistYUEpetraDataPool::getipindx()  { return ipindx_; }

Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> GLdistYUEpetraDataPool::getpcoords()  { return pcoords_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> GLdistYUEpetraDataPool::getpindx()  { return pindx_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> GLdistYUEpetraDataPool::gett()  { return t_; }

Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> GLdistYUEpetraDataPool::gete()  { return e_; }


void GLdistYUEpetraDataPool::computeNy( const Teuchos::RefCountPtr<const Epetra_MultiVector> & y )
{
  Epetra_Map overlapmap(-1, pindx_->M(), (int*)(pindx_)->A(), 1, *commptr_);
  Epetra_Map standardmap(A_->DomainMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> yoverlap = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Epetra_Import Importer(overlapmap, standardmap);
  yoverlap->Import(*y, Importer, Insert);
  nonlinvec(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, yoverlap, Ny_);
}


void GLdistYUEpetraDataPool::computeNpy( const Teuchos::RefCountPtr<const Epetra_MultiVector> & y )
{
  Epetra_Map overlapmap(-1, pindx_->M(), (int*)(pindx_)->A(), 1, *commptr_);
  Epetra_Map standardmap(A_->DomainMap());
  Teuchos::RefCountPtr<Epetra_MultiVector> yoverlap = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Epetra_Import Importer(overlapmap, standardmap);
  yoverlap->Import(*y, Importer, Insert);
  nonlinjac(*commptr_, *ipindx_, *ipcoords_, *pindx_, *pcoords_, *t_, yoverlap, Npy_);
}


void GLdistYUEpetraDataPool::computeAugmat()
{
  Epetra_Map standardmap(A_->DomainMap());
  Epetra_Map bdryctrlmap(B_->DomainMap());

  int indexBase = 1;

  int numstates = standardmap.NumGlobalElements();
  int nummystates = standardmap.NumMyElements();
  int nummycontrols = bdryctrlmap.NumMyElements();

  Epetra_IntSerialDenseVector KKTmapindx(2*nummystates+nummycontrols);
  
  
  // Build KKT map.
  Epetra_IntSerialDenseVector states(nummystates);
  Epetra_IntSerialDenseVector controls(nummycontrols);
  standardmap.MyGlobalElements(states.Values());
  bdryctrlmap.MyGlobalElements(controls.Values());
  for (int i=0; i<nummystates; i++) {
    KKTmapindx(i) = states(i);
    KKTmapindx(nummystates+nummycontrols+i) = 2*numstates+states(i);
  }
  for (int i=0; i<nummycontrols; i++) {
    KKTmapindx(nummystates+i) = numstates+controls(i);
  }
  Epetra_Map KKTmap(-1, 2*nummystates+nummycontrols, KKTmapindx.Values(), indexBase, *(commptr_));
  
  
  // Start building the KKT matrix.
  
  Augmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, KKTmap, 0));

  double one[1];
  one[0]=1.0;
  for (int i=0; i<nummystates+nummycontrols; i++) {
    Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], 1, one, KKTmapindx.Values()+i);
  }
  
  int checkentries=0;
  int nummyentries=0;
  Epetra_SerialDenseVector values(nummyentries);
  Epetra_IntSerialDenseVector indices(nummyentries);
  // Insert A and Npy into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = A_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    A_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
    nummyentries = Npy_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    Npy_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
  }
  // Insert B into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = B_->NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    B_->ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                             indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+2*numstates, nummyentries, values.Values(), 
                                  indices.Values());
  }
  
  bool MakeDataContiguous = false;
  EpetraExt::RowMatrix_Transpose transposer( MakeDataContiguous );
  Epetra_CrsMatrix & transA = dynamic_cast<Epetra_CrsMatrix&>(transposer(*A_));
  Epetra_CrsMatrix & transB = dynamic_cast<Epetra_CrsMatrix&>(transposer(*B_));
  Epetra_CrsMatrix & transNpy = dynamic_cast<Epetra_CrsMatrix&>(transposer(*Npy_));
  // Insert transpose of A and Npy into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = transA.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transA.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], nummyentries, values.Values(), 
                                  indices.Values());
    nummyentries = transNpy.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transNpy.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                  indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    if (nummyentries > 0)
      Augmat_->InsertGlobalValues(KKTmapindx.Values()[i], nummyentries, values.Values(), 
                                  indices.Values());
  }
  // Insert transpose of B into Augmat.
  for (int i=0; i<nummystates; i++) {
    nummyentries = transB.NumMyEntries(i);
    values.Resize(nummyentries);
    indices.Resize(nummyentries);
    transB.ExtractGlobalRowCopy(KKTmapindx.Values()[i], nummyentries, checkentries, values.Values(),
                                indices.Values());
    for (int j=0; j<nummyentries; j++)
      indices[j] = indices[j]+2*numstates;
    int err = 0;
    if (nummyentries > 0)
      err = Augmat_->InsertGlobalValues(KKTmapindx.Values()[i]+numstates, nummyentries,
                                        values.Values(), indices.Values());
    // This will give a nasty message if something goes wrong with the insertion of B transpose.
    if (err < 0) {
      cout << "Insertion of entries failed:\n";
      cout << indices;
      cout << nummyentries << endl;
      cout << "at row: " << KKTmapindx.Values()[i]+numstates << endl << endl;
    }
  }

  Augmat_->FillComplete(KKTmap, KKTmap);
  // End building the KKT matrix.

}



/*int GLdistYUEpetraDataPool::computePrec()
{
  Teuchos::ParameterList MLList;

  // =============================== //
  // Create the ML + IFPACK smoother //
  // =============================== //
                                                                                
  ML_Epetra::SetDefaults("DD",MLList);
  MLList.set("smoother: pre or post", "post");
  MLList.set("PDE equations", 1);

  // now we have to specify which IFPACK preconditioner should be
  // built. Any value that is valid for the IFPACK factory. We also need
  // to define the overlap (>= 0).
                                                                                
  MLList.set("smoother: ifpack type", "Amesos");
  MLList.set("smoother: ifpack overlap", 4);

  // Then, all parameters can will control the definition of the IFPACK
  // smoother are inserted in IFPACKList. In this case, we specify the fill-in
  // factor. For a list of supported parameters, please consult the IFPACK
  // documentation. For example, IFPACK preconditioner "Amesos" or
  // "Amesos stand-alone" can be used to solve with an LU
  // factorization on each domain.
  MLList.sublist("smoother: ifpack list").set("amesos: solver type", "Amesos_Umfpack");

  // we can now build the preconditioner...
                                                                                
  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner( &(*Augmat_), MLList, true);
}*/



int GLdistYUEpetraDataPool::computePrec()
{

  const RCP<FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if (useStratimikos) {
    
    //*out << "\nUsing Stratimikos!\n";

    using Teuchos::RCP; using Teuchos::FancyOStream;
    using Teuchos::describe; using Teuchos::rcpFromRef;

    if (is_null(Augmat_lowsFactory_)) {

      RCP<Teuchos::ParameterList> pl =
        Teuchos::getParametersFromXmlFile(stratimikosXmlFile);

      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      linearSolverBuilder.setParameterList(pl);

      Augmat_lowsFactory_ = linearSolverBuilder.createLinearSolveStrategy("");
    
      Augmat_lowsFactory_->setOStream(out);
      Augmat_lowsFactory_->setVerbLevel(Teuchos::VERB_LOW);

    }

    Augmat_lows_ = Augmat_lowsFactory_->createOp();
    thyra_Augmat_ = Thyra::epetraLinearOp(Augmat_);

    Thyra::initializeOp<double>(*Augmat_lowsFactory_, thyra_Augmat_,
      Augmat_lows_.ptr());
    *out << "\nlows = " << describe(*Augmat_lows_, Teuchos::VERB_MEDIUM);

  }
  else {

    // Create Ifpack preconditioner for AztecOO!

    Teuchos::ParameterList List;

    // allocates an IFPACK factory. No data is associated
    // to this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    string PrecType = "Amesos";
    int OverlapLevel = 4; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.

    //if (Prec_ != 0)
    //  delete Prec_;
    Prec_ = rcp(Factory.Create(PrecType, &(*Augmat_), OverlapLevel));
    assert(Prec_ != null);

    // specify the Amesos solver to be used.
    // If the selected solver is not available,
    // IFPACK will try to use Amesos' KLU (which is usually always
    // compiled). Amesos' serial solvers are:
    // "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
    // List.set("amesos: solver type", "Amesos_Umfpack");
    List.set("amesos: solver type", "Amesos_Klu");

    // sets the parameters
    IFPACK_CHK_ERR(Prec_->SetParameters(List));

    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    // At this call, Amesos will perform the symbolic factorization.
    IFPACK_CHK_ERR(Prec_->Initialize());

    // Builds the preconditioners, by looking for the values of
    // the matrix. At this call, Amesos will perform the
    // numeric factorization.
    IFPACK_CHK_ERR(Prec_->Compute());

  }

  return 0;

}



void GLdistYUEpetraDataPool::PrintSolutionMatlab( const Teuchos::RefCountPtr<const Epetra_Vector> & x )
{
  Vector2MATLAB(*x, cout);
}  


void GLdistYUEpetraDataPool::PrintSolutionVTK( const Teuchos::RefCountPtr<const Epetra_Vector> & y )
{
  Epetra_Map standardmap(A_->DomainMap());
  int numstates = standardmap.NumGlobalElements();
  int IndexBase = 1;

  Teuchos::RefCountPtr<Epetra_Map> printmap;

  if (commptr_->MyPID() == 0)
    printmap = Teuchos::rcp(new Epetra_Map(numstates, numstates, IndexBase, *commptr_));
  else
    printmap = Teuchos::rcp(new Epetra_Map(numstates, 0, IndexBase, *commptr_));

  Epetra_Import importer(*printmap, standardmap);

  Epetra_Vector yprint(*printmap);

  yprint.Import(*y, importer, Insert);

  char FileStates[120];

  if (commptr_->MyPID() == 0) {
    Epetra_SerialDenseMatrix ipcoords;
    Epetra_IntSerialDenseMatrix t;
    int numip = 0, numelems = 0;

    char FileNode[120];
    sprintf(FileNode, "%s.node", geomfile_);
    char FileEle[120];
    sprintf(FileEle, "%s.ele", geomfile_);

    FILE* nodefile = fopen(FileNode, "r");
    if (!nodefile)
      cerr << "cannot open input file " << FileNode << endl;
    fscanf(nodefile, "%d", &numip);
    fscanf(nodefile, "%*[^\n]");   // Skip to the End of the Line 
    fscanf(nodefile, "%*1[\n]");   // Skip One Newline 
    ipcoords.Shape(numip, 2);

    for (int i=0; i<numip; i++) {
      fscanf(nodefile, "%*d %lf %lf", &ipcoords(i,0), &ipcoords(i,1));
      fscanf(nodefile, "%*[^\n]");   // Skip to the End of the Line
      fscanf(nodefile, "%*1[\n]");   // Skip One Newline
    }

    fclose(nodefile);

    FILE* elefile = fopen(FileEle, "r");
    if (!elefile)
      cerr << "cannot open input file " << FileEle << endl;
    fscanf(elefile, "%d", &numelems);
    fscanf(elefile, "%*[^\n]");   // Skip to the End of the Line
    fscanf(elefile, "%*1[\n]");   // Skip One Newline
    t.Shape(numelems, 3);

    for (int i=0; i<numelems; i++) {
      fscanf(elefile, "%*d %d %d %d", &t(i,0), &t(i,1), &t(i,2));
      // Vtk starts at zero.
      t(i,0)--; t(i,1)--; t(i,2)--;
      fscanf(elefile, "%*[^\n]");   // Skip to the End of the Line
      fscanf(elefile, "%*1[\n]");   // Skip One Newline
    }

    fclose(elefile);

    sprintf(FileStates,   "%s%04d.vtk", "./states",   0);

    // Open files for writing.
    FILE* statefile = fopen(FileStates, "w");
    if (!statefile)
      cerr << "cannot open output file " << FileStates << endl;

    // Write headers.
    fprintf(statefile,   "# vtk DataFile Version 2\nHeatApp States\nASCII\nDATASET UNSTRUCTURED_GRID\n\n");

    // Write POINTS headers.
    fprintf(statefile,   "POINTS %d double\n", numip);

    // Write POINTS data, including the z-coordinate.
    for (int Row=0; Row<numip; Row++) {
      fprintf(statefile,   "%16.14e %16.14e %16.14e\n", ipcoords(Row,0), ipcoords(Row,1),
              yprint[Row]);
    }
    fprintf(statefile,   "\n");

    // Write CELL headers.
    fprintf(statefile,   "CELLS %d %d\n", numelems, 4*numelems);

    // Write CELL data.
    for (int Row=0; Row<numelems; Row++) {
      fprintf(statefile,   "%d %d %d %d\n", 3, t(Row,0), t(Row,1), t(Row,2));
    }
    fprintf(statefile,   "\n");

    // Write CELL_TYPES headers.
    fprintf(statefile,   "CELL_TYPES %d\n", numelems);
    for (int Row=0; Row<numelems; Row++) {
      fprintf(statefile,   "%d\n", 5);
    }
    fprintf(statefile,   "\n");

    // Write POINT_DATA headers.
    fprintf(statefile,   "POINT_DATA %d\nSCALARS scalars double 1\nLOOKUP_TABLE default\n", numip);

    // Write state snapshot.
    for(int Row=0; Row < numstates; Row++ )
      fprintf(statefile, "%16.14e\n", yprint[Row]);

    // Close files.
    fclose(statefile);

  }

}


} //namespace GLdistApp
