/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_solver.H"
#include "mrtr_utils.H"
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/05|
 *----------------------------------------------------------------------*/
MOERTEL::Solver::Solver(Epetra_Comm& comm, int outlevel) :
outlevel_(outlevel),
comm_(comm),
matrix_(null),
matrixisnew_(true),
x_(null),
b_(null),
linearproblem_(null),
amesossolver_(null),
mlprec_(null),
aztecsolver_(null)
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 12/05|
 *----------------------------------------------------------------------*/
MOERTEL::Solver::~Solver()
{
}

/*----------------------------------------------------------------------*
 |  set a linear system (public)                             mwgee 12/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Solver::SetSystem(RefCountPtr<Epetra_CrsMatrix> matrix,
                                RefCountPtr<Epetra_Vector> x,
                                RefCountPtr<Epetra_Vector> b)
{
  matrix_ = matrix;
  x_      = x;
  b_      = b;
  return;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve(RefCountPtr<Teuchos::ParameterList> params,
                            RefCountPtr<Epetra_CrsMatrix> matrix,
                            RefCountPtr<Epetra_Vector> x,
                            RefCountPtr<Epetra_Vector> b,
                            Epetra_CrsMatrix& I,
                            Epetra_CrsMatrix& BWT,
                            Epetra_CrsMatrix& B,
                            Epetra_CrsMatrix& WT)
{
  SetParameters(params.get());
  SetSystem(matrix,x,b);
  I_ = &I;
  BWT_ = &BWT;
  B_ = &B;
  WT_ = &WT;
  return Solve();
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve()
{
  bool ok = false;
  
  //---------------------------------------------------------------------------
  // check the linear system  
  if (x_==null || b_==null || matrix_==null)
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** matrix and/or rhs and/or solution vector are Teuchos::null\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //---------------------------------------------------------------------------
  // check for parameters
  if (params_==NULL)
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** solver parameters are Teuchos::null\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }  

  //---------------------------------------------------------------------------
  // (re)create a linear problem
  if (linearproblem_==null)
    linearproblem_ = rcp(new Epetra_LinearProblem());

  linearproblem_->SetLHS(x_.get());
  linearproblem_->SetRHS(b_.get());
  if (matrixisnew_)
    linearproblem_->SetOperator(matrix_.get());
    
  linearproblem_->CheckInput();
  //---------------------------------------------------------------------------
  // get type of solver to be used
  string solver = params_->get("Solver","None");
  
  //---------------------------------------------------------------------------
  // time the solution process
  Epetra_Time time(Comm());
  time.ResetStartTime();

  //---------------------------------------------------------------------------
  // use Amesos
  if (solver=="Amesos" || solver=="amesos" || solver=="AMESOS")
  {
    ParameterList& amesosparams = params_->sublist("Amesos");
    ok = Solve_Amesos(amesosparams);
    if (!ok)
    {
      cout << "***WRN*** MOERTEL::Solver::Solve:\n"
           << "***WRN*** Solve_Amesos returned false\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
  }

  //---------------------------------------------------------------------------
  // use ML/Aztec
  else if (solver=="ML/Aztec" || solver=="ml/aztec" || solver=="ML" ||
           solver=="ml"       || solver=="Ml"       || solver=="Aztec" ||
           solver=="AZTEC"    || solver=="aztec")
  {
    // see whether we have a spd system
    string system = params_->get("System","None");
    if (system!="SPDSystem"  && system!="spdsystem" && system!="spd_system" && 
        system!="SPD_System" && system!="SPDSYSTEM" && system!="SPD_SYSTEM")
    {
      cout << "***ERR*** MOERTEL::Solver::Solve:\n"
           << "***ERR*** To use ML?Aztec for solution, parameter \"System\" hast to be \"SPDSystem\" \n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    ParameterList& mlparams    = params_->sublist("ML");
    ParameterList& aztecparams = params_->sublist("Aztec");
    ok = Solve_MLAztec(mlparams,aztecparams);
    if (!ok)
    {
      cout << "***WRN*** MOERTEL::Solver::Solve:\n"
           << "***WRN*** Solve_MLAztec returned false\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
  }


  //---------------------------------------------------------------------------
  // unknown solver
  else
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** solver type is unknown\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    ok = false;
  }


  //---------------------------------------------------------------------------
  // time the solution process
  double t = time.ElapsedTime();
  if (OutLevel()>5 && Comm().MyPID()==0)
    cout << "MOERTEL (Proc 0): Solution of system of equations in " << t << " sec\n";

  return ok;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (private)                          mwgee 12/05|
 |  using Amesos                                                        |
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve_Amesos(ParameterList& amesosparams)
{
  int ok = 0;

  //---------------------------------------------------------------------------
  // which amesos solver
  string solver       = amesosparams.get("Solver","None");
  bool   usetranspose = amesosparams.get("UseTranspose",false);
  if (solver=="None")
  {
    cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
         << "***ERR*** No Amesos solver chosen\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //---------------------------------------------------------------------------
  // new amesos solver
  if (matrixisnew_ || amesossolver_==null)  
  {
    amesossolver_ = null;
    Amesos Factory;
    if (!Factory.Query(solver))
    {
      cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
           << "***ERR*** Amesos solver '" << solver << "' not supported\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    amesossolver_ = rcp(Factory.Create(solver,*linearproblem_));
    if (amesossolver_.get()==0)
    {
      cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
           << "***ERR*** Could not create Amesos solver '" << solver << "'\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    amesossolver_->SetUseTranspose(usetranspose);
    ok += amesossolver_->SetParameters(amesosparams);
    ok += amesossolver_->SymbolicFactorization();
    ok += amesossolver_->NumericFactorization();
    matrixisnew_ = false;
    ok += amesossolver_->Solve();
  }
  // neither the matrix is new nor the solver
  else
  {
    ok = amesossolver_->Solve();
  }
  
  if (ok==0) return true;
  else
  {
    cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
         << "***ERR*** Amesos returned " << ok << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  return false;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (private)                          mwgee 01/06|
 |  using ML and AztecOO                                                |
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve_MLAztec(ParameterList& mlparams, 
                                    ParameterList& aztecparams)
{
  
  // create ML preconditioner if aztec parameter indicates user preconditioner
  string preconditioner = aztecparams.get("AZ_precond","none");
  if (preconditioner=="none")
  {
    cout << "***ERR*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "***ERR*** Aztec parameter \"AZ_precond\" is not set\n"
         << "***ERR*** set to \"AZ_user_precond\" to use ML or to some Aztec internal method\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (preconditioner=="AZ_user_precond")
    if (mlprec_==null || matrixisnew_);
    mlprec_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*matrix_,mlparams),true);
  
#if 1
  // serial and on 1 level only
  const ML* ml = mlprec_->GetML();
  int nlevel = ml->ML_num_actual_levels; 
  Epetra_CrsMatrix* P;
  int maxnnz=0;
  double cputime;
  ML_Operator2EpetraCrsMatrix(&(ml->Pmat[1]),P,maxnnz,false,cputime);
  //cout << *P;
  
  // form coarse I
  Epetra_CrsMatrix* RI       = MOERTEL::MatMatMult(*P,true,*I_,false,OutLevel());  
  Epetra_CrsMatrix* tmp      = MOERTEL::MatMatMult(*RI,false,*P,false,OutLevel());  
  delete RI; RI = NULL;
  Epetra_CrsMatrix* Icoarse  = MOERTEL::StripZeros(*tmp,1.0e-12);
  delete tmp; tmp = NULL;
  //cout << *Icoarse; // correct

  // form BWTcoarse
  // padd WT to be of full size and square
  Epetra_CrsMatrix* WTsquare = new Epetra_CrsMatrix(Copy,I_->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*WT_,false,1.0,*WTsquare,0.0);
  WTsquare->FillComplete(I_->OperatorDomainMap(),I_->OperatorRangeMap());
  // form WTcoarse = R WTsquare P
  Epetra_CrsMatrix* RWTsquare     = MOERTEL::MatMatMult(*P,true,*WTsquare,false,OutLevel());
  Epetra_CrsMatrix* WTcoarse    = MOERTEL::MatMatMult(*RWTsquare,false,*P,false,OutLevel());
  delete RWTsquare;
  delete WTsquare;
  // padd B to be fullsize and square
  Epetra_CrsMatrix* Bsquare = new Epetra_CrsMatrix(Copy,I_->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*B_,false,1.0,*Bsquare,0.0);
  Bsquare->FillComplete(I_->OperatorDomainMap(),I_->OperatorRangeMap());
  // form Bcoarse = R B P
  Epetra_CrsMatrix* RB      = MOERTEL::MatMatMult(*P,true,*Bsquare,false,OutLevel());
  Epetra_CrsMatrix* Bcoarse = MOERTEL::MatMatMult(*RB,false,*P,false,OutLevel());
  delete Bsquare;
  delete RB;
  // form BWTcoarse2 = WTcoarse * Bcoarse
  Epetra_CrsMatrix* BWTcoarse2 = MOERTEL::MatMatMult(*WTcoarse,false,*Bcoarse,false,OutLevel());
  //cout << *BWTcoarse2;


  // form fine grid ImBWT
  Epetra_CrsMatrix* ImBWT = new Epetra_CrsMatrix(Copy,I_->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*I_,false,1.0,*ImBWT,0.0);
  MOERTEL::MatrixMatrixAdd(*BWT_,false,-1.0,*ImBWT,1.0);
  ImBWT->FillComplete();
  //cout << *ImBWT;
  
  // for coarse grid ImBWTcoarse
  Epetra_CrsMatrix* ImBWTcoarse = new Epetra_CrsMatrix(Copy,Icoarse->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*Icoarse,false,1.0,*ImBWTcoarse,0.0);
  MOERTEL::MatrixMatrixAdd(*BWTcoarse2,false,-1.0,*ImBWTcoarse,1.0);
  ImBWTcoarse->FillComplete();
  //cout << *ImBWTcoarse;
  
  // form modified prolongator Pmod
  Epetra_CrsMatrix* tmp1     = MOERTEL::MatMatMult(*ImBWTcoarse,false,*P,true,OutLevel());  
  Epetra_CrsMatrix* tmp2     = MOERTEL::MatMatMult(*tmp1,false,*ImBWT,false,OutLevel());  
  delete tmp1; tmp1 = NULL;
  Epetra_CrsMatrix* tmp3     = MOERTEL::MatMatMult(*BWTcoarse2,false,*P,true,OutLevel());  
  Epetra_CrsMatrix* tmp4     = MOERTEL::MatMatMult(*tmp3,false,*BWT_,false,OutLevel());  
  delete tmp3; tmp3 = NULL;
  Epetra_CrsMatrix* Rmod = new Epetra_CrsMatrix(Copy,P->OperatorDomainMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*tmp2,false,1.0,*Rmod,0.0);
  delete tmp2; tmp2 = NULL;
  MOERTEL::MatrixMatrixAdd(*tmp4,false,1.0,*Rmod,1.0);
  delete tmp4; tmp4 = NULL;
  Rmod->FillComplete(P->OperatorRangeMap(),P->OperatorDomainMap());
  Epetra_CrsMatrix* tmp5 = new Epetra_CrsMatrix(Copy,P->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*Rmod,true,1.0,*tmp5,0.0);
  delete Rmod; Rmod = NULL;
  tmp5->FillComplete(P->OperatorDomainMap(),P->OperatorRangeMap());
  Epetra_CrsMatrix* Pmod = MOERTEL::StripZeros(*tmp5,1.0e-12);
  //cout << *Pmod;
  
  // fine grid matrix is matrix_, fine grid lhs, rhs are b_, x_

  //Constraints are satisfied on the fine grid by the modified system
  //   A x = b
  //     r = b - A x (initial guess is zero here, though and therefore r = b)
  //   A e = r  
  // to satisfy constraints, we have to have WT r = 0
  // then we also have 
  //   B^T e = 0
  
  // make r
  Epetra_Vector* r = new Epetra_Vector(b_->Map());
  matrix_->Multiply(false,*x_,*r);
  r->Update(1.0,*b_,-1.0);
  //cout << *r;

  // evaluate WT r = 0, ok this is true
  //Epetra_Vector* WTr = new Epetra_Vector(b_->Map(),true);
  //WT_->Multiply(false,*r,*WTr);
  //cout << *WT_;
  //cout << *WTr;

  // evaluate RP = I and PR != I
  //Epetra_CrsMatrix* tmp6 = MOERTEL::MatMatMult(*P,true,*P,false,OutLevel());  
  //Epetra_CrsMatrix* RP   = MOERTEL::StripZeros(*tmp6,1.0e-12); delete tmp6;
  //cout << *RP;
  //Epetra_CrsMatrix* tmp7 = MOERTEL::MatMatMult(*P,false,*P,true,OutLevel());
  //Epetra_CrsMatrix* PR   = MOERTEL::StripZeros(*tmp7,1.0e-12); delete tmp7;
  //cout << *PR;
  
  // do WTcoarse = R WT P
  // padd WT to be of full size and square
  /* previously done
  Epetra_CrsMatrix* tmp8 = new Epetra_CrsMatrix(Copy,I_->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*WT_,false,1.0,*tmp8,0.0);
  WT_ = tmp8;
  WT_->FillComplete(I_->RowMap(),I_->RowMap());
  Epetra_CrsMatrix* tmp9    = MOERTEL::MatMatMult(*P,true,*WT_,false,OutLevel());
  Epetra_CrsMatrix* WTcoarse = MOERTEL::MatMatMult(*tmp9,false,*P,false,OutLevel());  
  delete tmp9;*/
  //cout << *WTcoarse;
  
  // restrict rcoarse = R r
  Epetra_Vector* rcoarse = new Epetra_Vector(Icoarse->RowMap(),true);
  P->Multiply(true,*r,*rcoarse);
  //cout << *rcoarse;
  
  // see that WTcoarse rcoarse != 0 ( with standard prolongator )
  Epetra_Vector* WTcoarse_rcoarse = new Epetra_Vector(Icoarse->RowMap(),true);
  WTcoarse->Multiply(false,*rcoarse,*WTcoarse_rcoarse);
  //cout << *WTcoarse_rcoarse;

  // see that WTcoarse2 rcoarse = 0 ( with modified prolongator )
  Pmod->Multiply(true,*r,*rcoarse);
  //cout << *rcoarse;
  WTcoarse->Multiply(false,*rcoarse,*WTcoarse_rcoarse);
  //cout << *WTcoarse_rcoarse;

#endif



#if 1
  // create the Aztec solver
  aztecsolver_ = rcp(new AztecOO());  
  aztecsolver_->SetAztecDefaults();
  aztecsolver_->SetProblem(*linearproblem_);
  aztecsolver_->SetParameters(aztecparams,true);
  if (mlprec_ != null)
    aztecsolver_->SetPrecOperator(mlprec_.get());
  
  // solve it
  double tol  = aztecparams.get("AZ_tol",1.0e-05);
  int maxiter = aztecparams.get("AZ_max_iter",1000);
  aztecsolver_->Iterate(maxiter,tol);
  matrixisnew_ = false;
  const double* azstatus = aztecsolver_->GetAztecStatus();
/*
  if (azstatus[AZ_why] == AZ_normal)
    return true;
  else if (azstatus[AZ_why] == AZ_breakdown)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_breakdown\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (azstatus[AZ_why] == AZ_loss)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_loss\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (azstatus[AZ_why] == AZ_ill_cond)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_ill_cond\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (azstatus[AZ_why] == AZ_maxits)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_maxits\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned unknown status: " << azstatus[AZ_why] << "\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
*/
#endif  


#if 1

  // test whether the solution e (or x as initial guess was zero) satisfies
  // B^T e = 0 ok this is true
  //Epetra_Vector* BTe = new Epetra_Vector(x_->Map(),true);
  //B_->Multiply(true,*x_,*BTe);
  //cout << *BTe;

#endif





  return true;
}





















