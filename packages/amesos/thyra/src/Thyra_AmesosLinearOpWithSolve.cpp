/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Thyra_AmesosLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_MultiVector.h"

namespace Thyra {

// Constructors/initializers/accessors

AmesosLinearOpWithSolve::AmesosLinearOpWithSolve()
{}

AmesosLinearOpWithSolve::AmesosLinearOpWithSolve(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,const Teuchos::RefCountPtr<Epetra_LinearProblem>          &epetraLP
  ,const Teuchos::RefCountPtr<Amesos_BaseSolver>             &amesosSolver
  ,const ETransp                                             amesosSolverTransp
  ,const double                                              amesosSolverScalar
  )
{
  this->initialize(fwdOp,epetraLP,amesosSolver,amesosSolverTransp,amesosSolverScalar);
}

void AmesosLinearOpWithSolve::initialize(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,const Teuchos::RefCountPtr<Epetra_LinearProblem>          &epetraLP
  ,const Teuchos::RefCountPtr<Amesos_BaseSolver>             &amesosSolver
  ,const ETransp                                             amesosSolverTransp
  ,const double                                              amesosSolverScalar
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEST_FOR_EXCEPT(epetraLP.get()==NULL);
  TEST_FOR_EXCEPT(amesosSolver.get()==NULL);
  TEST_FOR_EXCEPT(epetraLP->GetLHS()!=NULL);
  TEST_FOR_EXCEPT(epetraLP->GetRHS()!=NULL);
#endif
  fwdOp_              = fwdOp;
  epetraLP_           = epetraLP;
  amesosSolver_       = amesosSolver;
  amesosSolverTransp_ = amesosSolverTransp;
  amesosSolverScalar_ = amesosSolverScalar;
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AmesosLinearOpWithSolve::extract_fwdOp()
{
  Teuchos::RefCountPtr<const LinearOpBase<double> > _fwdOp = fwdOp_;
  fwdOp_ = Teuchos::null;
  return _fwdOp;
}

void AmesosLinearOpWithSolve::reset_fwdOp( const Teuchos::RefCountPtr<const LinearOpBase<double> > &fwdOp )
{
  fwdOp_ = fwdOp;
}

void AmesosLinearOpWithSolve::uninitialize(
  Teuchos::RefCountPtr<const LinearOpBase<double> >    *fwdOp
  ,Teuchos::RefCountPtr<Epetra_LinearProblem>          *epetraLP
  ,Teuchos::RefCountPtr<Amesos_BaseSolver>             *amesosSolver
  ,ETransp                                             *amesosSolverTransp
  ,double                                              *amesosSolverScalar
  )
{

  if(fwdOp)              *fwdOp              = fwdOp_;
  if(epetraLP)           *epetraLP           = epetraLP_;
  if(amesosSolver)       *amesosSolver       = amesosSolver_;
  if(amesosSolverTransp) *amesosSolverTransp = amesosSolverTransp_;
  if(amesosSolverScalar) *amesosSolverScalar = amesosSolverScalar_;

  fwdOp_              = Teuchos::null;
  epetraLP_           = Teuchos::null;
  amesosSolver_       = Teuchos::null;
  amesosSolverTransp_ = NOTRANS;
  amesosSolverScalar_ = 0.0;

}

// Overridden from LinearOpBase

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AmesosLinearOpWithSolve::range() const
{
  return ( fwdOp_.get() ? fwdOp_->range() : Teuchos::null );
}

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AmesosLinearOpWithSolve::domain() const
{
  return  ( fwdOp_.get() ? fwdOp_->domain() : Teuchos::null );
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AmesosLinearOpWithSolve::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// Overridden from Teuchos::Describable

std::string AmesosLinearOpWithSolve::description() const
{
  std::ostringstream oss;
  oss << "Thyra::AmesosLinearOpWithSolve";
  if(amesosSolver_.get()) {
    oss << "(fwdOp=\'"<<fwdOp_->description()<<"\'"
        << ",amesosSolver=\'"<<typeid(*amesosSolver_).name()<<"\')";
  }
  return oss.str();
}

// protected

// Overridden from SingleScalarLinearOpBase

bool AmesosLinearOpWithSolve::opSupported(ETransp M_trans) const
{
  return ::Thyra::opSupported(*fwdOp_,M_trans);
}

void AmesosLinearOpWithSolve::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<double>    &X
  ,MultiVectorBase<double>          *Y
  ,const double                     alpha
  ,const double                     beta
  ) const
{
  Thyra::apply( *fwdOp_, M_trans, X, Y, alpha, beta );
}

// Overridden from SingleScalarLinearOpWithSolveBase

bool AmesosLinearOpWithSolve::solveSupportsTrans(ETransp M_trans) const
{
  return true; // ToDo: Determine if the solver supports adjoints or not!
}

bool AmesosLinearOpWithSolve::solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const
{
  return true; // I am a direct solver so I should be able to do it all!
}

// Overridden from SingleRhsLinearOpWithSolveBase

void AmesosLinearOpWithSolve::solve(
  const ETransp                              M_trans
  ,const MultiVectorBase<double>             &B
  ,MultiVectorBase<double>                   *X
  ,const int                                 numBlocks
  ,const BlockSolveCriteria<double>          blockSolveCriteria[]
  ,SolveStatus<double>                       blockSolveStatus[]
  ) const
{
  typedef SolveCriteria<double>  SC;
  typedef SolveStatus<double>    SS;
#ifdef _DEBUG
  TEST_FOR_EXCEPT(X==NULL);
  TEST_FOR_EXCEPT(blockSolveCriteria==NULL && blockSolveStatus!=NULL);
#endif
  //
  // Get the op(...) range and domain maps
  //
  const ETransp amesosOpTransp = real_trans(trans_trans(amesosSolverTransp_,M_trans));
  const Epetra_Operator *amesosOp = epetraLP_->GetOperator();
  const Epetra_Map
    &opRangeMap  = ( amesosOpTransp == NOTRANS ? amesosOp->OperatorRangeMap()  : amesosOp->OperatorDomainMap() ),
    &opDomainMap = ( amesosOpTransp == NOTRANS ? amesosOp->OperatorDomainMap() : amesosOp->OperatorRangeMap()  );
  //
  // Get Epetra_MultiVector views of B and X
  //
  Teuchos::RefCountPtr<const Epetra_MultiVector>
    epetra_B = get_Epetra_MultiVector(opRangeMap,Teuchos::rcp(&B,false));
  Teuchos::RefCountPtr<Epetra_MultiVector>
    epetra_X = get_Epetra_MultiVector(opDomainMap,Teuchos::rcp(X,false));
  //
  // Set B and X in the linear problem
  //
  epetraLP_->SetLHS(&*epetra_X);
  epetraLP_->SetRHS(const_cast<Epetra_MultiVector*>(&*epetra_B)); // Should be okay but cross your fingers!
  //
  // Solve the linear system
  //
  const bool oldUseTranspose = amesosSolver_->UseTranspose();
  amesosSolver_->SetUseTranspose(amesosOpTransp==TRANS);
  TEST_FOR_EXCEPTION(
    0!=amesosSolver_->Solve(), CatastrophicSolveFailure
    ,"Error, the Amesos solver of type \'"<<typeid(*amesosSolver_).name()<<"\' could not perform the solve!"
    );
  amesosSolver_->SetUseTranspose(oldUseTranspose);
  //
  // Unset B and X
  //
  epetraLP_->SetLHS(NULL);
  epetraLP_->SetRHS(NULL);
  epetra_X = Teuchos::null;
  epetra_B = Teuchos::null;
  //
  // Scale X if needed
  //
  if(amesosSolverScalar_!=1.0)
    Thyra::scale(1.0/amesosSolverScalar_,X);
  //
  // Set the solve status if requested
  //
  if(numBlocks && blockSolveStatus) {
    for( int i = 0; i < numBlocks; ++i ) {
      blockSolveStatus[i].solveStatus
        = (blockSolveCriteria[i].solveCriteria.solveTolType!=SOLVE_TOL_DEFAULT
           ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNKNOWN );
      blockSolveStatus[i].achievedTol = SS::unknownTolerance();
      blockSolveStatus[i].iterations = 1;
    }
  }
}

}	// end namespace Thyra
