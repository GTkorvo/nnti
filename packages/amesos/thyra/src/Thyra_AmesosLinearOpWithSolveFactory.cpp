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

#ifndef __sun

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolve.hpp"
#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Teuchos_dyn_cast.hpp"

#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_PASTIX
#include "Amesos_Pastix.h"
#endif
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_PARDISO
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
#include "Amesos_Paraklete.h"
#endif

namespace {
const std::string epetraFwdOp_str = "epetraFwdOp";
} // namespace

namespace Thyra {


// Parameter names for Paramter List

const std::string AmesosLinearOpWithSolveFactory::SolverType_name = "Solver Type";

const std::string AmesosLinearOpWithSolveFactory::RefactorizationPolicy_name = "Refactorization Policy";

const std::string AmesosLinearOpWithSolveFactory::ThrowOnPreconditionerInput_name = "Throw on Preconditioner Input";

const std::string AmesosLinearOpWithSolveFactory::AMESOS_name = "AMESOS";

// Constructors/initializers/accessors

AmesosLinearOpWithSolveFactory::AmesosLinearOpWithSolveFactory(
  const Amesos::ESolverType                            solverType
  ,const Amesos::ERefactorizationPolicy                refactorizationPolicy
  ,const bool                                          throwOnPrecInput
    )
  :solverType_(solverType)
  ,refactorizationPolicy_(refactorizationPolicy)
  ,throwOnPrecInput_(throwOnPrecInput)
  ,epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
{}

// Overridden from LinearOpWithSolveFactoryBase

bool AmesosLinearOpWithSolveFactory::isCompatible(
  const LinearOpBase<double> &fwdOp
  ) const
{
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    Teuchos::rcp(&fwdOp,false)
    ,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,&epetraFwdOpAdjointSupport,&epetraFwdOpScalar
    );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}

Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
AmesosLinearOpWithSolveFactory::createOp() const
{
  return Teuchos::rcp(new AmesosLinearOpWithSolve());
}

void AmesosLinearOpWithSolveFactory::initializeOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,LinearOpWithSolveBase<double>                             *Op
  ,const ESupportSolveUse                                    supportSolveUse
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  //
  // Unwrap and get the forward Epetra_Operator object
  //
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,&epetraFwdOpAdjointSupport,&epetraFwdOpScalar
    );
  // Get the AmesosLinearOpWithSolve object
  AmesosLinearOpWithSolve
    *amesosOp = &Teuchos::dyn_cast<AmesosLinearOpWithSolve>(*Op);
  //
  // Determine if we must start over or not
  //
  bool startOver = ( amesosOp->get_amesosSolver()==Teuchos::null );
  if(!startOver) {
    startOver =
      (
        epetraFwdOpTransp != amesosOp->get_amesosSolverTransp() ||
        epetraFwdOp.get() != amesosOp->get_epetraLP()->GetOperator()
        // We must start over if the matrix object changes.  This is a
        // weakness of Amesos but there is nothing I can do about this right
        // now!
        );
  }
  //
  // Update the amesos solver
  //
  if(startOver) {
    //
    // This LOWS object has not be initialized yet or is not compatible with the existing
    // 
    // so this is where we setup everything from the ground up.
    //
    // Create the linear problem and set the operator with memory of RCP to Epetra_Operator view!
    Teuchos::RefCountPtr<Epetra_LinearProblem>
      epetraLP = Teuchos::rcp(new Epetra_LinearProblem());
    epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp));
    Teuchos::set_extra_data< Teuchos::RefCountPtr<const Epetra_Operator> >( epetraFwdOp, epetraFwdOp_str, &epetraLP );
    // Create the concrete solver
    Teuchos::RefCountPtr<Amesos_BaseSolver>
      amesosSolver;
    switch(solverType()) {
      case Thyra::Amesos::LAPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Lapack(*epetraLP));
        break;
#ifdef HAVE_AMESOS_KLU
      case Thyra::Amesos::KLU :
        amesosSolver = Teuchos::rcp(new Amesos_Klu(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_PASTIX
      case Thyra::Amesos::PASTIX :
        amesosSolver = Teuchos::rcp(new Amesos_Pastix(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_MUMPS
      case Thyra::Amesos::MUMPS :
        amesosSolver = Teuchos::rcp(new Amesos_Mumps(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_SCALAPACK
      case Thyra::Amesos::SCALAPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Scalapack(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_UMFPACK
      case Thyra::Amesos::UMFPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Umfpack(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
      case Thyra::Amesos::SUPERLUDIST :
        amesosSolver = Teuchos::rcp(new Amesos_Superludist(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_SUPERLU
      case Thyra::Amesos::SUPERLU :
        amesosSolver = Teuchos::rcp(new Amesos_Superlu(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_DSCPACK
      case Thyra::Amesos::DSCPACK :
        amesosSolver = Teuchos::rcp(new Amesos_Dscpack(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_PARDISO
      case Thyra::Amesos::PARDISO :
        amesosSolver = Teuchos::rcp(new Amesos_Pardiso(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_TAUCS
      case Thyra::Amesos::TAUCS :
        amesosSolver = Teuchos::rcp(new Amesos_Taucs(*epetraLP));
        break;
#endif
#ifdef HAVE_AMESOS_PARAKLETE
      case Thyra::Amesos::PARAKLETE :
        amesosSolver = Teuchos::rcp(new Amesos_Paraklete(*epetraLP));
        break;
#endif
      default:
        TEST_FOR_EXCEPTION(
          true, std::logic_error
          ,"Error, the solver type ID = " << solverType() << " is invalid!"
          );
    }
    // Set the parameters
    if(paramList_.get()) amesosSolver->SetParameters(paramList_->sublist("AMESOS"));
    // Do the initial factorization
    amesosSolver->SymbolicFactorization();
    amesosSolver->NumericFactorization();
    // Initialize the LOWS object and we are done!
    amesosOp->initialize(fwdOp,epetraLP,amesosSolver,epetraFwdOpTransp,epetraFwdOpScalar);
  }
  else {
    //
    // This LOWS object has already be initialized once so we must just reset
    // the matrix and refactor it.
    //
    // Get non-const pointers to the linear problem and the amesos solver.
    // These const-casts are just fine since the amesosOp in non-const.
    Teuchos::RefCountPtr<Epetra_LinearProblem>
      epetraLP = Teuchos::rcp_const_cast<Epetra_LinearProblem>(amesosOp->get_epetraLP());
    Teuchos::RefCountPtr<Amesos_BaseSolver>
      amesosSolver = amesosOp->get_amesosSolver();
    // Reset the forward operator with memory of RCP to Epetra_Operator view!
    epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp));
    Teuchos::get_extra_data< Teuchos::RefCountPtr<const Epetra_Operator> >(epetraLP,epetraFwdOp_str) = epetraFwdOp;
    // Reset the parameters
    if(paramList_.get()) amesosSolver->SetParameters(paramList_->sublist(AMESOS_name));
    // Repivot if asked
    if(refactorizationPolicy()==Amesos::REPIVOT_ON_REFACTORIZATION)
      amesosSolver->SymbolicFactorization();
    amesosSolver->NumericFactorization();
    // Reinitialize the LOWS object and we are done! (we must do this to get the
    // possibly new transpose and scaling factors back in)
    amesosOp->initialize(fwdOp,epetraLP,amesosSolver,epetraFwdOpTransp,epetraFwdOpScalar);
  }
}

bool AmesosLinearOpWithSolveFactory::supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const
{
  return false;
}

void AmesosLinearOpWithSolveFactory::initializePreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >             &fwdOp
  ,const Teuchos::RefCountPtr<const PreconditionerBase<double> >      &prec
  ,LinearOpWithSolveBase<double>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(
    this->throwOnPrecInput(), std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' does not support precondtioners "
    "and has been configured to throw this exception when the  initializePreconditionedOp(...) function is called!"
    );
  this->initializeOp(fwdOp,Op,supportSolveUse); // Ignore the precondtioner!
}

void AmesosLinearOpWithSolveFactory::initializePreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >             &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<double> >            &approxFwdOp
  ,LinearOpWithSolveBase<double>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(
    this->throwOnPrecInput(), std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' does not support precondtioners "
    "and has been configured to throw this exception when the  initializePreconditionedOp(...) function is called!"
    );
  this->initializeOp(fwdOp,Op,supportSolveUse); // Ignore the precondtioner!
}

void AmesosLinearOpWithSolveFactory::uninitializeOp(
  LinearOpWithSolveBase<double>                               *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >          *fwdOp
  ,Teuchos::RefCountPtr<const PreconditionerBase<double> >    *prec
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >          *approxFwdOp
  ,ESupportSolveUse                                           *supportSolveUse
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  AmesosLinearOpWithSolve
    *amesosOp = &Teuchos::dyn_cast<AmesosLinearOpWithSolve>(*Op);
  Teuchos::RefCountPtr<const LinearOpBase<double> >
    _fwdOp = amesosOp->extract_fwdOp(); // Will be null if uninitialized!
  if(_fwdOp.get()) {
    // Erase the Epetra_Operator view of the forward operator!
    Teuchos::RefCountPtr<Epetra_LinearProblem> epetraLP = amesosOp->get_epetraLP();
    Teuchos::get_extra_data< Teuchos::RefCountPtr<const Epetra_Operator> >(epetraLP,epetraFwdOp_str) = Teuchos::null;
    // Note, we did not erase the address of the operator in epetraLP->GetOperator() since
    // it seems that the amesos solvers do not recheck the value of GetProblem()->GetOperator()
    // so you had better not rest this!
  }
  if(fwdOp) *fwdOp = _fwdOp; // It is fine if the client does not want this object back!
  if(prec) *prec = Teuchos::null; // We never keep a preconditioner!
  if(approxFwdOp) *approxFwdOp = Teuchos::null; // We never keep a preconditioner!
}

// Overridden from ParameterListAcceptor

void AmesosLinearOpWithSolveFactory::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),0); // Only validate this level for now!
  paramList_ = paramList;
  solverType(
    Amesos::solverTypeNameToEnumMap.get<Amesos::ESolverType>(
      paramList_->get(
        SolverType_name
        ,Amesos::toString(solverType())
        )
      ,paramList_->name()+"->"+SolverType_name
      )
    );
  refactorizationPolicy(
    Amesos::refactorizationPolicyNameToEnumMap.get<Amesos::ERefactorizationPolicy>(
      paramList_->get(
        RefactorizationPolicy_name
        ,Amesos::toString(refactorizationPolicy())
        )
      ,paramList_->name()+"->"+RefactorizationPolicy_name
      )
    );
  throwOnPrecInput(paramList_->get(ThrowOnPreconditionerInput_name,throwOnPrecInput()));
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::getValidParameters() const
{
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

std::string AmesosLinearOpWithSolveFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::AmesosLinearOpWithSolveFactory{";
  oss << "solverType=" << toString(solverType_);
  oss << "}";
  return oss.str();
}

// private

Teuchos::RefCountPtr<const Teuchos::ParameterList>
AmesosLinearOpWithSolveFactory::generateAndGetValidParameters()
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("AmesosLinearOpWithSolveFactory"));
    validParamList->set(
      SolverType_name
#ifdef HAVE_AMESOS_KLU
      ,Amesos::toString(Amesos::KLU)
#else
      ,Amesos::toString(Amesos::LAPACK)
#endif
      );
    validParamList->set(RefactorizationPolicy_name,Amesos::toString(Amesos::REPIVOT_ON_REFACTORIZATION));
    validParamList->set(ThrowOnPreconditionerInput_name,bool(true));
    validParamList->sublist(AMESOS_name);
  }
  return validParamList;
}

} // namespace Thyra

#endif // __sun
