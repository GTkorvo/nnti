/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#include "Thyra_MLPreconditionerFactory.hpp"

#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_implicit_cast.hpp"

namespace {


enum EMLProblemType {
  ML_PROBTYPE_NONE,
  ML_PROBTYPE_SMOOTHED_AGGREGATION, 
  ML_PROBTYPE_DOMAIN_DECOMPOSITION,
  ML_PROBTYPE_DOMAIN_DECOMPOSITION_ML,
  ML_PROBTYPE_MAXWELL
};
const Teuchos::Array<std::string> BaseMethodDefaults_valueNames
= Teuchos::tuple<std::string>(
  "none",
  "SA", 
  "DD",
  "DD-ML",
  "maxwell"
  );

const std::string BaseMethodDefaults_name = "Base Method Defaults";
const std::string BaseMethodDefaults_default = "DD";
Teuchos::RefCountPtr<
  Teuchos::StringToIntegralParameterEntryValidator<EMLProblemType>
  >
BaseMethodDefaults_validator;
  
const std::string MLSettings_name = "ML Settings";


} // namespace


namespace Thyra {


using Teuchos::RefCountPtr;
using Teuchos::ParameterList;


// Constructors/initializers/accessors

  
MLPreconditionerFactory::MLPreconditionerFactory()
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
{}


// Overridden from PreconditionerFactoryBase


bool MLPreconditionerFactory::isCompatible(
  const LinearOpSourceBase<double> &fwdOpSrc
  ) const
{
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp epetraFwdOpTransp;
  EApplyEpetraOpAs epetraFwdOpApplyAs;
  EAdjointEpetraOp epetraFwdOpAdjointSupport;
  double epetraFwdOpScalar;
  Teuchos::RefCountPtr<const LinearOpBase<double> >
    fwdOp = fwdOpSrc.getOp();
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,
    &epetraFwdOp,&epetraFwdOpTransp,
    &epetraFwdOpApplyAs,
    &epetraFwdOpAdjointSupport,
    &epetraFwdOpScalar
    );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}


bool MLPreconditionerFactory::applySupportsConj(EConj conj) const
{
  return true;
}


bool MLPreconditionerFactory::applyTransposeSupportsConj(EConj conj) const
{
  return false; // See comment below
}


RefCountPtr<PreconditionerBase<double> >
MLPreconditionerFactory::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<double>());
}


void MLPreconditionerFactory::initializePrec(
  const Teuchos::RefCountPtr<const LinearOpSourceBase<double> > &fwdOpSrc,
  PreconditionerBase<double> *prec,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::set_extra_data;
  using Teuchos::get_optional_extra_data;
  using Teuchos::implicit_cast;
  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);
  const RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::MLPreconditionerFactory::initializePrec(...) ...\n";

  Teuchos::RefCountPtr<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();
#ifdef _DEBUG
  TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEST_FOR_EXCEPT(prec==NULL);
#endif
  //
  // Unwrap and get the forward Epetra_Operator object
  //
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp epetraFwdOpTransp;
  EApplyEpetraOpAs epetraFwdOpApplyAs;
  EAdjointEpetraOp epetraFwdOpAdjointSupport;
  double epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,
    &epetraFwdOpAdjointSupport,&epetraFwdOpScalar
                                             );
  // Validate what we get is what we need
  RefCountPtr<const Epetra_RowMatrix>
    epetraFwdRowMat = rcp_dynamic_cast<const Epetra_RowMatrix>(epetraFwdOp,true);
  TEST_FOR_EXCEPTION(
    epetraFwdOpApplyAs != EPETRA_OP_APPLY_APPLY, std::logic_error
    ,"Error, incorrect apply mode for an Epetra_RowMatrix"
    );
  //
  // Get the concrete precondtioner object
  //
  DefaultPreconditioner<double>
    *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);
  //
  // Get the EpetraLinearOp object that is used to implement the preconditoner linear op
  //
  RefCountPtr<EpetraLinearOp>
    epetra_precOp = rcp_dynamic_cast<EpetraLinearOp>(defaultPrec->getNonconstUnspecifiedPrecOp(),true);
  //
  // Get the embedded ML_Epetra::MultiLevelPreconditioner object if it exists
  //
  Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner> ml_precOp;
  if(epetra_precOp.get())
    ml_precOp = rcp_dynamic_cast<ML_Epetra::MultiLevelPreconditioner>(epetra_precOp->epetra_op(),true);
  //
  // Get the attached forward operator if it exists and make sure that it matches
  //
  if(ml_precOp.get()) {
    // ToDo: Get the forward operator and make sure that it matches what is
    // already being used!
  }
  //
  // Permform initialization if needed
  //
  const bool startingOver = (ml_precOp.get() == NULL);
  if(startingOver) 
  {
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      *out << "\nCreating the initial ML_Epetra::MultiLevelPreconditioner object...\n";
    timer.start(true);
    // Create the initial preconditioner
    ml_precOp = rcp(
      new ML_Epetra::MultiLevelPreconditioner(
        *epetraFwdRowMat, paramList_->sublist(MLSettings_name)
        )
      );
    
    timer.stop();
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      OSTab(out).o() <<"\n=> Creation time = "<<timer.totalElapsedTime()<<" sec\n";
    // RAB: Above, I am just passing a string to ML::Create(...) in order
    // get this code written.  However, in the future, it would be good to
    // copy the contents of what is in ML::Create(...) into a local
    // function and then use switch(...) to create the initial
    // ML_Epetra::MultiLevelPreconditioner object.  This would result in better validation
    // and faster code.
    // Set parameters if the list exists
    if(paramList_.get())
      TEST_FOR_EXCEPT(
        0!=ml_precOp->SetParameterList(paramList_->sublist(MLSettings_name))
        );
    // Initailize the structure for the preconditioner
    //        TEST_FOR_EXCEPT(0!=ml_precOp->Initialize());
  }
  //
  // Attach the epetraFwdOp to the ml_precOp to guarantee that it will not go away
  //
  set_extra_data(epetraFwdOp,"IFPF::epetraFwdOp",&ml_precOp,Teuchos::POST_DESTROY,false);
  //
  // Update the factorization
  //
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the factorization of the preconditioner ...\n";
  timer.start(true);
  TEST_FOR_EXCEPT(0!=ml_precOp->ComputePreconditioner());
  timer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() <<"\n=> Factorization time = "<<timer.totalElapsedTime()<<" sec\n";
  //
  // Compute the conditioner number estimate if asked
  //

  // ToDo: Implement

  //
  // Attach fwdOp to the ml_precOp
  //
  set_extra_data(fwdOp,"IFPF::fwdOp",&ml_precOp,Teuchos::POST_DESTROY,false);
  //
  // Initialize the output EpetraLinearOp
  //
  if(startingOver) {
    epetra_precOp = rcp(new EpetraLinearOp);
  }
  epetra_precOp->initialize(
    ml_precOp
    ,epetraFwdOpTransp
    ,EPETRA_OP_APPLY_APPLY_INVERSE
    ,EPETRA_OP_ADJOINT_UNSUPPORTED  // ToDo: Look into adjoints again.
    );
  //
  // Initialize the preconditioner
  //
  defaultPrec->initializeUnspecified(
    Teuchos::rcp_implicit_cast<LinearOpBase<double> >(epetra_precOp)
    );
  totalTimer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal time = "<<totalTimer.totalElapsedTime()<<" sec\n"
      << "\nLeaving Thyra::MLPreconditionerFactory::initializePrec(...) ...\n";
}


void MLPreconditionerFactory::uninitializePrec(
  PreconditionerBase<double> *prec,
  Teuchos::RefCountPtr<const LinearOpSourceBase<double> > *fwdOp,
  ESupportSolveUse *supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(true);
}


// Overridden from ParameterListAcceptor


void MLPreconditionerFactory::setParameterList(
  Teuchos::RefCountPtr<ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),0);
  paramList_ = paramList;
  const EMLProblemType
    defaultType = BaseMethodDefaults_validator->getIntegralValue(
      *paramList_,BaseMethodDefaults_name,BaseMethodDefaults_default
      );
  if( ML_PROBTYPE_NONE != defaultType ) {
    const std::string
      defaultTypeStr = BaseMethodDefaults_valueNames[defaultType];
    Teuchos::ParameterList defaultParams;
    TEST_FOR_EXCEPTION(
      0!=ML_Epetra::SetDefaults(defaultTypeStr,defaultParams)
      ,Teuchos::Exceptions::InvalidParameterValue
      ,"Error, the ML problem type \"" << defaultTypeStr << "\' is not recongnised by ML!"
      );
    // Note, the only way the above exception message could be generated is if
    // a default problem type was removed from ML_Epetra::SetDefaults(...).
    // When a new problem type is added to this function, it must be added to
    // our enum EMLProblemType along with associated objects ...  In other
    // words, this adapter must be maintained as ML is maintained.  An
    // alternative design would be to just pass in whatever string to this
    // function.  This would improve maintainability but it would not generate
    // very good error messages when a bad string was passed in.  Currenly,
    // the error message attached to the exception will contain the list of
    // valid problem types.
    paramList_->sublist(MLSettings_name).setParametersNotAlreadySet(
      defaultParams);
  }
#ifdef TEUCHOS_DEBUG
  paramList->validateParameters(*this->getValidParameters(),0);
#endif
}


RefCountPtr<ParameterList>
MLPreconditionerFactory::getParameterList()
{
  return paramList_;
}


RefCountPtr<ParameterList>
MLPreconditionerFactory::unsetParameterList()
{
  Teuchos::RefCountPtr<ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


RefCountPtr<const ParameterList>
MLPreconditionerFactory::getParameterList() const
{
  return paramList_;
}


RefCountPtr<const ParameterList>
MLPreconditionerFactory::getValidParameters() const
{
  using Teuchos::tuple;
  static RefCountPtr<const ParameterList> validPL;
  if(is_null(validPL)) {
    RefCountPtr<ParameterList>
      pl = rcp(new ParameterList());
    BaseMethodDefaults_validator = rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<EMLProblemType>(
        tuple<std::string>(
        "none",
        "SA", 
        "DD",
        "DD-ML",
        "maxwell"
        ),
      tuple<std::string>(
        "Do not set any default parameters",
        "Set default parameters for a smoothed aggregation method",
        "Set default parameters for a domain decomposition method",
        "Set default parameters for a domain decomposition method special to ML",
        "Set default parameters for a Maxwell-type of linear operator"
        ),
      tuple<EMLProblemType>(
        ML_PROBTYPE_NONE,
        ML_PROBTYPE_SMOOTHED_AGGREGATION, 
        ML_PROBTYPE_DOMAIN_DECOMPOSITION,
        ML_PROBTYPE_DOMAIN_DECOMPOSITION_ML,
        ML_PROBTYPE_MAXWELL
        ),
      BaseMethodDefaults_name
      )
    );
    pl->set(BaseMethodDefaults_name,BaseMethodDefaults_default,
      "Select the default method type which also sets parameter defaults\n"
      "in the sublist \"" + MLSettings_name + "\"!",
      BaseMethodDefaults_validator
      );
    pl->sublist(MLSettings_name);
    validPL = pl;
  }
  return validPL;
}


// Public functions overridden from Teuchos::Describable


std::string MLPreconditionerFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::MLPreconditionerFactory";
  return oss.str();
}


} // namespace Thyra
