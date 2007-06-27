
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief <tt>LinearOpWithSolveFactoryBase</tt> subclass implemented in terms
 * of <tt>Belos</tt>.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Belos_Thyra_adapters_grp
 */
template<class Scalar>
class BelosLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<Scalar> {
public:

  /** \name Public types */
  //@{
  /** \brief . */

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType  MagnitudeType;

  //@}

  /** \name Parameter names for Paramter List */
  //@{

  /** \brief . */
  static const std::string  SolverType_name;
  /** \brief . */           
  static const std::string  SolverType_default;
  /** \brief . */           
  static const std::string  MaxIters_name;
  /** \brief . */           
  static const int          MaxIters_default;
  /** \brief . */           
  static const std::string  MaxRestarts_name;
  /** \brief . */           
  static const int          MaxRestarts_default;
  /** \brief . */           
  static const std::string  BlockSize_name;
  /** \brief . */           
  static const int          BlockSize_default;
  /** \brief . */           
  static const std::string  AdjustableBlockSize_name;
  /** \brief . */           
  static const bool         AdjustableBlockSize_default;
  /** \brief . */           
  static const std::string  DefaultRelResNorm_name;
  /** \brief . */
  static const MagnitudeType DefaultRelResNorm_default;
  /** \brief . */
  static const std::string  OrthoType_name;
  /** \brief . */           
  static const std::string  OrthoType_default;
  /** \brief . */           
  static const std::string  Restart_Timers_name;
  /** \brief . */           
  static const bool         Restart_Timers_default;
  /** \brief . */
  static const std::string  GMRES_name;
  /** \brief . */           
  static const std::string  GMRES_MaxNumberOfKrylovVectors_name;
  /** \brief . */           
  static const int          GMRES_MaxNumberOfKrylovVectors_default;
  /** \brief . */           
  static const std::string  GMRES_Variant_name;
  /** \brief . */           
  static const std::string  GMRES_Variant_default;
  /** \brief . */           
  static const std::string  Outputter_name;
  /** \brief . */           
  static const std::string  Outputter_OutputFrequency_name;
  /** \brief . */           
  static const int          Outputter_OutputFrequency_default;
  /** \brief . */           
  static const std::string  Outputter_OutputMaxResOnly_name;
  /** \brief . */           
  static const bool         Outputter_OutputMaxResOnly_default;

  //@}

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct without preconditioner factory. */
  BelosLinearOpWithSolveFactory();

  /** \brief Calls <tt>this->setPreconditionerFactory(precFactory)</tt.  . */
  BelosLinearOpWithSolveFactory(
    const Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  &precFactory
    );

  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{
  /** \brief Returns true . */
  bool acceptsPreconditionerFactory() const;
  /** \brief . */
  void setPreconditionerFactory(
    const Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  &precFactory
    ,const std::string                                              &precFactoryName
    );
  /** \brief . */
  Teuchos::RCP<PreconditionerFactoryBase<Scalar> > getPreconditionerFactory() const;
  /** \brief . */
  void unsetPreconditionerFactory(
    Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  *precFactory
    ,std::string                                              *precFactoryName
    );
  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<Scalar> &fwdOpSrc ) const;
  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > createOp() const;
  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                                   *Op
    ,const ESupportSolveUse                                          supportSolveUse
    ) const;
  /** \brief . */
  void initializeAndReuseOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                                   *Op
    ) const;
  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<Scalar>                               *Op
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >    *fwdOpSrc
    ,Teuchos::RCP<const PreconditionerBase<Scalar> >    *prec
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >    *approxFwdOpSrc
    ,ESupportSolveUse                                           *supportSolveUse
    ) const;
  /** \brief . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;
  /** \brief . */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;
  /** \brief . */
  void initializeApproxPreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;
  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RCP<PreconditionerFactoryBase<Scalar> >  precFactory_;
  std::string                                               precFactoryName_;
  Teuchos::RCP<Teuchos::ParameterList>              thisValidParamList_;
  Teuchos::RCP<Teuchos::ParameterList>              paramList_;
  bool                                                      useGmres_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RCP<const Teuchos::ParameterList> generateAndGetValidParameters();

  void updateThisValidParamList();

  void initializeOpImpl(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec
    ,const bool                                                         reusePrec
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;

};

//@}

} // namespace Thyra

#endif // THYRA_BELOS_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
