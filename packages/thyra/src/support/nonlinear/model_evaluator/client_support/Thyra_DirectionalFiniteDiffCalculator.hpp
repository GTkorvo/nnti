
#ifndef THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_HPP
#define THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Thyra {

namespace DirectionalFiniteDiffCalculatorTypes {

/** \brief .
 * \relates DirectionalFiniteDiffCalculator
 */
enum EFDMethodType {
  FD_ORDER_ONE           ///< Use O(eps) one sided finite differences (cramped bounds)
  ,FD_ORDER_TWO          ///< Use O(eps^2) one sided finite differences (cramped bounds)
  ,FD_ORDER_TWO_CENTRAL  ///< Use O(eps^2) two sided central finite differences
  ,FD_ORDER_TWO_AUTO     ///< Use FD_ORDER_TWO_CENTRAL when not limited by bounds, otherwise use FD_ORDER_TWO
  ,FD_ORDER_FOUR         ///< Use O(eps^4) one sided finite differences (cramped bounds)
  ,FD_ORDER_FOUR_CENTRAL ///< Use O(eps^4) two sided central finite differences
  ,FD_ORDER_FOUR_AUTO    ///< Use FD_ORDER_FOUR_CENTRAL when not limited by bounds, otherwise use FD_ORDER_FOUR
};

/** \brief .
 * \relates DirectionalFiniteDiffCalculator
 */
enum EFDStepSelectType {
  FD_STEP_ABSOLUTE      ///< Use absolute step size <tt>fd_step_size</tt>
  ,FD_STEP_RELATIVE     ///< Use relative step size <tt>fd_step_size * ||xo||inf</tt>
};

/** \brief Simple utility class used to select finite difference
 * derivatives for OutArgs object.
 *
 * \relates DirectionalFiniteDiffCalculator
 */
class SelectedDerivatives {
public:
  /** \brief . */
  SelectedDerivatives() {}
  /** \brief . */
  SelectedDerivatives& supports( ModelEvaluatorBase::EOutArgsDfDp arg, int l )
    { supports_DfDp_.push_back(l); return *this; }
  /** \brief . */
  SelectedDerivatives& supports( ModelEvaluatorBase::EOutArgsDgDp arg, int j, int l )
    { supports_DgDp_.push_back(std::pair<int,int>(j,l)); return *this; }
  // These should be private but I am too lazy to deal with the porting
  // issues of friends ...
  typedef std::list<int> supports_DfDp_t;
  typedef std::list<std::pair<int,int> > supports_DgDp_t;
  supports_DfDp_t supports_DfDp_;
  supports_DgDp_t supports_DgDp_;
};

} // namespace DirectionalFiniteDiffCalculatorTypes

/** \brief Utility calss for computing directional finite differences of a
 * model.
 *
 * This class computes finite difference approximations to the variations:
 * <ul>
 * <li> <tt>df = DfDx*delta_x + sum(DfDp(l)*delta_p(l),l=0...Np)</tt>
 * <li> <tt>dg(j) = sum(DgDx(j)*delta_x,j=0...Ng) + sum(DfDp(j,l)*delta_p(l),j=0...Ng,l=0...Np)</tt>
 * </ul>
 *
 * The client can leave any of the <tt>delta_x</tt> or <tt>delta_p(l)</tt>
 * directions as <tt>NULL</tt> and they will be assumed to be zero.
 *
 * <b>Warning!</b The client should only set parameters using either the
 * parameter list function <tt>setParameterList()</tt> or the more typesafe
 * functions but not a mixture of the two.  The behavior of setting options in
 * two different ways is undefined and is likely to change.
 *
 * ToDo: Finish documentaton!
 */
template<class Scalar>
class DirectionalFiniteDiffCalculator
  : public Teuchos::VerboseObject<DirectionalFiniteDiffCalculator<Scalar> >
  , public Teuchos::ParameterListAcceptor
{
public:
  
  /** \name Public Types */
  //@{

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief . */
  typedef DirectionalFiniteDiffCalculatorTypes::EFDMethodType EFDMethodType;

  /** \brief . */
  typedef DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType EFDStepSelectType;

  /** \brief . */
  typedef DirectionalFiniteDiffCalculatorTypes::SelectedDerivatives SelectedDerivatives;

  //@}

  /** \name Constructors/setup */
  //@{

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDMethodType, fd_method_type );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDStepSelectType, fd_step_select_type );

  /** \brief Pick the size of the finite difference step.
   *
   * If <tt>fd_step_size < 0</tt> then the implementation will try to
   * select it based on the order of method <tt>fd_method_type()</tt>
   * that is selected.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, fd_step_size );

  /** \brief Pick the minimum step size under which the finite difference
   * product will not be computed.
   *
   * If <tt>fd_step_size_min == 0</tt> then the finite difference computation
   * will always be performed.  If <tt>fd_step_size_min < 0</tt> then the
   * minimum step size will be determined internally.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, fd_step_size_min );

  /** \brief . */
  DirectionalFiniteDiffCalculator(
    EFDMethodType               fd_method_type        = DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO
    ,EFDStepSelectType          fd_step_select_type   = DirectionalFiniteDiffCalculatorTypes::FD_STEP_ABSOLUTE
    ,ScalarMag                  fd_step_size          = -1.0
    ,ScalarMag                  fd_step_size_min      = -1.0
    );

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

  /** \name Finite difference functions. */
  //@{

  /** \brief Create an augmented out args object for holding finite difference
   *  objects.
   *
   * <b>Warning!</b> The returned object must only be used with the below
   * functions <tt>calcVariations()</tt> and <tt>calcDerivatives()</tt> and
   * not with the original <tt>model</tt> object directly.
   */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs(
    const ModelEvaluator<Scalar> &model
    ,const SelectedDerivatives   &fdDerivatives
    );

  /** \brief Compute variations using directional finite differences..
   *
   * The computation may fail if a <tt>NaN</tt> or <tt>Inf</tt> is encountered
   * during any of the computations in which case a <tt>NaNInfException</tt>
   * exception will be thrown.  Otherwise the computation should be completed
   * successfully.
   *
   * If the finite difference could not be computed because of cramped bounds
   * then a <tt>CrampedBoundsException</tt> object will be thrown.
   *
   * ToDo: Discuss options!
   */
  void calcVariations(
    const ModelEvaluator<Scalar>                 &model
    ,const ModelEvaluatorBase::InArgs<Scalar>    &basePoint
    ,const ModelEvaluatorBase::InArgs<Scalar>    &directions
    ,const ModelEvaluatorBase::OutArgs<Scalar>   &baseFunctionValues
    ,const ModelEvaluatorBase::OutArgs<Scalar>   &variations
    ) const;

  /** \brief Compute entire derivative objects using finite differences
   */
  void calcDerivatives(
    const ModelEvaluator<Scalar>                 &model
    ,const ModelEvaluatorBase::InArgs<Scalar>    &basePoint
    ,const ModelEvaluatorBase::OutArgs<Scalar>   &baseFunctionValues
    ,const ModelEvaluatorBase::OutArgs<Scalar>   &derivatives
    ) const;

  //@}

private:

  Teuchos::RCP<Teuchos::ParameterList>  paramList_;

  // //////////////////////////////
  // Private static data members

  static const std::string FDMethod_name;
  static const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<EFDMethodType> >
  fdMethodValidator;
  static const std::string FDMethod_default;

  static const std::string FDStepSelectType_name;
  static const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<EFDStepSelectType> >
  fdStepSelectTypeValidator;
  static const std::string FDStepSelectType_default;

  static const std::string FDStepLength_name;
  static const double FDStepLength_default;

};

} // namespace Thyra

// //////////////////////////////
// Implementations

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace Thyra {

namespace DirectionalFiniteDiffCalculatorTypes {

// Undocumented utility class for setting new derivative objects on an OutArgs
// object!  Warning, users should not attempt to play these tricks on their
// own!
//
// Note that because of the design of the OutArgs and OutArgsSetup classes,
// you can only change the list of supported arguments in a subclass of
// ModelEvaluatorBase since OutArgsSetup is a protected type.
template<class Scalar>
class OutArgsCreator : public StateFuncModelEvaluatorBase<Scalar>
{
public:
  // Public functions overridden from ModelEvaulator.
  Teuchos::RCP<const VectorSpaceBase<Scalar> > get_x_space() const
    { TEST_FOR_EXCEPT(true); return Teuchos::null; }
  Teuchos::RCP<const VectorSpaceBase<Scalar> > get_f_space() const
    { TEST_FOR_EXCEPT(true); return Teuchos::null; }
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const
    { TEST_FOR_EXCEPT(true); return ModelEvaluatorBase::InArgs<Scalar>(); }
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const
    { TEST_FOR_EXCEPT(true); return ModelEvaluatorBase::OutArgs<Scalar>(); }
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const
    { TEST_FOR_EXCEPT(true); }
  // Static function that does the magic!
  static ModelEvaluatorBase::OutArgs<Scalar> createOutArgs(
    const ModelEvaluator<Scalar> &model,
    const SelectedDerivatives   &fdDerivatives
    )
    {
      typedef ModelEvaluatorBase MEB;

      const MEB::OutArgs<Scalar> wrappedOutArgs = model.createOutArgs();
      const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
      MEB::OutArgsSetup<Scalar> outArgs;

      outArgs.setModelEvalDescription(
        "DirectionalFiniteDiffCalculator: " + model.description()
        );

      // Start off by supporting everything that the underlying model supports
      // computing!

      outArgs.set_Np_Ng(Np,Ng);
      outArgs.setSupports(wrappedOutArgs);

      // Add support for finite difference DfDp(l) if asked

      const SelectedDerivatives::supports_DfDp_t
        &supports_DfDp = fdDerivatives.supports_DfDp_;
      for(
        SelectedDerivatives::supports_DfDp_t::const_iterator
          itr = supports_DfDp.begin();
        itr != supports_DfDp.end();
        ++itr
        )
      {
        const int l = *itr;
        outArgs.setSupports(MEB::OUT_ARG_DfDp,l,MEB::DERIV_MV_BY_COL);
      }

      const SelectedDerivatives::supports_DgDp_t
        &supports_DgDp = fdDerivatives.supports_DgDp_;
      for(
        SelectedDerivatives::supports_DgDp_t::const_iterator
          itr = supports_DgDp.begin();
        itr != supports_DgDp.end();
        ++itr
        )
      {
        const int j = itr->first;
        const int l = itr->second;
        const bool g_space_j_is_in_core = model.get_g_space(j)->hasInCoreView();
        const bool p_space_l_is_in_core = model.get_p_space(l)->hasInCoreView();
        MEB::DerivativeSupport derivSupport;
        if (g_space_j_is_in_core)
          derivSupport.plus(MEB::DERIV_TRANS_MV_BY_ROW);
        if (p_space_l_is_in_core)
          derivSupport.plus(MEB::DERIV_MV_BY_COL);
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPTION(
          derivSupport.none(), std::logic_error,
          "Error, for the model " << model.description()
          << ", at lease one of the spaces g_space("<<j<<")"
          " or p_space("<<l<<") must be in-core so that they can"
          " act as the domain space for the multi-vector derivative!"
          );
#endif
        outArgs.setSupports(MEB::OUT_ARG_DgDp,j,l,derivSupport);
      }

      return outArgs;

    }
};

} // namespace DirectionalFiniteDiffCalculatorTypes

// Private static data members

template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDMethod_name = "FD Method";
template<class Scalar>
const Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType
  >
>
DirectionalFiniteDiffCalculator<Scalar>::fdMethodValidator
= rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType>(
    Teuchos::tuple<std::string>(
      "order-one"
      ,"order-two"
      ,"order-two-central"
      ,"order-two-auto"
      ,"order-four"
      ,"order-four-central"
      ,"order-four-auto"
      )
    ,Teuchos::tuple<std::string>(
      "Use O(eps) one sided finite differences (cramped bounds)"
      ,"Use O(eps^2) one sided finite differences (cramped bounds)"
      ,"Use O(eps^2) two sided central finite differences"
      ,"Use \"order-two-central\" when not cramped by bounds, otherwise use \"order-two\""
      ,"Use O(eps^4) one sided finite differences (cramped bounds)"
      ,"Use O(eps^4) two sided central finite differences"
      ,"Use \"order-four-central\" when not cramped by bounds, otherwise use \"order-four\""
      )
    ,Teuchos::tuple<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType>(
      Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_ONE
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_CENTRAL
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_AUTO
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_CENTRAL
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO
      )
    ,""
    )
  );
template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDMethod_default = "order-one";

template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDStepLength_name = "FD Step Length";
template<class Scalar>
const double DirectionalFiniteDiffCalculator<Scalar>::FDStepLength_default = -1.0;

template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDStepSelectType_name = "FD Step Select Type";
template<class Scalar>
const Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    Thyra::DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType
    >
  >
DirectionalFiniteDiffCalculator<Scalar>::fdStepSelectTypeValidator
= rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType>(
    Teuchos::tuple<std::string>(
      "Absolute"
      ,"Relative"
      )
    ,Teuchos::tuple<std::string>(
      "Use absolute step size \""+FDStepLength_name+"\""
      ,"Use relative step size \""+FDStepLength_name+"\"*||xo||inf"
      )
    ,Teuchos::tuple<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType>(
      Thyra::DirectionalFiniteDiffCalculatorTypes::FD_STEP_ABSOLUTE
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_STEP_RELATIVE
      )
    ,""
    )
  );
template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDStepSelectType_default = "Absolute";

// Constructors/initializer

template<class Scalar>
DirectionalFiniteDiffCalculator<Scalar>::DirectionalFiniteDiffCalculator(
  EFDMethodType               fd_method_type
  ,EFDStepSelectType          fd_step_select_type
  ,ScalarMag                  fd_step_size
  ,ScalarMag                  fd_step_size_min
  )
  :fd_method_type_(fd_method_type)
  ,fd_step_select_type_(fd_step_select_type)
  ,fd_step_size_(fd_step_size)
  ,fd_step_size_min_(fd_step_size_min)
{}

// Overriden from ParameterListAcceptor

template<class Scalar>
void DirectionalFiniteDiffCalculator<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(paramList.get()==0);
  paramList->validateParameters(*getValidParameters());
  paramList_ = paramList;
  fd_method_type_ = fdMethodValidator->getIntegralValue(
    *paramList_,FDMethod_name,FDMethod_default);
  fd_step_select_type_ = fdStepSelectTypeValidator->getIntegralValue(
    *paramList_,FDStepSelectType_name,FDStepSelectType_default);
  fd_step_size_ = paramList_->get(
    FDStepLength_name,FDStepLength_default);
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::getParameterList()
{
  return paramList_;
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> pl;
  if(pl.get()==NULL) {
    pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set(
      FDMethod_name,FDMethod_default
      ,"The method used to compute the finite differences."
      ,fdMethodValidator
      );
    pl->set(
      FDStepSelectType_name,FDStepSelectType_default
      ,"Method used to select the finite difference step length."
      ,fdStepSelectTypeValidator
      );
    pl->set(
      FDStepLength_name,FDStepLength_default
      ,"The length of the finite difference step to take.\n"
      "A value of < 0.0 means that the step length will be determined automatically."
      );
  }
  return pl;
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DirectionalFiniteDiffCalculator<Scalar>::createOutArgs(
  const ModelEvaluator<Scalar> &model
  ,const SelectedDerivatives   &fdDerivatives
  )
{
  return DirectionalFiniteDiffCalculatorTypes::OutArgsCreator<Scalar>::createOutArgs(
    model,fdDerivatives
    );
}

template<class Scalar>
void DirectionalFiniteDiffCalculator<Scalar>::calcVariations(
  const ModelEvaluator<Scalar>                 &model
  ,const ModelEvaluatorBase::InArgs<Scalar>    &bp         // basePoint
  ,const ModelEvaluatorBase::InArgs<Scalar>    &dir        // directions
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &bfunc      // baseFunctionValues
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &var        // variations
  ) const
{

  using std::string;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_FUNC_TIME_MONITOR(
    string("Thyra::DirectionalFiniteDiffCalculator<")+ST::name()+">::calcVariations(...)"
    );
  
  using std::setw;
  using std::endl;
  using std::right;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ModelEvaluatorBase MEB;
  namespace DFDCT = DirectionalFiniteDiffCalculatorTypes;
  typedef VectorBase<Scalar> V;
  typedef Teuchos::RCP<VectorBase<Scalar> > VectorPtr;
  
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool trace = (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM));
  Teuchos::OSTab tab(out);

  if(out.get() && trace)
    *out << "\nEntering DirectionalFiniteDiffCalculator<Scalar>::calcVariations(...)\n";

  if(out.get() && trace)
    *out
      << "\nbasePoint=\n" << describe(bp,verbLevel)
      << "\ndirections=\n" << describe(dir,verbLevel)
      << "\nbaseFunctionValues=\n" << describe(bfunc,verbLevel)
      << "\nvariations=\n" << describe(var,Teuchos::VERB_LOW);

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    var.isEmpty(), std::logic_error,
    "Error, all of the variations can not be null!"
    );
#endif
  
  //
  // To illustrate the theory behind this implementation consider
  // the generic multi-variable function h(z) <: R^n -> R.  Now let's
  // consider we have the base point zo and the vector v to
  // perturb h(z) along.  First form the function h(zo+a*v) and then
  // let's compute d(h)/d(a) at a = 0:
  // 
  // (1) d(h(zo+a*v))/d(a) at a = 0
  //         = sum( d(h)/d(x(i)) * d(x(i))/d(a), i = 1...n)
  //         = sum( d(h)/d(x(i)) * v(i), i = 1...n)
  //         = d(h)/d(a) * v
  //
  // Now we can approximate (1) using, for example, central differences as:
  // 
  // (2) d(h(zo+a*v))/d(a) at a = 0
  //          \approx ( h(zo+h*v) - h(zo+h*v) ) / (2*h)
  //
  // If we equate (1) and (2) we have the approximation:
  // 
  // (3) d(h)/d(a) * v \approx ( h(zo+h*v) - g(zo+h*v) ) / (2*h)
  // 
  // 

  // /////////////////////////////////////////
  // Validate the input

  // ToDo: Validate input!

  switch(this->fd_method_type()) {
    case DFDCT::FD_ORDER_ONE:
      if(out.get()&&trace) *out<<"\nUsing one-sided, first-order finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_TWO:
      if(out.get()&&trace) *out<<"\nUsing one-sided, second-order finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_TWO_CENTRAL:
      if(out.get()&&trace) *out<<"\nUsing second-order central finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_TWO_AUTO:
      if(out.get()&&trace) *out<<"\nUsing auto selection of some second-order finite difference method ...\n";
      break;
    case DFDCT::FD_ORDER_FOUR:
      if(out.get()&&trace) *out<<"\nUsing one-sided, fourth-order finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_FOUR_CENTRAL:
      if(out.get()&&trace) *out<<"\nUsing fourth-order central finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_FOUR_AUTO:
      if(out.get()&&trace) *out<<"\nUsing auto selection of some fourth-order finite difference method ...\n";
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not get here!
  }

  // ////////////////////////
  // Find the step size

  //
  // Get defaults for the optimal step sizes
  //

  const ScalarMag
    sqrt_epsilon = ST::squareroot(ST::eps()),
    u_optimal_1  = sqrt_epsilon,
    u_optimal_2  = ST::squareroot(sqrt_epsilon),
    u_optimal_4  = ST::squareroot(u_optimal_2);

  ScalarMag
    bp_norm = ST::zero();
  // ToDo above: compute a reasonable norm somehow based on the base-point vector(s)!
  
  ScalarMag
    uh_opt = 0.0;
  switch(this->fd_method_type()) {
    case DFDCT::FD_ORDER_ONE:
      uh_opt = u_optimal_1 * ( fd_step_select_type() == DFDCT::FD_STEP_ABSOLUTE ? 1.0 : bp_norm + 1.0 );
      break;
    case DFDCT::FD_ORDER_TWO:
    case DFDCT::FD_ORDER_TWO_CENTRAL:
    case DFDCT::FD_ORDER_TWO_AUTO:
      uh_opt = u_optimal_2 * ( fd_step_select_type() == DFDCT::FD_STEP_ABSOLUTE ? 1.0 : bp_norm + 1.0 );
      break;
    case DFDCT::FD_ORDER_FOUR:
    case DFDCT::FD_ORDER_FOUR_CENTRAL:
    case DFDCT::FD_ORDER_FOUR_AUTO:
      uh_opt = u_optimal_4 * ( fd_step_select_type() == DFDCT::FD_STEP_ABSOLUTE ? 1.0 : bp_norm + 1.0 );
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not get here!
  }

  if(out.get()&&trace) *out<<"\nDefault optimal step length uh_opt = " << uh_opt << " ...\n";

  //
  // Set the step sizes used.
  //

  ScalarMag
    uh      = this->fd_step_size();

  if( uh < 0 )
    uh = uh_opt;
  else if( fd_step_select_type() == DFDCT::FD_STEP_RELATIVE )
    uh *= (bp_norm + 1.0);

  if(out.get()&&trace) *out<<"\nStep size to be used uh="<<uh<<"\n";

  //
  // Find step lengths that stay in bounds!
  //

  // ToDo: Add logic for variable bounds when needed!

  //
  // Set the actual method being used
  //
  // ToDo: Consider cramped bounds and method order.
  //
  
  DFDCT::EFDMethodType fd_method_type = this->fd_method_type();
  switch(fd_method_type) {
    case DFDCT::FD_ORDER_TWO_AUTO:
      fd_method_type = DFDCT::FD_ORDER_TWO_CENTRAL;
      break;
    case DFDCT::FD_ORDER_FOUR_AUTO:
      fd_method_type = DFDCT::FD_ORDER_FOUR_CENTRAL;
      break;
    default:
      break; // Okay
  }

  //if(out.get()&&trace) *out<<"\nStep size to fit in bounds: uh="<<uh"\n";

  int p_saved;
  if(out.get())
    p_saved = out->precision();

  // ///////////////////////////////////////////////
  // Compute the directional derivatives

  const int Np = var.Np(), Ng = var.Ng();
  
  // Setup storage for perturbed variables
  VectorPtr per_x, per_x_dot;
  std::vector<VectorPtr> per_p(Np);
  MEB::InArgs<Scalar> pp = model.createInArgs();
  pp.setArgs(bp); // Set all args to start with
  if( bp.supports(MEB::IN_ARG_x) ) {
    if( dir.get_x().get() )
      pp.set_x(per_x=createMember(model.get_x_space()));
    else
      pp.set_x(bp.get_x());
  }
  if( bp.supports(MEB::IN_ARG_x_dot) ) {
    if( dir.get_x_dot().get() )
      pp.set_x_dot(per_x_dot=createMember(model.get_x_space()));
    else
      pp.set_x_dot(bp.get_x_dot());
  }
  for( int l = 0; l < Np; ++l ) {
    if( dir.get_p(l).get() )
      pp.set_p(l,per_p[l]=createMember(model.get_p_space(l)));
    else
      pp.set_p(l,bp.get_p(l));
  }
  if(out.get() && trace)
    *out
      << "\nperturbedPoint after initial setup (with some unintialized vectors) =\n"
      << describe(pp,verbLevel);
  
  // Setup storage for perturbed functions
  bool all_funcs_at_base_computed = true;
  MEB::OutArgs<Scalar> pfunc = model.createOutArgs();
  {
    VectorPtr f;
    if( var.supports(MEB::OUT_ARG_f) && (f=var.get_f()).get() ) {
      pfunc.set_f(createMember(model.get_f_space()));
      assign(&*f,ST::zero());
      if(!bfunc.get_f().get()) all_funcs_at_base_computed = false;
    }
    for( int j = 0; j < Ng; ++j ) {
      VectorPtr g_j;
      if( (g_j=var.get_g(j)).get() ) {
        pfunc.set_g(j,createMember(model.get_g_space(j)));
        assign(&*g_j,ST::zero());
        if(!bfunc.get_g(j).get()) all_funcs_at_base_computed = false;
      }
    }
  }
  if(out.get() && trace)
    *out
      << "\nperturbedFunctions after initial setup (with some unintialized vectors) =\n"
      << describe(pfunc,verbLevel);
  
  const int dbl_p = 15;
  if(out.get())
    *out << std::setprecision(dbl_p);
    
  //
  // Compute the weighted sum of the terms
  //
    
  int         num_evals  = 0;
  ScalarMag   dwgt       = ST::zero();
  switch(fd_method_type) {
    case DFDCT::FD_ORDER_ONE: // may only need one eval if f(xo) etc is passed in
      num_evals = 2;
      dwgt      = ScalarMag(1.0);
      break;
    case DFDCT::FD_ORDER_TWO: // may only need two evals if c(xo) etc is passed in
      num_evals = 3;
      dwgt      = ScalarMag(2.0);
      break;
    case DFDCT::FD_ORDER_TWO_CENTRAL:
      num_evals = 2;
      dwgt      = ScalarMag(2.0);
      break;
    case DFDCT::FD_ORDER_FOUR:
      num_evals = 5;
      dwgt      = ScalarMag(12.0);
      break;
    case DFDCT::FD_ORDER_FOUR_CENTRAL:
      num_evals = 4;
      dwgt      = ScalarMag(12.0);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not get here!
  }
  for( int eval_i = 1; eval_i <= num_evals; ++eval_i ) {
    // Set the step constant and the weighting constant
    ScalarMag
      uh_i   = 0.0,
      wgt_i  = 0.0;
    switch(fd_method_type) {
      case DFDCT::FD_ORDER_ONE: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -1.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_TWO: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -3.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +4.0;
            break;
          case 3:
            uh_i  = +2.0;
            wgt_i = -1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_TWO_CENTRAL: {
        switch(eval_i) {
          case 1:
            uh_i  = -1.0;
            wgt_i = -1.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_FOUR: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -25.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +48.0;
            break;
          case 3:
            uh_i  = +2.0;
            wgt_i = -36.0;
            break;
          case 4:
            uh_i  = +3.0;
            wgt_i = +16.0;
            break;
          case 5:
            uh_i  = +4.0;
            wgt_i = -3.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_FOUR_CENTRAL: {
        switch(eval_i) {
          case 1:
            uh_i  = -2.0;
            wgt_i = +1.0;
            break;
          case 2:
            uh_i  = -1.0;
            wgt_i = -8.0;
            break;
          case 3:
            uh_i  = +1.0;
            wgt_i = +8.0;
            break;
          case 4:
            uh_i  = +2.0;
            wgt_i = -1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_TWO_AUTO:
      case DFDCT::FD_ORDER_FOUR_AUTO:
        break; // Okay
      default:
        TEST_FOR_EXCEPT(true);
    }

    if(out.get() && trace)
      *out << "\neval_i="<<eval_i<<", uh_i="<<uh_i<<", wgt_i="<<wgt_i<<"\n";
    Teuchos::OSTab tab(out);
    
    // Compute the weighted term and add it to the sum
    if(uh_i == ST::zero()) {
      MEB::OutArgs<Scalar> bfuncall;
      if(!all_funcs_at_base_computed) {
        // Compute missing functions at the base point
        bfuncall = model.createOutArgs();
        if( pfunc.supports(MEB::OUT_ARG_f) && pfunc.get_f().get() && !bfunc.get_f().get() ) {
          bfuncall.set_f(createMember(model.get_f_space()));
        }
        for( int j = 0; j < Ng; ++j ) {
          if( pfunc.get_g(j).get() && !bfunc.get_g(j).get() ) {
            bfuncall.set_g(j,createMember(model.get_g_space(j)));
          }
        }
        model.evalModel(bp,bfuncall);
        bfuncall.setArgs(bfunc);
      }
      else {
        bfuncall = bfunc;
      }
      // Use functions at the base point
      if(out.get() && trace)
        *out << "\nSetting variations = wgt_i * basePoint ...\n";
      VectorPtr f;
      if( pfunc.supports(MEB::OUT_ARG_f) && (f=var.get_f()).get() ) {
        V_StV(&*f,wgt_i,*bfuncall.get_f());
      }
      for( int j = 0; j < Ng; ++j ) {
        VectorPtr g_j;
        if( (g_j=var.get_g(j)).get() ) {
          V_StV(&*g_j,wgt_i,*bfuncall.get_g(j));
        }
      }
    }
    else {
      if(out.get() && trace)
        *out << "\nSetting perturbedPoint = basePoint + uh_i*uh*direction ...\n";
      // z = zo + uh_i*uh*v
      {
        if( dir.supports(MEB::IN_ARG_x) && dir.get_x().get() )
          V_StVpV(&*per_x,Scalar(uh_i*uh),*dir.get_x(),*bp.get_x());
        if( dir.supports(MEB::IN_ARG_x_dot) && dir.get_x_dot().get() )
          V_StVpV(&*per_x_dot,Scalar(uh_i*uh),*dir.get_x_dot(),*bp.get_x_dot());
        for( int l = 0; l < Np; ++l ) {
          if( dir.get_p(l).get() )
            V_StVpV(&*per_p[l],Scalar(uh_i*uh),*dir.get_p(l),*bp.get_p(l));
        }
      }
      if(out.get() && trace)
        *out << "\nperturbedPoint=\n" << describe(pp,verbLevel);
      // Compute perturbed function values h(zo+uh_i*uh)
      if(out.get() && trace)
        *out << "\nCompute perturbedFunctions at perturbedPoint...\n";
      model.evalModel(pp,pfunc);
      if(out.get() && trace)
        *out << "\nperturbedFunctions=\n" << describe(pfunc,verbLevel);
      // Sum perturbed function values into total variation
      {
        // var_h += wgt_i*perturbed_h
        if(out.get() && trace)
          *out << "\nComputing variations += wgt_i*perturbedfunctions ...\n";
        VectorPtr f;
        if( pfunc.supports(MEB::OUT_ARG_f) && (f=pfunc.get_f()).get() )
          Vp_StV(&*var.get_f(),wgt_i,*f);
        for( int j = 0; j < Ng; ++j ) {
          VectorPtr g_j;
          if( (g_j=pfunc.get_g(j)).get() )
            Vp_StV(&*var.get_g(j),wgt_i,*g_j);
        }
      }
    }
    if(out.get() && trace)
      *out << "\nvariations=\n" << describe(var,verbLevel);
  }
  //
  // Multiply by the scaling factor!
  //
  
  {
    // var_h *= 1.0/(dwgt*uh)
    const Scalar alpha = ST::one()/(dwgt*uh);
    if(out.get() && trace)
      *out
        << "\nComputing variations *= (1.0)/(dwgt*uh),"
        << " where (1.0)/(dwgt*uh) = (1.0)/("<<dwgt<<"*"<<uh<<") = "<<alpha<<" ...\n";
    VectorPtr f;
    if( var.supports(MEB::OUT_ARG_f) && (f=var.get_f()).get() )
      Vt_S(&*f,alpha);
    for( int j = 0; j < Ng; ++j ) {
      VectorPtr g_j;
      if( (g_j=var.get_g(j)).get() )
        Vt_S(&*g_j,alpha);
    }
    if(out.get() && trace)
      *out << "\nFinal variations=\n" << describe(var,verbLevel);
  }
  
  if(out.get())
    *out << std::setprecision(p_saved);

  if(out.get() && trace)
    *out << "\nLeaving DirectionalFiniteDiffCalculator<Scalar>::calcVariations(...)\n";
  
}

template<class Scalar>
void DirectionalFiniteDiffCalculator<Scalar>::calcDerivatives(
  const ModelEvaluator<Scalar>                 &model
  ,const ModelEvaluatorBase::InArgs<Scalar>    &bp      // basePoint
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &bfunc   // baseFunctionValues
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &deriv   // derivatives
  ) const
{

  using std::string;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_FUNC_TIME_MONITOR(
    string("Thyra::DirectionalFiniteDiffCalculator<")+ST::name()+">::calcDerivatives(...)"
    );

  typedef ModelEvaluatorBase MEB;
  typedef Teuchos::RCP<VectorBase<Scalar> > VectorPtr;
  typedef Teuchos::RCP<MultiVectorBase<Scalar> > MultiVectorPtr;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool trace = (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM));
  Teuchos::OSTab tab(out);

  if(out.get() && trace)
    *out << "\nEntering DirectionalFiniteDiffCalculator<Scalar>::calcDerivatives(...)\n";

  if(out.get() && trace)
    *out
      << "\nbasePoint=\n" << describe(bp,verbLevel)
      << "\nbaseFunctionValues=\n" << describe(bfunc,verbLevel)
      << "\nderivatives=\n" << describe(deriv,Teuchos::VERB_LOW);
  
  //
  // We will only compute finite differences w.r.t. p(l) for now
  //
  const int Np = bp.Np(), Ng = bfunc.Ng();
  MEB::InArgs<Scalar> dir = model.createInArgs();
  MEB::OutArgs<Scalar> var = model.createOutArgs();
  MultiVectorPtr DfDp_l;
  std::vector<MEB::DerivativeMultiVector<Scalar> > DgDp_l(Ng);
  std::vector<VectorPtr> var_g(Ng);
  for( int l = 0; l < Np; ++l ) {
    if(out.get() && trace)
      *out << "\nComputing derivatives for parameter subvector p("<<l<<") ...\n";
    Teuchos::OSTab tab(out);
    //
    // Load up OutArgs var object of derivative variations to compute
    //
    bool hasDerivObject = false;
    // DfDp(l)
    if (
      !deriv.supports(MEB::OUT_ARG_DfDp,l).none()
      && !deriv.get_DfDp(l).isEmpty()
      )
    {
      hasDerivObject = true;
      std::ostringstream name; name << "DfDp("<<l<<")";
      DfDp_l = get_mv(deriv.get_DfDp(l),name.str(),MEB::DERIV_MV_BY_COL);
    }
    else {
      DfDp_l = Teuchos::null;
    }
    // DgDp(j=1...Ng,l)
    for ( int j = 0; j < Ng; ++j ) {
      if (
        !deriv.supports(MEB::OUT_ARG_DgDp,j,l).none()
        &&
        !deriv.get_DgDp(j,l).isEmpty()
        )
      {
        hasDerivObject = true;
        std::ostringstream name; name << "DgDp("<<j<<","<<l<<")";
        DgDp_l[j] = get_dmv(deriv.get_DgDp(j,l),name.str());
        if( DgDp_l[j].getMultiVector().get() && !var_g[j].get() )
        {
          // Need a temporary vector for the variation for g(j)
          var_g[j] = createMember(model.get_g_space(j));
        }
      }
      else{
        DgDp_l[j] = MEB::DerivativeMultiVector<Scalar>();
        var_g[j] = Teuchos::null;
      }
    }
    //
    // Compute derivative variations by finite differences
    //
    if (hasDerivObject) {
      VectorPtr e_i = createMember(model.get_p_space(l));
      dir.set_p(l,e_i);
      assign(&*e_i,ST::zero());
      const int np_l = e_i->space()->dim();
      for( int i = 0 ; i < np_l; ++ i ) {
        if(out.get() && trace)
          *out << "\nComputing derivatives for single variable p("<<l<<")("<<i<<") ...\n";
        Teuchos::OSTab tab(out);
        if(DfDp_l.get()) var.set_f(DfDp_l->col(i)); // Compute in place!
        for(int j = 0; j < Ng; ++j) {
          MultiVectorPtr DgDp_j_l;
          if( (DgDp_j_l=DgDp_l[j].getMultiVector()).get() ) {
            var.set_g(j,var_g[j]); // Computes d(g(j))/d(p(l)(i))
          }
        }
        set_ele(i,ST::one(),&*e_i);
        this->calcVariations(
          model,bp,dir,bfunc,var
          );
        set_ele(i,ST::zero(),&*e_i);
        if (DfDp_l.get()) var.set_f(Teuchos::null);
        for (int j = 0; j < Ng; ++j) {
          MultiVectorPtr DgDp_j_l;
          if ( (DgDp_j_l=DgDp_l[j].getMultiVector()).get() ) {
            ConstDetachedVectorView<Scalar> _var_g_j(var_g[j]); //d(g(j))/d(p(l)(i))
            if (DgDp_l[j].getOrientation()==MEB::DERIV_TRANS_MV_BY_ROW) {
              DetachedMultiVectorView<Scalar> dv_DgDp_trans_j_l_i(*DgDp_j_l,Range1D(i,i),Range1D());
              const int ng = _var_g_j.subDim();
              for( int k = 0; k < ng; ++k )
                dv_DgDp_trans_j_l_i(0,k) = _var_g_j(k);
            }
            else {
              DetachedMultiVectorView<Scalar> dv_DgDp_j_l_i(*DgDp_j_l,Range1D(),Range1D(i,i));
              const int ng = _var_g_j.subDim();
              for( int k = 0; k < ng; ++k )
                dv_DgDp_j_l_i(k,0) = _var_g_j(k);
            }
          }
        }
      }
      dir.set_p(l,Teuchos::null);
    }
  }
    
  if(out.get() && trace)
    *out
      << "\nderivatives=\n" << describe(deriv,verbLevel);
  
  if(out.get() && trace)
    *out << "\nLeaving DirectionalFiniteDiffCalculator<Scalar>::calcDerivatives(...)\n";

}

} // end namespace Thyra

#endif	// THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_HPP
