#ifndef MUELU_AMESOS_SMOOTHER_HPP
#define MUELU_AMESOS_SMOOTHER_HPP

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_AMESOS
#include "Amesos_BaseSolver.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"

namespace MueLu {

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class AmesosSmoother
  @brief Class that encapsulates Amesos direct solvers.

  This class creates an Amesos preconditioner factory.  The factory creates a direct solver based on the
  type and ParameterList passed into the constructor.  See the constructor for more information.
*/

  template<class ScalarType,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class AmesosSmoother : public SmootherPrototype<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    //! amesos-specific key phrase that denote smoother type (same as SmootherBase::Type_)
    std::string amesosType_;
    //! pointer to Amesos solver object
    RCP<Amesos_BaseSolver> prec_;
    //! matrix operator 
    Teuchos::RCP<Operator> A_;
    //! parameter list that is used by Amesos internally
    Teuchos::ParameterList list_;
    //! Problem that Amesos uses internally.
    RCP<Epetra_LinearProblem> AmesosLinearProblem;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

    */
    AmesosSmoother(std::string const & type, Teuchos::ParameterList & list)
      : amesosType_(type), list_(list), out_(this->getOStream())
    {
      MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating a new Amesos smoother" << std::endl;
      SmootherBase::SetType(type);
    }

    //! Destructor
    virtual ~AmesosSmoother() {}
    //@}

    //! @name Set/Get methods

    //@{
    void SetNIts(LO const &nIts) {
      throw(Exceptions::RuntimeError("Only one iteration of Amesos solve is supported."));
    }

    LO GetNIts() {
      return 1;
    }
    //@}

    void Setup(Level &level) {
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "AmesosSmoother::Setup()" << std::endl;
      SmootherPrototype::IsSetup(true);
      A_ = level.GetA();
      RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
      AmesosLinearProblem = rcp(new Epetra_LinearProblem());
      AmesosLinearProblem->SetOperator(&*epA); //FIXME RCP probably has a safer way to do this
      Amesos factory;
      prec_ = rcp(factory.Create(amesosType_, *AmesosLinearProblem));
      if (prec_ == Teuchos::null) {
        std::string msg = "Amesos::Create: factorization type '" + amesosType_ + "' is not supported";
        throw(Exceptions::RuntimeError(msg));
      }
      prec_->SetParameters(list_);
      int rv = prec_->NumericFactorization();
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "Amesos_BaseSolver::NumericFactorization return value of " + buf.str(); 
        throw(Exceptions::RuntimeError(msg));
      }
    }

    void Apply(RCP<MultiVector> x, RCP<MultiVector> const rhs, bool InitialGuessIsZero=false)
    {
      if (!SmootherPrototype::IsSetup())
        throw(Exceptions::RuntimeError("Setup has not been called"));

      RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(x);
      RCP<Epetra_MultiVector> epRhs = Utils::MV2NonConstEpetraMV(rhs);//FIXME Amesos wants a nonconst rhs!

      AmesosLinearProblem->SetLHS(&*epX); //FIXME RCP probably has a safer way to do this
      AmesosLinearProblem->SetRHS(&*epRhs); //FIXME RCP probably has a safer way to do this

      prec_->Solve();
    }

    void Print(std::string prefix) {
      throw(Exceptions::NotImplemented("AmesosSmoother::Print is not implemented"));
    }

    RCP<SmootherPrototype> Copy()
    {
      return rcp(new AmesosSmoother(*this) );
    }

    void CopyParameters(RCP<SmootherPrototype> source)
    {
      RCP<AmesosSmoother> amesosSmoo = rcp_dynamic_cast<AmesosSmoother>(source);
      //TODO check if dynamic cast fails
      amesosType_ = amesosSmoo->amesosType_;
      prec_ = amesosSmoo->prec_;
      A_ = amesosSmoo->A_;
      list_ = amesosSmoo->list_;
    }

  }; //class AmesosSmoother

} //namespace MueLu

#define MUELU_AMESOS_SMOOTHER_SHORT

#endif //ifdef HAVE_MUELU_AMESOS

#endif //ifndef MUELU_AMESOS_SMOOTHER_HPP
