// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceFunctional.hpp"

  %}


namespace Sundance
{
class Functional
  {
  public:
    /** */
    Functional(){;}

    /** */
    Functional(const Sundance::Mesh& mesh, 
      const Sundance::Expr& integral, 
      const TSFExtended::VectorType<double>& vecType);

    /** */
    Functional(const Sundance::Mesh& mesh, 
      const Sundance::Expr& integral, 
      const Sundance::Expr& essentialBC,
      const TSFExtended::VectorType<double>& vecType);

    /** */
    NonlinearProblem
    nonlinearVariationalProb(const Sundance::Expr& var,
                             const Sundance::Expr& varEvalPts,
                             const Sundance::Expr& unk,
                             const Sundance::Expr& unkEvalPts,
                             const Sundance::Expr& fixed,
                             const Sundance::Expr& fixedEvalPts) const ;


    /** */
    FunctionalEvaluator evaluator(const Sundance::Expr& var,
                                  const Sundance::Expr& varEvalPts,
                                  const Sundance::Expr& fixed,
                                  const Sundance::Expr& fixedEvalPts) const ;


    /** */
    FunctionalEvaluator evaluator(const Sundance::Expr& var,
                                  const Sundance::Expr& varEvalPts) const ;

    /** */
    const Sundance::Mesh& mesh() const ;
};

class FunctionalEvaluator 
  : public TSFExtended::ParameterControlledObjectWithVerbosity<FunctionalEvaluator>
{
public:
  /** */
  FunctionalEvaluator();

  /** */
  FunctionalEvaluator(const Sundance::Mesh& mesh, 
    const Sundance::Expr& integral,
    const Teuchos::ParameterList& verbParams = *defaultVerbParams());
  /** */
  FunctionalEvaluator(const Sundance::Mesh& mesh, 
    const Sundance::Expr& integral,
    const Sundance::Expr& bcs,
    const Sundance::Expr& var,
    const Sundance::Expr& varEvalPts,
    const TSFExtended::VectorType<double>& vectorType,
    const Teuchos::ParameterList& verbParams = *defaultVerbParams());
  /** */
  FunctionalEvaluator(const Sundance::Mesh& mesh, 
    const Sundance::Expr& integral,
    const Sundance::Expr& bcs,
    const Sundance::Expr& vars,
    const Sundance::Expr& varEvalPts,
    const Sundance::Expr& fields,
    const Sundance::Expr& fieldValues,
    const TSFExtended::VectorType<double>& vectorType,
    const Teuchos::ParameterList& verbParams = *defaultVerbParams());


  /** */
  double evaluate() const ;

  /** */
  Sundance::Expr evalGradient(double& value) const ;

  /** */
  double fdGradientCheck(double h) const ;
};


}
