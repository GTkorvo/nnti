/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STDMATHFUNCTORS_H
#define SUNDANCE_STDMATHFUNCTORS_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceUnaryFunctor.hpp"

namespace SundanceUtils
{
  using namespace Teuchos;

  /** */
  class PowerFunctor : public UnaryFunctor
  {
  public:
    /** */
    PowerFunctor(const double& p);
    
    /** Evaluate power function and deriv at an array of values */ 
    virtual void eval1(const double* const x, 
              int nx, 
              double* f, 
              double* df) const ;
    /** Evaluate power function at an array of values */ 
    virtual void eval0(const double* const x, int nx, double* f) const ;

    /** Evaluate power function and first two derivs at an array of values */
    virtual void eval2(const double* const x, 
                      int nx, 
                      double* f, 
                      double* df_dx,
                      double* d2f_dxx) const ;
  private:
    double p_;
  };


  using std::string;
  using std::ostream;

  SUNDANCE_UNARY_FUNCTOR(reciprocal, StdReciprocal, "reciprocal function", 
                         NonzeroDomain(), 1.0/x[i], -f[i]*f[i], -2.0*df[i]/x[i]);

  SUNDANCE_UNARY_FUNCTOR(fabs, StdFabs, "absolute value", UnboundedDomain(), ::fabs(x[i]), ((x[i]>=0.0) ? x[i] : -x[i]), 0.0);

  SUNDANCE_UNARY_FUNCTOR(sign, StdSign, "sign function", UnboundedDomain(), 
                         ((x[i]>0.0) ? 1.0 : ( (x[i]<0.0) ? -1.0 : 0.0)), 
                         0.0, 0.0);

  SUNDANCE_UNARY_FUNCTOR(exp, StdExp, "exponential function", UnboundedDomain(), ::exp(x[i]), f[i], f[i]);

  SUNDANCE_UNARY_FUNCTOR(log, StdLog, "logarithm", PositiveDomain(), ::log(x[i]), 1.0/x[i], -df[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(sqrt, StdSqrt, "square root", PositiveDomain(), ::sqrt(x[i]), 0.5/f[i], -0.5*df[i]/x[i]);

  SUNDANCE_UNARY_FUNCTOR(sin, StdSin, "sine function", UnboundedDomain(), ::sin(x[i]), ::cos(x[i]), -f[i]);

  SUNDANCE_UNARY_FUNCTOR(cos, StdCos, "cosine function", UnboundedDomain(), ::cos(x[i]), -::sin(x[i]), -f[i]);

  SUNDANCE_UNARY_FUNCTOR(tan, StdTan, "tangent function", UnboundedDomain(),
                         ::tan(x[i]), 1.0 + f[i]*f[i], 2.0*f[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(asin, StdASin, "inverse sine", 
                         BoundedDomain(-1.0, 1.0),
                         ::asin(x[i]), 1.0/::sqrt(1.0-x[i]*x[i]),
                         x[i]*df[i]*df[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(acos, StdACos, "inverse cosine",
                         BoundedDomain(-1.0, 1.0), 
                         ::acos(x[i]), -1.0/::sqrt(1.0-x[i]*x[i]),
                         x[i]*df[i]*df[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(atan, StdATan, "inverse tangent", 
                         UnboundedDomain(),
                         ::atan(x[i]), 1.0/(1.0 + x[i]*x[i]),
                         -2.0*x[i]*df[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(sinh, StdSinh, "hyperbolic sine",
                         UnboundedDomain(),
                         ::sinh(x[i]), ::cosh(x[i]), f[i]);

  SUNDANCE_UNARY_FUNCTOR(cosh, StdCosh, "hyperbolic cosine",
                         UnboundedDomain(),
                         ::cosh(x[i]), ::sinh(x[i]), f[i]);

  SUNDANCE_UNARY_FUNCTOR(tanh, StdTanh, "hyperbolic tangent",
                         UnboundedDomain(),
                         ::tanh(x[i]), 1.0 - f[i]*f[i], -2.0*f[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(asinh, StdASinh, "inverse hyperbolic sine",
                         UnboundedDomain(),
                         ::asinh(x[i]), 1.0/::sqrt(1.0 + x[i]*x[i]),
                         -x[i]*df[i]*df[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(acosh, StdACosh, "inverse hyperbolic cosine",
                         LowerBoundedDomain(1.0),
                         ::acosh(x[i]), 1.0/::sqrt(x[i]*x[i]-1.0),
                         -x[i]*df[i]*df[i]*df[i]);

  SUNDANCE_UNARY_FUNCTOR(atanh, StdATanh, "inverse hyperbolic tangent",
                         BoundedDomain(-1.0, 1.0), 
                         ::atanh(x[i]), 1.0/(1.0 - x[i]*x[i]),
                         2.0*x[i]*df[i]*df[i]);


}

#endif
