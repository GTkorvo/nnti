#ifndef ML_LINEARCOMBINATION_H
#define ML_LINEARCOMBINATION_H

#include "MLAPI_BaseLinearCombination.h"

namespace MLAPI {

class BaseOperator;
class MultiVector;

// ============================================================================ 
class LinearCombinationAdd : public BaseLinearCombination
{
public:
  LinearCombinationAdd(const BaseLinearCombination& left,
                       const BaseLinearCombination& right) :
    left_(left),
    right_(right)
  {}

  const Space GetVectorSpace() const
  {
    return(left_.GetVectorSpace());
  }

  void Update(MultiVector& v) const
  {
    left_.Update(v);
    right_.Update(v);
  }

  void Set(MultiVector& v) const
  {
    left_.Set(v);
    right_.Update(v);
  }

private:
  const BaseLinearCombination& left_;
  const BaseLinearCombination& right_;
};

// ============================================================================ 
class LinearCombinationMixed : public BaseLinearCombination
{
public:

  LinearCombinationMixed(const BaseLinearCombination& left, 
                         const MultiVector& right, double alpha) :
    left_(left),
    right_(right),
    alpha_(alpha)
  {}

  const Space GetVectorSpace() const;

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;
  
private:
  const BaseLinearCombination& left_;
  const MultiVector            right_;
  double                       alpha_;
};

// ============================================================================
class LinearCombinationScaled : public BaseLinearCombination 
{
public:
  LinearCombinationScaled(const BaseLinearCombination& left, double scalar) :
    left_(left),
    scalar_(scalar)
  {}

  const Space GetVectorSpace() const;

  void Set(MultiVector& v) const;

  void Update(MultiVector& v) const;

private:
  const BaseLinearCombination& left_;
  double                       scalar_;
};

// ============================================================================ 
// scaled vector, ScaledMultiVector = alpha * MultiVector
// ============================================================================ 
class MultiVectorScaled : public BaseLinearCombination
{
public:
  MultiVectorScaled(const MultiVector& vector, const double alpha) :
    vector_(vector),
    alpha_(alpha)
  {}

  const Space GetVectorSpace() const;

  const MultiVector& GetMultiVector() const
  {
    return(vector_);
  }

  const double GetScalar() const
  {
    return(alpha_);
  }

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:
  const MultiVector vector_;
  double            alpha_;
};

// ============================================================================ 
class MultiVectorCombination : public BaseLinearCombination
{
public:
  MultiVectorCombination(const double alpha, 
                         const MultiVector x,
                         const double beta,
                         const MultiVector y) :
    alpha_(alpha),
    x_(x),
    beta_(beta),
    y_(y)
  {}

  const Space GetVectorSpace() const;

  const MultiVector GetLeftMultiVector() const
  {
    return(x_);
  }

  const double GetLeftScalar() const
  {
    return(alpha_);
  }

  const MultiVector GetRightMultiVector() const
  {
    return(y_);
  }

  const double GetRightScalar() const
  {
    return(beta_);
  }

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:
  const MultiVector x_;
  const MultiVector y_;
  double alpha_, beta_;
};

// ============================================================================ 
// v = A * x
// ============================================================================ 
class BaseOperatorTimesMultiVector : public BaseLinearCombination
{
public:
  BaseOperatorTimesMultiVector(const BaseOperator& A,
                               const MultiVector& x) :
    A_(A),
    x_(x)
  {}

  const Space GetVectorSpace() const;

  const BaseOperator& GetBaseOperator() const
  {
    return(A_);
  }

  const MultiVector& GetMultiVector() const
  {
    return(x_);
  }

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private: 
  const BaseOperator& A_;
  const MultiVector   x_;
};

// ============================================================================ 
// v += alpha * b + beta * A * x
// ============================================================================ 
class Residual : public BaseLinearCombination
{
public:

  Residual(double alpha, const MultiVector& b, double beta, 
           const BaseOperator& A, const MultiVector& x) :
    alpha_(alpha),
    beta_(beta),
    b_(b),
    A_(A),
    x_(x)
  {}

  const Space GetVectorSpace() const;

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:

  const BaseOperator& A_;
  const MultiVector  b_;
  const MultiVector  x_;
  double             alpha_;
  double             beta_;
};

} // namespace MLAPI

#endif 
