/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LINEARCOMBINATIONIMPL_HPP
#define PLAYA_LINEARCOMBINATIONIMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationDecl.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif


namespace PlayaExprTemplates
{

using Playa::Out;
using Playa::Tabs;


/* -------- methods of OpTimesLC ------ */

template <class Scalar, class Node> inline
OpTimesLC<Scalar, Node>::OpTimesLC(const Scalar& alpha, 
  const Node& x)
  : alpha_(alpha), op_(), x_(x) 
{;}

template <class Scalar, class Node> inline
OpTimesLC<Scalar, Node>
::OpTimesLC(const Scalar& alpha,
  const Playa::LinearOperator<Scalar>& op, 
  const Node& x)
  : alpha_(alpha), op_(op), x_(x) 
{;}
  
  
template <class Scalar, class Node> inline
void OpTimesLC<Scalar, Node>::evalInto(Playa::Vector<Scalar>& result) const
{
  if (op_.ptr().get() != 0)
  {
    op_.apply(x_.eval(), result);
  }
  else
  {
    x_.evalInto(result);
  }
  if (alpha_ != one()) result.scale(alpha_);
}

template <class Scalar, class Node> inline
void OpTimesLC<Scalar, Node>::addInto(Playa::Vector<Scalar>& result,
  LCSign sign) const
{
  if (op_.ptr().get() != 0)
  {
    Vector<Scalar> tmp;
    op_.apply(x_.eval(), tmp);
    result.update(sign*alpha_, tmp);
  }
  else
  {
    result.update(sign*alpha_, x_.eval());
  }
} 

template <class Scalar, class Node> inline
Playa::Vector<Scalar> OpTimesLC<Scalar, Node>::eval() const 
{
  Playa::Vector<Scalar> result;
  if (op_.ptr().get() != 0)
  {
    result = op_.range().createMember();
    op_.apply(x_.eval(), result);
  }
  else
  {
    result = x_.eval();
  }

  if (alpha_ != one()) result.scale(alpha_);    
  return result;
}


template <class Scalar, class Node> inline
bool OpTimesLC<Scalar, Node>::containsVector(const VectorBase<Scalar>* vec) const 
{return x_.containsVector(vec);}





  
/* ------------------------ methods of LC2 --------------------------- */
  
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, Node1, Node2>::LC2(const Node1& x1, const Node2& x2, LCSign sign)
  : x1_(x1), x2_(x2), sign_(sign) 
{;}

template <class Scalar, class Node1, class Node2> inline
bool LC2<Scalar, Node1, Node2>::containsVector(const VectorBase<Scalar>* vec) const
{return x1_.containsVector(vec) || x2_.containsVector(vec);}

template <class Scalar, class Node1, class Node2> inline
void LC2<Scalar, Node1, Node2>::evalInto(Playa::Vector<Scalar>& result) const
{
  Tabs tab;
  x1_.evalInto(result);
  x2_.addInto(result, sign_);
} 

template <class Scalar, class Node1, class Node2> inline
void LC2<Scalar, Node1, Node2>::addInto(Playa::Vector<Scalar>& result,
  PlayaExprTemplates::LCSign sign) const
{
  x1_.addInto(result, sign);
  if (sign_*sign < 0) x2_.addInto(result, LCSubtract);
  else x2_.addInto(result, LCAdd);
}

template <class Scalar, class Node1, class Node2> inline
Playa::Vector<Scalar> LC2<Scalar, Node1, Node2>::eval() const
{
  Playa::Vector<Scalar> result = x1_.eval();
  x2_.addInto(result, sign_);
  return result;
}
}

namespace Playa
{
using PlayaExprTemplates::OpTimesLC;
using PlayaExprTemplates::LC2;
using PlayaExprTemplates::LCAdd;
using PlayaExprTemplates::LCSubtract;

/* ------------------------ global methods ----------------------- */


/*======================================================================
 *
 *    scalar times vector
 *
 *======================================================================*/

/* scalar * vec */
template <class Scalar> inline
OpTimesLC<Scalar, Vector<Scalar> > operator*(const Scalar& alpha, 
  const Vector<Scalar>& x)
{
  return OpTimesLC<Scalar, Vector<Scalar> >(alpha, x);
} 

/* vec * scalar */
template <class Scalar> inline
OpTimesLC<Scalar, Vector<Scalar> > operator*(const Vector<Scalar>& x, 
  const Scalar& alpha)
{
  return OpTimesLC<Scalar, Vector<Scalar> >(alpha, x);
}


/*======================================================================
 *
 *    scalar times OpTimesLC
 *
 *======================================================================*/

/* scalar * OpTimesLC */
template <class Scalar, class Node> inline
OpTimesLC<Scalar, Node> 
operator*(const Scalar& alpha, 
  const OpTimesLC<Scalar, Node>& x)
{
  return OpTimesLC<Scalar, Node>(alpha * x.alpha(), x.op(), x.node());
}

/* OpTimesLC * scalar */
template <class Scalar, class Node> inline
OpTimesLC<Scalar, Node> 
operator*(const OpTimesLC<Scalar, Node>& x, const Scalar& alpha)
{
  return alpha * x;
}


/*======================================================================
 *
 *    scalar times LC2
 *
 *======================================================================*/

/* scalar * LC2 */
template <class Scalar, class Node1, class Node2> inline
OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
operator*(const Scalar& alpha, 
  const LC2<Scalar, Node1, Node2>& x)
{
  return OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> >(alpha, x);
}

/* LC2 * scalar */
template <class Scalar, class Node1, class Node2> inline
OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
operator*(const LC2<Scalar, Node1, Node2>& x, const Scalar& alpha)
{
  return alpha * x;
}
  


/*======================================================================
 *
 *    operator times [vectors, OpTimesLC, LC2]
 *
 *======================================================================*/

/* op * vec */
template <class Scalar> inline
OpTimesLC<Scalar, Vector<Scalar> > 
operator*(const LinearOperator<Scalar>& op, 
  const Vector<Scalar>& x)
{
  return OpTimesLC<Scalar, Vector<Scalar> >(Teuchos::ScalarTraits<Scalar>::one(), op, x);
}


/* op * OpTimesLC */
template <class Scalar, class Node> inline
OpTimesLC<Scalar, Node> 
operator*(const LinearOperator<Scalar>& op, 
  const OpTimesLC<Scalar, Node>& x)
{
  TEST_FOR_EXCEPTION(op.ptr().get()==0, std::runtime_error,
    "null operator in LinearOperator * ( OpTimesLC )");
  if (x.op().ptr().get()==0)
  {
    return OpTimesLC<Scalar, Node>(x.alpha(), op, x.node());
  }
  else
  {
    return OpTimesLC<Scalar, Node>(x.alpha(), op * x.op(), x.node());
  }
}


/* op * LC2 */
template <class Scalar, class Node1, class Node2> inline
OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
operator*(const LinearOperator<Scalar>& op, 
  const LC2<Scalar, Node1, Node2>& x)
{
  return OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> >(Teuchos::ScalarTraits<Scalar>::one(), op, x);
}


/*======================================================================
 *
 *    add/subtract vector, vector
 *
 *======================================================================*/
  
/* vec + vec */
template <class Scalar> inline
LC2<Scalar, Vector<Scalar>, Vector<Scalar> >
operator+(const Vector<Scalar>& x1, 
  const Vector<Scalar>& x2)
{
  return LC2<Scalar, Vector<Scalar>, Vector<Scalar> >(x1, x2);
}
  
/* vec - vec */
template <class Scalar> inline
LC2<Scalar, Vector<Scalar>, Vector<Scalar> >
operator-(const Vector<Scalar>& x1, 
  const Vector<Scalar>& x2)
{
  return LC2<Scalar, Vector<Scalar>, Vector<Scalar> >(x1, x2, LCSubtract);
}

/*======================================================================
 *
 *    add/subtract vector, OpTimesLC
 *
 *======================================================================*/


/* vec + OpTimesLC */
template <class Scalar, class Node> inline
LC2<Scalar, Vector<Scalar>, OpTimesLC<Scalar, Node> >
operator+(const Vector<Scalar>& x1, 
  const OpTimesLC<Scalar, Node>& x2)
{
  return LC2<Scalar, Vector<Scalar>, OpTimesLC<Scalar, Node> >(x1, x2);
}


/* vec - OpTimesLC */
template <class Scalar, class Node> inline
LC2<Scalar, Vector<Scalar>, OpTimesLC<Scalar, Node> >
operator-(const Vector<Scalar>& x1, 
  const OpTimesLC<Scalar, Node>& x2)
{
  return LC2<Scalar, Vector<Scalar>, OpTimesLC<Scalar, Node> >(x1, x2, 
    LCSubtract);
}

/* OpTimesLC + vec */
template <class Scalar, class Node> inline
LC2<Scalar, OpTimesLC<Scalar, Node>, Vector<Scalar> >
operator+(const OpTimesLC<Scalar, Node>& x1, 
  const Vector<Scalar>& x2)
{
  return LC2<Scalar, OpTimesLC<Scalar, Node>, Vector<Scalar> >(x1, x2);
}
  
/* OpTimesLC - vec */
template <class Scalar, class Node> inline
LC2<Scalar, OpTimesLC<Scalar, Node>, Vector<Scalar> >
operator-(const OpTimesLC<Scalar, Node>& x1, 
  const Vector<Scalar>& x2)
{
  return LC2<Scalar, OpTimesLC<Scalar, Node>, Vector<Scalar> >(x1, x2,
    LCSubtract);
}

  
/*======================================================================
 *
 *    add/subtract OpTimesLC, OpTimesLC
 *
 *======================================================================*/
  
/* OpTimesLC + OpTimesLC */
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, OpTimesLC<Scalar, Node1>, OpTimesLC<Scalar, Node2> >
operator+(const OpTimesLC<Scalar, Node1>& x1, 
  const OpTimesLC<Scalar, Node2>& x2)
{
  return LC2<Scalar, OpTimesLC<Scalar, Node1>, 
    OpTimesLC<Scalar, Node2> >(x1, x2);
}
  
/* OpTimesLC - OpTimesLC */
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, OpTimesLC<Scalar, Node1>, OpTimesLC<Scalar, Node2> >
operator-(const OpTimesLC<Scalar, Node1>& x1, 
  const OpTimesLC<Scalar, Node2>& x2)
{
  return LC2<Scalar, OpTimesLC<Scalar, Node1>, 
    OpTimesLC<Scalar, Node2> >(x1, x2, LCSubtract);
}
  

  
/*======================================================================
 *
 *    add/subtract Vector, LC2
 *
 *======================================================================*/

  
/* vec + LC2 */
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, Vector<Scalar>, LC2<Scalar, Node1, Node2> >
operator+(const Vector<Scalar>& x1, 
  const LC2<Scalar, Node1, Node2>& x2)
{
  return LC2<Scalar, Vector<Scalar>, LC2<Scalar, Node1, Node2> >(x1, x2);
}

/* vec - LC2 */
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, Vector<Scalar>, LC2<Scalar, Node1, Node2> >
operator-(const Vector<Scalar>& x1, 
  const LC2<Scalar, Node1, Node2>& x2)
{
  return LC2<Scalar, Vector<Scalar>, LC2<Scalar, Node1, Node2> >(x1, x2,
    LCSubtract);
}


/* LC2 + vec */
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, LC2<Scalar, Node1, Node2>, Vector<Scalar> >
operator+(const LC2<Scalar, Node1, Node2>& x1, 
  const Vector<Scalar>& x2)
{
  return LC2<Scalar, LC2<Scalar, Node1, Node2>, Vector<Scalar> >(x1, x2);
}

/* LC2 - vec */
template <class Scalar, class Node1, class Node2> inline
LC2<Scalar, LC2<Scalar, Node1, Node2>, Vector<Scalar> >
operator-(const LC2<Scalar, Node1, Node2>& x1, 
  const Vector<Scalar>& x2)
{
  return LC2<Scalar, LC2<Scalar, Node1, Node2>, Vector<Scalar> >(x1, x2,
    LCSubtract);
}


/*======================================================================
 *
 *    add/subtract OpTimesLC, LC2
 *
 *======================================================================*/


/* OpTimesLC + LC2 */
template <class Scalar, class Node0, class Node1, class Node2> inline
LC2<Scalar, OpTimesLC<Scalar, Node0>, LC2<Scalar, Node1, Node2> > 
operator+(const OpTimesLC<Scalar, Node0>& x1, 
  const LC2<Scalar, Node1, Node2>& x2)
{
  return LC2<Scalar, OpTimesLC<Scalar, Node0>,
    LC2<Scalar, Node1, Node2> >(x1, x2);
}

/* OpTimesLC - LC2 */
template <class Scalar, class Node0, class Node1, class Node2> inline
LC2<Scalar, OpTimesLC<Scalar, Node0>, LC2<Scalar, Node1, Node2> > 
operator-(const OpTimesLC<Scalar, Node0>& x1, 
  const LC2<Scalar, Node1, Node2>& x2)
{
  return LC2<Scalar, OpTimesLC<Scalar, Node0>,
    LC2<Scalar, Node1, Node2> >(x1, x2, LCSubtract);
}


/* LC2 + OpTimesLC */
template <class Scalar, class Node1, class Node2, class Node3> inline
LC2<Scalar, LC2<Scalar, Node1, Node2>, OpTimesLC<Scalar, Node3> > 
operator+(const LC2<Scalar, Node1, Node2>& x1, 
  const OpTimesLC<Scalar, Node3>& x2)
{
  return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
    OpTimesLC<Scalar, Node3> >(x1, x2);
}

/* LC2 - OpTimesLC */
template <class Scalar, class Node1, class Node2, class Node3> inline
LC2<Scalar, LC2<Scalar, Node1, Node2>, OpTimesLC<Scalar, Node3> > 
operator-(const LC2<Scalar, Node1, Node2>& x1, 
  const OpTimesLC<Scalar, Node3>& x2)
{
  return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
    OpTimesLC<Scalar, Node3> >(x1, x2, LCSubtract);
}


/*======================================================================
 *
 *    add/subtract LC2, LC2
 *
 *======================================================================*/
  
/* LC2 + LC2 */
template <class Scalar, class Node1, class Node2, 
          class Node3, class Node4> inline
LC2<Scalar, LC2<Scalar, Node1, Node2>, LC2<Scalar, Node3, Node4> >
operator+(const LC2<Scalar, Node1, Node2>& x1, 
  const LC2<Scalar, Node3, Node4>& x2)
{
  return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
    LC2<Scalar, Node3, Node4> >(x1, x2);
}

/* LC2 - LC2 */
template <class Scalar, class Node1, class Node2, 
          class Node3, class Node4> inline
LC2<Scalar, LC2<Scalar, Node1, Node2>, LC2<Scalar, Node3, Node4> >
operator-(const LC2<Scalar, Node1, Node2>& x1, 
  const LC2<Scalar, Node3, Node4>& x2)
{
  return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
    LC2<Scalar, Node3, Node4> >(x1, x2, LCSubtract);
}


/*======================================================================
 *
 *    assignment of [OpTimesLC, LC2] to vector
 *
 *======================================================================*/
  
  
/* definition of assignment from 1-term linear combination to a vector */
template <class Scalar> 
template <class Node> inline
Vector<Scalar>& Vector<Scalar>::operator=(const PlayaExprTemplates::OpTimesLC<Scalar, Node>& x)
{
  if (this->ptr().get()==0)
  {
    *this = x.eval();
  }
  else if (x.containsVector(this->ptr().get()))
  {
    Vector<Scalar> rtn = x.eval();
    acceptCopyOf(rtn);
  }
  else
  {
    x.evalInto(*this);
  }
  return *this;
}

 
/* definition of assignment from N-term linear combination to a vector */
template <class Scalar>
template <class Node1, class Node2> inline
Vector<Scalar>& Vector<Scalar>::operator=(const PlayaExprTemplates::LC2<Scalar, Node1, Node2>& x)
{
  if (this->ptr().get()==0)
  {
    *this = x.eval();
  }
  else if (x.containsVector(this->ptr().get()))
  {
    Vector<Scalar> rtn = x.eval();
    acceptCopyOf(rtn);
  }
  else
  {
    x.evalInto(*this);
  }
  return *this;
}

 
/* definition of sum-assignment from 1-term linear combination to a vector */
template <class Scalar> 
template <class Node> inline
Vector<Scalar>& Vector<Scalar>::operator+=(const PlayaExprTemplates::OpTimesLC<Scalar, Node>& x)
{
  if (this->ptr().get()==0)
  {
    *this = x.eval();
  }
  else if (x.containsVector(this->ptr().get()))
  {
    Vector<Scalar> rtn = x.eval();
    update(1.0, rtn);
  }
  else
  {
    x.addInto(*this);
  }
  return *this;
}


/* definition of subtraction-assignment from 
 * 1-term linear combination to a vector */
template <class Scalar> 
template <class Node> inline
Vector<Scalar>& Vector<Scalar>::operator-=(const PlayaExprTemplates::OpTimesLC<Scalar, Node>& x)
{
  if (this->ptr().get()==0)
  {
    *this = -x.eval();
  }
  else if (x.containsVector(this->ptr().get()))
  {
    Vector<Scalar> rtn = x.eval();
    update(-1.0, rtn);
  }
  else
  {
    x.addInto(*this, LCSubtract);
  }
  return *this;
}


 
/* definition of sum-assignment from N-term linear combination to a vector */
template <class Scalar>
template <class Node1, class Node2> inline
Vector<Scalar>& Vector<Scalar>::operator+=(const PlayaExprTemplates::LC2<Scalar, Node1, Node2>& x)
{
  if (this->ptr().get()==0)
  {
    *this = x.eval();
  }
  else if (x.containsVector(this->ptr().get()))
  {
    Vector<Scalar> rtn = x.eval();
    update(1.0, rtn);
  }
  else
  {
    x.addInto(*this);
  }
  return *this;
}


 
/* definition of subtract-assignment from N-term linear combination 
 * to a vector */
template <class Scalar>
template <class Node1, class Node2> inline
Vector<Scalar>& Vector<Scalar>::operator-=(const PlayaExprTemplates::LC2<Scalar, Node1, Node2>& x)
{
  if (this->ptr().get()==0)
  {
    *this = -x.eval();
  }
  else if (x.containsVector(this->ptr().get()))
  {
    Vector<Scalar> rtn = x.eval();
    update(-1.0, rtn);
  }
  else
  {
    x.addInto(*this, LCSubtract);
  }
  return *this;
}




/*======================================================================
 *
 *    construction of vectors from [OpTimesLC, LC2]
 *
 *======================================================================*/
   

template <class Scalar>
template <class Node1, class Node2> inline
Vector<Scalar>::Vector(const PlayaExprTemplates::LC2<Scalar, Node1, Node2>& x)
  : Handle<VectorBase<Scalar> >(x.eval().ptr())
{;}

template <class Scalar> 
template <class Node> inline
Vector<Scalar>::Vector(const PlayaExprTemplates::OpTimesLC<Scalar, Node>& x)
  : Handle<VectorBase<Scalar> >(x.eval().ptr())
{;}




  

}



#endif
